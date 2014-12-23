package candidategene;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import net.sf.samtools.util.CloseableIterator;

import org.apache.log4j.Logger;

import expression.DifferentialExpressionCuffdiff;
import score.DifferentialRegionScore;
import score.DifferentialRegionScoreIntersection;
import score.DifferentialRegionScoreUnion;
import score.GenericDifferentialRegionScore;
import score.GenericRegionScore;
import score.RegionScore;
import score.RegionScoreFactory;
import score.RegionScoreIntersection;
import score.ScoreType;
import score.SignificanceType;
import translation.DifferentialTranslationalEfficiency;
import translation.TranslationalEfficiency;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;
import guttmanlab.core.util.StringParser;

/**
 * Find genes that are significant according to multiple scores
 * @author prussell
 *
 */
public class CandidateFinderCombinedScores implements CandidateFinder<Gene> {
	
	/**
	 * Regular scores, not comparisons of two samples
	 */
	private Map<RegionScore<Gene>, SignificanceType> scores;
	
	/**
	 * Scores that are comparisons of two samples
	 */
	private Map<DifferentialRegionScore<Gene>, SignificanceType> diffScores;
	
	private static Logger logger = Logger.getLogger(CandidateFinderCombinedScores.class.getName());
	
	private CandidateFinderCombinedScores(String configFile) throws IOException {
		logger.info("");
		logger.info("Instantiating combined score candidate finder with config file " + configFile + "...");
		initializeScores(configFile);
		logger.info("");
		logger.info("Done instantiating combined score candidate finder.");
	}
	
	/**
	 * Get the scores from the config file
	 * @param configFile Config file
	 * @throws IOException
	 */
	private void initializeScores(String configFile) throws IOException {
		scores = new HashMap<RegionScore<Gene>, SignificanceType>();
		diffScores = new HashMap<DifferentialRegionScore<Gene>, SignificanceType>();
		BufferedReader b = new BufferedReader(new FileReader(configFile));
		// Map to keep track of scores that will comprise an intersection of multiple scores
		Map<String, Map<SignificanceType, Collection<RegionScore<Gene>>>> intersectionScores = new HashMap<String, Map<SignificanceType, Collection<RegionScore<Gene>>>>();
		// Map to keep track of scores that will comprise a union of multiple scores
		Map<String, Map<SignificanceType, Collection<RegionScore<Gene>>>> unionScores = new HashMap<String, Map<SignificanceType, Collection<RegionScore<Gene>>>>();
		// Map to keep track of comparison scores that will comprise an intersection of multiple scores
		Map<String, Map<SignificanceType, Collection<DifferentialRegionScore<Gene>>>> intersectionDiffScores = new HashMap<String, Map<SignificanceType, Collection<DifferentialRegionScore<Gene>>>>();
		// Map to keep track of comparison scores that will comprise a union of multiple scores
		Map<String, Map<SignificanceType, Collection<DifferentialRegionScore<Gene>>>> unionDiffScores = new HashMap<String, Map<SignificanceType, Collection<DifferentialRegionScore<Gene>>>>();
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			try {
				ScoreType scoreType = ScoreType.fromString(s.asString(0));
				SignificanceType sigType = SignificanceType.fromString(s.asString(1));
				if(scoreType.equals(ScoreType.INTERSECTION_DIFFERENTIAL)) {
					String expID = s.asString(2);
					@SuppressWarnings("unchecked")
					DifferentialRegionScore<Gene> score = (DifferentialRegionScore<Gene>) RegionScoreFactory.createDiffScoreFromConfigFileLine(line, 3);
					if(!intersectionDiffScores.containsKey(expID)) {
						intersectionDiffScores.put(expID, new HashMap<SignificanceType, Collection<DifferentialRegionScore<Gene>>>());
						intersectionDiffScores.get(expID).put(sigType, new ArrayList<DifferentialRegionScore<Gene>>());
					}
					intersectionDiffScores.get(expID).get(sigType).add(score);
					continue;
				}
				if(scoreType.equals(ScoreType.INTERSECTION_REGULAR)) {
					String expID = s.asString(2);
					@SuppressWarnings("unchecked")
					RegionScore<Gene> score = (RegionScore<Gene>) RegionScoreFactory.createScoreFromConfigFileLine(line, 3);
					if(!intersectionScores.containsKey(expID)) {
						intersectionScores.put(expID, new HashMap<SignificanceType, Collection<RegionScore<Gene>>>());
						intersectionScores.get(expID).put(sigType, new ArrayList<RegionScore<Gene>>());
					}
					intersectionScores.get(expID).get(sigType).add(score);
					continue;
				}
				if(scoreType.equals(ScoreType.SINGLE_DIFFERENTIAL)) {
					@SuppressWarnings("unchecked")
					DifferentialRegionScore<Gene> score = (DifferentialRegionScore<Gene>) RegionScoreFactory.createDiffScoreFromConfigFileLine(line, 2);
					diffScores.put(score, sigType); // Store the score
					continue;
				}
				if(scoreType.equals(ScoreType.SINGLE_REGULAR)) {
					@SuppressWarnings("unchecked")
					RegionScore<Gene> score = (RegionScore<Gene>) RegionScoreFactory.createScoreFromConfigFileLine(line, 2);
					scores.put(score, sigType); // Store the score
					continue;
				}
				if(scoreType.equals(ScoreType.UNION_DIFFERENTIAL)) {
					String expID = s.asString(2);
					@SuppressWarnings("unchecked")
					DifferentialRegionScore<Gene> score = (DifferentialRegionScore<Gene>) RegionScoreFactory.createDiffScoreFromConfigFileLine(line, 3);
					if(!unionDiffScores.containsKey(expID)) {
						unionDiffScores.put(expID, new HashMap<SignificanceType, Collection<DifferentialRegionScore<Gene>>>());
						unionDiffScores.get(expID).put(sigType, new ArrayList<DifferentialRegionScore<Gene>>());
					}
					unionDiffScores.get(expID).get(sigType).add(score);
					continue;
				}
				if(scoreType.equals(ScoreType.UNION_REGULAR)) {
					String expID = s.asString(2);
					@SuppressWarnings("unchecked")
					RegionScore<Gene> score = (RegionScore<Gene>) RegionScoreFactory.createScoreFromConfigFileLine(line, 3);
					if(!unionScores.containsKey(expID)) {
						unionScores.put(expID, new HashMap<SignificanceType, Collection<RegionScore<Gene>>>());
						unionScores.get(expID).put(sigType, new ArrayList<RegionScore<Gene>>());
					}
					unionScores.get(expID).get(sigType).add(score);
					continue;
				}
			} catch(Exception e) {
				crashWithHelpMessage(line);
			}
		}
		b.close();
		
		// Store intersection and union scores after instantiating them
		
		for(String id : intersectionScores.keySet()) {
			for(SignificanceType sigType : intersectionScores.get(id).keySet()) {
				scores.put(new RegionScoreIntersection(intersectionScores.get(id).get(sigType), id), sigType);
			}
		}
		
		for(String id : unionScores.keySet()) {
			for(SignificanceType sigType : unionScores.get(id).keySet()) {
				scores.put(new RegionScoreIntersection(unionScores.get(id).get(sigType), id), sigType);
			}
		}
		
		for(String id : intersectionDiffScores.keySet()) {
			for(SignificanceType sigType : intersectionDiffScores.get(id).keySet()) {
				diffScores.put(new DifferentialRegionScoreIntersection(intersectionDiffScores.get(id).get(sigType), id), sigType);
			}
		}
		
		for(String id : unionDiffScores.keySet()) {
			for(SignificanceType sigType : unionDiffScores.get(id).keySet()) {
				diffScores.put(new DifferentialRegionScoreUnion(unionDiffScores.get(id).get(sigType), id), sigType);
			}
		}
		

		
	}
	
	private static void crashWithHelpMessage(String line) {
		logger.error("");
		logger.error("Invalid config file line:");
		logger.error(line);
		logger.error("");
		printConfigFileDescription();
		System.exit(-1);
	}
	
	private static void printConfigFileDescription() {
		logger.info("\nline format:");
		logger.info("score_type\tscore_group_id_for_intersection_or_union[omit_if_not_combined]\tsignificance_type\tscore_name\tscore_info\n");
		logger.info("Score types:");
		logger.info(ScoreType.commaSeparatedList() + "\n");
		logger.info("Significance types:");
		logger.info(SignificanceType.commaSeparatedList() + "\n");
		logger.info("Score info formats:");
		logger.info(new GenericRegionScore().getConfigFileLineFormat());
		logger.info(new GenericDifferentialRegionScore().getConfigFileLineFormat());
		logger.info(new DifferentialExpressionCuffdiff().getConfigFileLineFormat());
		logger.info(new TranslationalEfficiency().getConfigFileLineFormat());
		logger.info(new DifferentialTranslationalEfficiency().getConfigFileLineFormat());
	}
	
	@Override
	public boolean isCandidate(Gene region) {
		for(RegionScore<Gene> score : scores.keySet()) {
			if(!score.isSignificant(region)) {
				return false;
			}
		}
		for(DifferentialRegionScore<Gene> diffScore : diffScores.keySet()) {
			SignificanceType sigType = diffScores.get(diffScore);
			if(!diffScore.isSignificant(region)) {
				return false;
			}
			if(sigType.equals(SignificanceType.SAMPLE_1_UP) && diffScore.experiment2IsUp(region)) {
				return false;
			}
			if(sigType.equals(SignificanceType.SAMPLE_2_UP) && !diffScore.experiment2IsUp(region)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public String getOutputTableLine(Gene region) {
		String rtrn = region.getName() + "\t";
		rtrn += region.toUCSC() + "\t";
		rtrn += isCandidate(region) + "\t";
		for(RegionScore<Gene> score : scores.keySet()) {
			rtrn += score.getScore(region) + "\t";
			rtrn += score.isSignificant(region) + "\t";
		}
		for(DifferentialRegionScore<Gene> score : diffScores.keySet()) {
			rtrn += score.getScore(region) + "\t";
			rtrn += score.experiment2IsUp(region) + "\t";
			rtrn += score.isSignificant(region) + "\t";
		}
		rtrn += "\n";
		return rtrn;
	}

	@Override
	public String getOutputBedLine(Gene region) {
		return region.toBED();
	}

	@Override
	public String getOutputTableHeader() {
		String rtrn = "gene_ID\t";
		rtrn += "coordinates\t";
		rtrn += "is_candidate\t";
		for(RegionScore<Gene> score : scores.keySet()) {
			rtrn += "score_" + score.getExperimentID() + "\t";
			rtrn += "is_significant_" + score.getExperimentID() + "\t";
		}
		for(DifferentialRegionScore<Gene> score : diffScores.keySet()) {
			rtrn += "score_" + score.getExperimentID() + "\t";
			rtrn += "experiment_2_is_up_" + score.getExperimentID() + "\t";
			rtrn += "is_significant_" + score.getExperimentID() + "\t";
		}
		return rtrn;
	}
	
	
	/**
	 * Write all genes to a table and candidates to a bed file
	 * @param outFilePrefix Output file prefix
	 * @throws IOException
	 */
	private void writeResults(String geneBed, String chrSizes, String outFilePrefix) throws IOException {
		Map<String, FeatureCollection<Gene>> genes = BEDFileIO.loadFromFileByReferenceName(geneBed, chrSizes);
		String outTable = outFilePrefix + ".out";
		String outBed = outFilePrefix + ".candidates.bed";
		logger.info("");
		logger.info("Writing candidate uORFs to table " + outTable + " and bed file " + outBed + "...");
		FileWriter wt = new FileWriter(outTable);
		FileWriter wb = new FileWriter(outBed);
		wt.write(getOutputTableHeader() + "\n");
		for(String chr : genes.keySet()) {
			logger.info(chr);
			int numGenes = genes.get(chr).getNumAnnotations();
			CountLogger cl = new CountLogger(numGenes, 10);
			CloseableIterator<Gene> iter = genes.get(chr).sortedIterator();
			while(iter.hasNext()) {
				cl.advance();
				Gene gene = iter.next();
				String line = getOutputTableLine(gene);
				wt.write(line);
				//wt.flush();
				if(isCandidate(gene)) wb.write(getOutputBedLine(gene));
				//wb.flush();
			}
			iter.close();
		}
		wt.close();
		wb.close();
		logger.info("Done writing file.");
	}


	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-cf", "Config file", true);
		p.addStringArg("-gb", "Bed file of genes to test for candidates", true);
		p.addStringArg("-cs", "Chromosome size file", true);
		p.addStringArg("-o", "Output file prefix", true);
		
		try {
			p.parse(args);
		} catch(Exception e) {
			System.out.println();
			p.printHelpMessage();
			System.out.println();
			printConfigFileDescription();
		}
		
		String configFile = p.getStringArg("-cf");
		String geneBed = p.getStringArg("-gb");
		String chrSizes = p.getStringArg("-cs");
		String outFile = p.getStringArg("-o");
		
		CandidateFinderCombinedScores cf = new CandidateFinderCombinedScores(configFile);
		
		cf.writeResults(geneBed, chrSizes, outFile);
		
	}

}
