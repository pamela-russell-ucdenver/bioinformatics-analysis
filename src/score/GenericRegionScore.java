package score;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.util.StringParser;

/**
 * A score that is read in from a table
 * Could be any type of score
 * @author prussell
 *
 */
public class GenericRegionScore extends AbstractRegionScore<Gene> {
	
	private Map<String, Double> scoresByGeneName;
	private double cutoff;
	private boolean scoreCutoffIsMax;
	private String expID;
	private static Logger logger = Logger.getLogger(GenericRegionScore.class.getName());
	
	public GenericRegionScore() {}
	
	/**
	 * @param scoreTable Score table where each line has gene name in column 0 and score in specified column (zero based)
	 * @param geneNameColNum Zero based column number containing the gene name
	 * @param scoreColNum Zero based column number containing the score (gene name must be in column 0)
	 * @param scoreCutoff Score significance cutoff
	 * @param cutoffIsMax True if the cutoff is a maximum, false otherwise
	 * @param experimentID Experiment ID
	 * @throws IOException
	 */
	public GenericRegionScore(String scoreTable, int geneNameColNum, int scoreColNum, double scoreCutoff, boolean cutoffIsMax, String experimentID) throws IOException {
		logger.info("");
		logger.info("Instantiating generic region score from column " + scoreColNum + " of table " + scoreTable + ". Cutoff is " + scoreCutoff + ".");
		cutoff = scoreCutoff;
		scoreCutoffIsMax = cutoffIsMax;
		expID = experimentID;
		initializeScores(scoreTable, geneNameColNum, scoreColNum);
		logger.info("");
		logger.info("Done instantiating generic region score.");
	}
	
	protected void initializeScores(String tableFile, int geneNameColNum, int scoreColNum) throws IOException {
		scoresByGeneName = new HashMap<String, Double>();
		BufferedReader b = new BufferedReader(new FileReader(tableFile));
		StringParser s = new StringParser();
		boolean firstLine = true; // Catch exceptions if there is a header line
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			try {
				String name = s.asString(geneNameColNum);
				if(scoresByGeneName.containsKey(name)) {
					logger.warn("SKIPPING LINE: Score map already contains key " + name + ".");
					continue;
				}
				scoresByGeneName.put(name, Double.valueOf(s.asDouble(scoreColNum)));
				firstLine = false;
			} catch(Exception e) {
				if(firstLine) {
					continue;
				}
				b.close();
				logger.error("Exception on line: " + line);
				throw e;
			}
		}
		b.close();
	}
		
	@Override
	public double getScore(Gene region) {
		String name = region.getName();
		Double rtrn = scoresByGeneName.get(name);
		if(rtrn == null) {
			throw new IllegalArgumentException("Score map does not contain key " + name + ".");
		}
		return rtrn.doubleValue();
	}

	@Override
	public boolean isSignificant(double score, SignificanceType significanceType) {
		switch(significanceType) {
		case EITHER_SAMPLE_UP:
			throw new IllegalArgumentException("Can't use two-sample significance type for single-sample generic region score");
		case SAMPLE_1_UP:
			throw new IllegalArgumentException("Can't use two-sample significance type for single-sample generic region score");
		case SAMPLE_2_UP:
			throw new IllegalArgumentException("Can't use two-sample significance type for single-sample generic region score");
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			if(scoreCutoffIsMax) return score > cutoff;
			return score <= cutoff;
		case SINGLE_SAMPLE_SIGNIFICANT:
			if(scoreCutoffIsMax) return score <= cutoff;
			return score >= cutoff;
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use two-sample significance type for single-sample generic region score");
		default:
			throw new UnsupportedOperationException("Significance type " + significanceType.toString() + " not implemented.");
		}
	}

	@Override
	public String getExperimentID() {
		return expID;
	}

	@Override
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		validateConfigFileLine(line);
		StringParser s = new StringParser();
		String table = s.asString(0);
		int geneCol = s.asInt(1);
		int col = s.asInt(2);
		double cutoff = s.asDouble(3);
		boolean max = s.asBoolean(4);
		String id = s.asString(5);
		try {
			return new GenericRegionScore(table, geneCol, col, cutoff, max, id);
		} catch (IOException e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	@Override
	public String getConfigFileLineFormat() {
		return GenericRegionScore.class.getSimpleName() + ":\tscore_table\tgene_name_col_num\tscore_col_num\tscore_cutoff\tcutoff_is_max\texperiment_ID";
	}
	
	@SuppressWarnings("unused")
	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 6) {
			logger.error("Field count is not 6: " + line);
			crashWithHelpMessage(line, logger);
		}
		try {
			String table = s.asString(0);
			int geneCol = s.asInt(1);
			int col = s.asInt(2);
			double cutoff = s.asDouble(3);
			boolean max = s.asBoolean(4);
			String id = s.asString(5);
		} catch(Exception e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			crashWithHelpMessage(line, logger);
		}
	}


}
