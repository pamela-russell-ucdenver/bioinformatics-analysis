package translation;

import java.io.IOException;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.util.StringParser;
import score.AbstractRegionScore;
import score.DifferentialRegionScore;
import score.GenericRegionScore;
import score.RegionScore;
import score.SignificanceType;

/**
 * Differential translational efficiency
 * @author prussell
 *
 */
public class DifferentialTranslationalEfficiency extends AbstractRegionScore<Gene> implements DifferentialRegionScore<Gene> {
	
	private TranslationalEfficiencyFromBam te1;
	private TranslationalEfficiencyFromBam te2;
	private double log2ratioCutoff;
	public static Logger logger = Logger.getLogger(DifferentialTranslationalEfficiency.class.getName());
	
	public DifferentialTranslationalEfficiency() {}
	
	/**
	 * @param translationalEfficiency1 Translational efficiency object 1
	 * @param translationalEfficiency2 Translational efficiency object 2
	 * @param cutoffLog2ratio Cutoff for absolute value of log2 TE ratio
	 */
	public DifferentialTranslationalEfficiency(TranslationalEfficiencyFromBam translationalEfficiency1, TranslationalEfficiencyFromBam translationalEfficiency2, double cutoffLog2ratio) {
		logger.info("");
		logger.info("Instantiating differential translational efficiency object...");
		te1 = translationalEfficiency1;
		te2 = translationalEfficiency2;
		log2ratioCutoff = cutoffLog2ratio;
		logger.info("");
		logger.info("Done instantiating differential translational efficiency object for experiments " + te1.getExperimentID() + " and " + te2.getExperimentID());
	}
	
	/**
	 * @param ribosomeBam1 Bam file of ribosome profiling sample 1
	 * @param ribosomeBam 2Bam file of ribosome profiling sample 2
	 * @param controlBam1 Bam file of control sample 1
	 * @param controlBam2 Bam file of control sample 2
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal1 Optional total number of ribosome reads mapped to genome (instead of computing from data), sample 1
	 * @param ribosomeGenomeTotal2 Optional total number of ribosome reads mapped to genome (instead of computing from data), sample 2
	 * @param controlGenomeTotal1 Optional total number of control reads mapped to genome (instead of computing from data), sample 1
	 * @param controlGenomeTotal2 Optional total number of control reads mapped to genome (instead of computing from data), sample 2
	 * @param ribosomeExonTotal1 Optional total number of ribosome reads mapped to exons (instead of computing from data), sample 1
	 * @param ribosomeExonTotal2 Optional total number of ribosome reads mapped to exons (instead of computing from data), sample 2
	 * @param controlExonTotal1 Optional total number of control reads mapped to exons (instead of computing from data), sample 1
	 * @param controlExonTotal2 Optional total number of control reads mapped to exons (instead of computing from data), sample 2
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @param cutoffLog2ratio
	 * @throws IOException
	 */
	public static DifferentialTranslationalEfficiency factory(String ribosomeBam1, String ribosomeBam2, String controlBam1, String controlBam2, String geneBed, String chrSizes, 
			double ribosomeGenomeTotal1, double ribosomeGenomeTotal2, double controlGenomeTotal1, double controlGenomeTotal2, double ribosomeExonTotal1, double ribosomeExonTotal2, 
			double controlExonTotal1, double controlExonTotal2, boolean isStrandSpecific, double cutoffLog2ratio) throws IOException {
		TranslationalEfficiencyFromBam te1 = new TranslationalEfficiencyFromBam(ribosomeBam1, controlBam1, geneBed, chrSizes, ribosomeGenomeTotal1, controlGenomeTotal1, ribosomeExonTotal1, controlExonTotal1, isStrandSpecific);
		TranslationalEfficiencyFromBam te2 = new TranslationalEfficiencyFromBam(ribosomeBam2, controlBam2, geneBed, chrSizes, ribosomeGenomeTotal2, controlGenomeTotal2, ribosomeExonTotal2, controlExonTotal2, isStrandSpecific);
		return new DifferentialTranslationalEfficiency(te1, te2, cutoffLog2ratio);
	}
	
	/**
	 * @return Translational efficiency object for control 1 and sample 1
	 */
	public TranslationalEfficiencyFromBam getTE1() {
		return te1;
	}
	
	/**
	 * @return Translational efficiency object for control 2 and sample 2
	 */
	public TranslationalEfficiencyFromBam getTE2() {
		return te2;
	}
	
	/**
	 * @return Sample 1 name
	 */
	public String getRibosome1Name() {
		return te1.getRibosomeName();
	}
	
	/**
	 * @return Sample 2 name
	 */
	public String getRibosome2Name() {
		return te2.getRibosomeName();
	}
	
	/**
	 * @return Control 1 name
	 */
	public String getControl1Name() {
		return te1.getControlName();
	}
	
	/**
	 * @return Control 2 name
	 */
	public String getControl2Name() {
		return te2.getControlName();
	}
	
	@Override
	public double getScore(Gene region) {
		double score1 = te1.getScore(region);
		double score2 = te2.getScore(region);
		double ratio = score2 / score1;
		double rtrn = Math.log(ratio) / Math.log(2);
		if(score1 == Double.NaN || score2 == Double.NaN) {
			logger.debug(region.getName() + "\tOne score is NaN\tlog_ratio=" + rtrn);
		}
		return rtrn;
	}

	@Override
	public boolean isSignificant(double score, SignificanceType significanceType) {
		switch(significanceType) {
		case EITHER_SAMPLE_UP:
			return Math.abs(score) >= log2ratioCutoff;
		case SAMPLE_1_UP:
			throw new IllegalArgumentException("Can't pass " + SignificanceType.SAMPLE_1_UP.toString() + " as significance type to method that just takes score");
		case SAMPLE_2_UP:
			throw new IllegalArgumentException("Can't pass " + SignificanceType.SAMPLE_2_UP.toString() + " as significance type to method that just takes score");
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case SINGLE_SAMPLE_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			return Math.abs(score) < log2ratioCutoff;
		default:
			throw new UnsupportedOperationException("Significance type " + significanceType.toString() + " not implemented.");
		}
	}

	@Override
	public boolean isSignificant(Gene region, SignificanceType significanceType) {
		double score = getScore(region);
		switch(significanceType) {
		case EITHER_SAMPLE_UP:
			return Math.abs(score) >= log2ratioCutoff;
		case SAMPLE_1_UP:
			return Math.abs(score) >= log2ratioCutoff && !experiment2IsUp(region);
		case SAMPLE_2_UP:
			return Math.abs(score) >= log2ratioCutoff && experiment2IsUp(region);
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case SINGLE_SAMPLE_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			return Math.abs(score) < log2ratioCutoff;
		default:
			throw new UnsupportedOperationException("Significance type " + significanceType.toString() + " not implemented.");
		}
	}
	
	@Override
	public String getExperimentID1() {
		return te1.getExperimentID();
	}

	@Override
	public String getExperimentID2() {
		return te2.getExperimentID();
	}

	@Override
	public boolean experiment2IsUp(Gene gene) {
		double score1 = te1.getScore(gene);
		if(score1 == Double.NaN) {
			throw new IllegalArgumentException("Score 1 is NaN");
		}
		double score2 = te2.getScore(gene);
		if(score2 == Double.NaN) {
			throw new IllegalArgumentException("Score 2 is NaN");
		}
		return score2 > score1;
	}

	@Override
	public String getExperimentID() {
		return "diff_TE_" + getExperimentID1() + "_" + getExperimentID2();
	}

	@Override
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		String rib1 = s.asString(0);
		String rib2 = s.asString(1);
		String con1 = s.asString(2);
		String con2 = s.asString(3);
		String geneBed = s.asString(4);
		String chrSizes = s.asString(5);
		double rgt1 = s.asDouble(6);
		double rgt2 = s.asDouble(7);
		double cgt1 = s.asDouble(8);
		double cgt2 = s.asDouble(9);
		double ret1 = s.asDouble(10);
		double ret2 = s.asDouble(11);
		double cet1 = s.asDouble(12);
		double cet2 = s.asDouble(13);
		boolean ss = s.asBoolean(14);
		double cutoff = s.asDouble(15);
		try {
			return factory(rib1, rib2, con1, con2, geneBed, chrSizes, rgt1, rgt2, cgt1, cgt2, ret1, ret2, cet1, cet2, ss, cutoff);
		} catch (IOException e) {
			logger.warn("Caught exception:");
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	@Override
	public String getConfigFileLineFormat() {
		String rtrn = DifferentialTranslationalEfficiency.class.getSimpleName() + ":\tribosomeBam1\tribosomeBam2\tcontrolBam1\tcontrolBam2\tgeneBed\tchrSizes\t";
		rtrn += "ribosomeGenomeTotal1\tribosomeGenomeTotal2\tcontrolGenomeTotal1\tcontrolGenomeTotal2\t";
		rtrn += "ribosomeExonTotal1\tribosomeExonTotal2\tcontrolExonTotal1\tcontrolExonTotal2\tisStrandSpecific\tcutoffLog2ratio";
		return rtrn;
	}

	@SuppressWarnings("unused")
	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 16) {
			logger.error("Field count is not 16: " + line);
			crashWithHelpMessage(line, logger);
		}
		try {
			String rib1 = s.asString(0);
			String rib2 = s.asString(1);
			String con1 = s.asString(2);
			String con2 = s.asString(3);
			String geneBed = s.asString(4);
			String chrSizes = s.asString(5);
			double rgt1 = s.asDouble(6);
			double rgt2 = s.asDouble(7);
			double cgt1 = s.asDouble(8);
			double cgt2 = s.asDouble(9);
			double ret1 = s.asDouble(10);
			double ret2 = s.asDouble(11);
			double cet1 = s.asDouble(12);
			double cet2 = s.asDouble(13);
			boolean ss = s.asBoolean(14);
			double cutoff = s.asDouble(15);
		} catch(Exception e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			crashWithHelpMessage(line, logger);
		}
	}


}
