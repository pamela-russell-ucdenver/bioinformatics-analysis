package translation;

import java.io.IOException;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;
import score.RegionScore;

/**
 * Differential translational efficiency
 * @author prussell
 *
 */
public class DifferentialTranslationalEfficiency implements RegionScore<Gene> {
	
	private TranslationalEfficiency te1;
	private TranslationalEfficiency te2;
	private double log2ratioCutoff;
	public static Logger logger = Logger.getLogger(DifferentialTranslationalEfficiency.class.getName());
	
	/**
	 * @param translationalEfficiency1 Translational efficiency object 1
	 * @param translationalEfficiency2 Translational efficiency object 2
	 * @param cutoffLog2ratio Cutoff for absolute value of log2 TE ratio
	 */
	public DifferentialTranslationalEfficiency(TranslationalEfficiency translationalEfficiency1, TranslationalEfficiency translationalEfficiency2, double cutoffLog2ratio) {
		te1 = translationalEfficiency1;
		te2 = translationalEfficiency2;
		log2ratioCutoff = cutoffLog2ratio;
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
		TranslationalEfficiency te1 = new TranslationalEfficiency(ribosomeBam1, controlBam1, geneBed, chrSizes, ribosomeGenomeTotal1, controlGenomeTotal1, ribosomeExonTotal1, controlExonTotal1, isStrandSpecific);
		TranslationalEfficiency te2 = new TranslationalEfficiency(ribosomeBam2, controlBam2, geneBed, chrSizes, ribosomeGenomeTotal2, controlGenomeTotal2, ribosomeExonTotal2, controlExonTotal2, isStrandSpecific);
		return new DifferentialTranslationalEfficiency(te1, te2, cutoffLog2ratio);
	}
	
	/**
	 * @return Translational efficiency object for control 1 and sample 1
	 */
	public TranslationalEfficiency getTE1() {
		return te1;
	}
	
	/**
	 * @return Translational efficiency object for control 2 and sample 2
	 */
	public TranslationalEfficiency getTE2() {
		return te2;
	}
	
	/**
	 * @return Sample 1 name
	 */
	public String getRibosome1Name() {
		return te1.getSampleName();
	}
	
	/**
	 * @return Sample 2 name
	 */
	public String getRibosome2Name() {
		return te2.getSampleName();
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
	public boolean isSignificant(double score) {
		return Math.abs(score) >= log2ratioCutoff;
	}

}
