package expression;

import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.util.StringParser;
import score.AbstractRegionScore;
import score.DifferentialRegionScore;
import score.GenericRegionScore;
import score.RegionScore;
import score.SignificanceType;

public class DifferentialExpressionCuffdiff extends AbstractRegionScore<Gene> implements DifferentialRegionScore<Gene> {
	
	private Map<String, CuffdiffRecord> recordsByID;
	public static double QVAL_CUTOFF = 0.05;
	private static Logger logger = Logger.getLogger(DifferentialExpressionCuffdiff.class.getName());
	
	/**
	 * @param cuffdiffOutputIsoformExpDiff Output file isoform_exp.diff from Cuffdiff
	 * @throws IOException
	 */
	public DifferentialExpressionCuffdiff(String configFileLine) throws IOException {
		StringParser s = new StringParser();
		s.parse(configFileLine);
		String cuffdiffOutputIsoformExpDiff = s.asString(0);
		QVAL_CUTOFF = s.asDouble(1);
		logger.info("");
		logger.info("Instantiating differential expression cuffdiff object with file " + cuffdiffOutputIsoformExpDiff + "...");
		recordsByID = CuffdiffRecord.loadRecordsById(cuffdiffOutputIsoformExpDiff);
		logger.info("");
		logger.info("Done instantiating differential expression cuffdiff object. Experiment IDs are " + getExperimentID1() + " and " + getExperimentID2() + ".");
	}
	
	public DifferentialExpressionCuffdiff() {}

	/**
	 * @param region Gene
	 * @return Cuffdiff output record for the gene or null if name doesn't exist
	 */
	public CuffdiffRecord getRecord(Gene region) {
		String name = region.getName();
		if(!recordsByID.containsKey(name)) {
			return null;
		}
		return recordsByID.get(name);
	}
	
	@Override
	public double getScore(Gene region) {
		return getRecord(region).getQval();
	}

	@Override
	public boolean isSignificant(double score, SignificanceType significanceType) {
		switch(significanceType) {
		case EITHER_SAMPLE_UP:
			return score < QVAL_CUTOFF;
		case SAMPLE_1_UP:
			throw new IllegalArgumentException("Can't pass " + SignificanceType.SAMPLE_1_UP.toString() + " as significance type to method that just takes score");
		case SAMPLE_2_UP:
			throw new IllegalArgumentException("Can't pass " + SignificanceType.SAMPLE_2_UP.toString() + " as significance type to method that just takes score");
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case SINGLE_SAMPLE_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			return score >= QVAL_CUTOFF;
		default:
			throw new UnsupportedOperationException("Significance type " + significanceType.toString() + " not implemented.");
		}
	}
	
	@Override
	public boolean isSignificant(Gene region, SignificanceType significanceType) {
		double score = getScore(region);
		switch(significanceType) {
		case EITHER_SAMPLE_UP:
			return score < QVAL_CUTOFF;
		case SAMPLE_1_UP:
			return score < QVAL_CUTOFF && !experiment2IsUp(region);
		case SAMPLE_2_UP:
			return score < QVAL_CUTOFF && experiment2IsUp(region);
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case SINGLE_SAMPLE_SIGNIFICANT:
			throw new IllegalArgumentException("Can't use single sample significance type for differential expression");
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			return score >= QVAL_CUTOFF;
		default:
			throw new UnsupportedOperationException("Significance type " + significanceType.toString() + " not implemented.");
		}
	}

	@Override
	public String getExperimentID1() {
		return recordsByID.values().iterator().next().getSample1();
	}

	@Override
	public String getExperimentID2() {
		return recordsByID.values().iterator().next().getSample2();
	}

	@Override
	public boolean experiment2IsUp(Gene gene) {
		CuffdiffRecord record = recordsByID.get(gene.getName());
		return record.getLog2fpkmRatio() > 0;
	}

	@Override
	public String getExperimentID() {
		return "diff_exp_cuffdiff_" + getExperimentID1() + "_" + getExperimentID2();
	}

	@Override
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		validateConfigFileLine(line);
		try {
			return new DifferentialExpressionCuffdiff(line);
		} catch(IOException e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	@Override
	public String getConfigFileLineFormat() {
		return DifferentialExpressionCuffdiff.class.getSimpleName() + ":\tcuffdiff_isoform_exp_diff_file\tmax_qval";
	}

	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 2) {
			throw new IllegalArgumentException("Invalid config file line: " + line + ". Format: " + getConfigFileLineFormat());
		}
	}

}
