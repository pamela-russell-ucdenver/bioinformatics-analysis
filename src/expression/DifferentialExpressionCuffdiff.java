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

public class DifferentialExpressionCuffdiff extends AbstractRegionScore<Gene> implements DifferentialRegionScore<Gene> {
	
	private Map<String, CuffdiffRecord> recordsByID;
	public static double QVAL_CUTOFF = 0.05;
	private static Logger logger = Logger.getLogger(DifferentialExpressionCuffdiff.class.getName());
	
	/**
	 * @param cuffdiffOutputIsoformExpDiff Output file isoform_exp.diff from Cuffdiff
	 * @throws IOException
	 */
	public DifferentialExpressionCuffdiff(String cuffdiffOutputIsoformExpDiff) throws IOException {
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
	public boolean isSignificant(double score) {
		return score < QVAL_CUTOFF;
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
		return DifferentialExpressionCuffdiff.class.getSimpleName() + ":\tcuffdiff_isoform_exp_diff_file";
	}

	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 1) {
			throw new IllegalArgumentException("Invalid config file line: " + line + ". Format: " + getConfigFileLineFormat());
		}
	}

}
