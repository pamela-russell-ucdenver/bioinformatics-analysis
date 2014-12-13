package expression;

import java.io.IOException;
import java.util.Map;

import guttmanlab.core.annotation.Gene;
import score.RegionScore;

public class DifferentialExpressionCuffdiff implements RegionScore<Gene> {
	
	private Map<String, CuffdiffRecord> recordsByID;
	public static double QVAL_CUTOFF = 0.05;
	
	/**
	 * @param cuffdiffOutputIsoformExpDiff Output file isoform_exp.diff from Cuffdiff
	 * @throws IOException
	 */
	public DifferentialExpressionCuffdiff(String cuffdiffOutputIsoformExpDiff) throws IOException {
		recordsByID = CuffdiffRecord.loadRecordsById(cuffdiffOutputIsoformExpDiff);
	}
	
	/**
	 * @param region Gene
	 * @return Cuffdiff output record for the gene
	 */
	public CuffdiffRecord getRecord(Gene region) {
		String name = region.getName();
		if(!recordsByID.containsKey(name)) {
			throw new IllegalArgumentException("Map of records does not contain key " + name + ".");
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

}
