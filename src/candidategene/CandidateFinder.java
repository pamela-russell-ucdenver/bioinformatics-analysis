package candidategene;

import guttmanlab.core.annotation.Annotation;

/**
 * Find candidate regions
 * @author prussell
 *
 * @param <T>
 */
public interface CandidateFinder<T extends Annotation> {
	
	/**
	 * @param region The region
	 * @return True iff the region is a candidate
	 */
	public boolean isCandidate(T region);
	
	/**
	 * Get line(s) to write to an output table
	 * Includes new line character
	 * @param region The region
	 * @return Data line for output table, or multiple lines, including new line character, or null
	 */
	public String getOutputTableLine(T region);
	
	/**
	 * Get line(s) to write to an output bed file of candidate regions
	 * Includes new line character
	 * @param region The region
	 * @return Data line for output bed file, or multiple lines, including new line character, or null
	 */
	public String getOutputBedLine(T region);
	
	/**
	 * @return A header line for the output table
	 */
	public String getOutputTableHeader();
	
}
