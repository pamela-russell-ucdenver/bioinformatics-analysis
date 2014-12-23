package score;

import guttmanlab.core.annotation.Annotation;

/**
 * A region score that represents a change between two samples
 * @author prussell
 *
 * @param <T>
 */
public interface DifferentialRegionScore<T extends Annotation> extends RegionScore<T> {
	
	/**
	 * @return Experiment ID 1
	 */
	public String getExperimentID1();
	
	/**
	 * @return Experiment ID 2
	 */
	public String getExperimentID2();
	
	/**
	 * Returns whether sample 2 is "up" with respect to sample 1
	 * Regardless of whether the difference is significant
	 * @return True iff sample 2 is "up" with respect to sample 1
	 */
	public boolean experiment2IsUp(T region);
	
}
