package score;

import guttmanlab.core.annotation.Annotation;

/**
 * A score of a region
 * @author prussell
 *
 * @param <T> Generic annotation type
 */
public interface RegionScore<T extends Annotation> {
	
	/**
	 * Get the score over a region
	 * @param region The region
	 * @return The score
	 */
	public double getScore(T region);
	
	/**
	 * Test whether a score is significant
	 * @param score The score
	 * @return True iff the score is significant
	 */
	public boolean isSignificant(double score);
	
}
