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
	 * @param significanceType Significance type
	 * @return True iff the score is significant
	 */
	public boolean isSignificant(double score, SignificanceType significanceType);
	
	/**
	 * Test whether the score for a region is significant
	 * @param region The region
	 * @param significanceType Significance type
	 * @return True iff the region's score is significant
	 */
	public boolean isSignificant(T region, SignificanceType significanceType);
	
	/**
	 * Get the name of the experiment
	 * @return Name of experiment
	 */
	public String getExperimentID();
	
	/**
	 * Instantiate this object using a line of information
	 * @param line Line describing how to instantiate
	 * @return The object of this class specified in the line
	 */
	public RegionScore<T> createFromConfigFileLine(String line);
	
	/**
	 * Get a description of a config file line to instantiate this class
	 * @return Config file format description
	 */
	public String getConfigFileLineFormat();
	
	
	/**
	 * Validate a config file line that is supposed to contain information to instantiate an object of this class
	 * Throw an exception if line is not valid
	 * @param line The line
	 */
	public void validateConfigFileLine(String line);
	
}
