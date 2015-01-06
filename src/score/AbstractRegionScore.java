package score;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Annotation;

public abstract class AbstractRegionScore<T extends Annotation> implements RegionScore<T> {
	
	public boolean isSignificant(T region, SignificanceType significanceType) {
		return isSignificant(getScore(region), significanceType);
	}
	
	/**
	 * When validating a config file line that supposedly contains information
	 * to instantiate the class, if the line is invalid, crash and print the
	 * correct format of the line
	 * @param line Actual provided config file line
	 * @param logger Logger to print message from
	 */
	public void crashWithHelpMessage(String line, Logger logger) {
		logger.error("Invalid config file line:");
		logger.error(line);
		logger.error("Format:");
		logger.error(getConfigFileLineFormat());
		System.exit(-1);
	}

}
