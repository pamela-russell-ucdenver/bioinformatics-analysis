package score;

import java.util.Collection;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;

/**
 * A union of scores
 * A region is significant if it is significant in all of the scores
 * @author prussell
 *
 */
public class RegionScoreIntersection implements RegionScore<Gene> {
	
	private Collection<RegionScore<Gene>> scores;
	private String expID;
	private static Logger logger = Logger.getLogger(RegionScoreIntersection.class.getName());
		
	/**
	 * @param regionScores Score objects
	 * @param experimentID Name of this experiment/test
	 */
	public RegionScoreIntersection(Collection<RegionScore<Gene>> regionScores, String experimentID) {
		logger.info("");
		logger.info("Instantiating region score intersection object from " + regionScores.size() + " scores for experiment " + experimentID + ".");
		scores = regionScores;
		expID = experimentID;
		logger.info("");
		logger.info("Done instantiating score intersection for " + experimentID + ".");
	}
	
	@Override
	public double getScore(Gene region) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isSignificant(double score) {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getExperimentID() {
		return expID;
	}

	@Override
	public boolean isSignificant(Gene region) {
		for(RegionScore<Gene> score : scores) {
			if(!score.isSignificant(region)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getConfigFileLineFormat() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void validateConfigFileLine(String line) {
		throw new UnsupportedOperationException();
	}

}
