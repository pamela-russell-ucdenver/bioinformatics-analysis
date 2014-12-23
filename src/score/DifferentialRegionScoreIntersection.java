package score;

import java.util.Collection;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;

/**
 * All scores must be significant in the same direction in order to be called significant
 * @author prussell
 *
 */
public class DifferentialRegionScoreIntersection implements DifferentialRegionScore<Gene> {

	private Collection<DifferentialRegionScore<Gene>> scores;
	private String expID;
	private static Logger logger = Logger.getLogger(DifferentialRegionScoreIntersection.class.getName());
	
	/**
	 * @param regionScores Score objects
	 * @param experimentID Name of this experiment/test
	 */
	public DifferentialRegionScoreIntersection(Collection<DifferentialRegionScore<Gene>> regionScores, String experimentID) {
		logger.info("");
		logger.info("Instantiating differential region score intersection object from " + regionScores.size() + " scores for experiment " + experimentID + ".");
		scores = regionScores;
		expID = experimentID;
		logger.info("");
		logger.info("Done instantiating differential score intersection for " + experimentID + ".");
	}
	
	@Override
	public double getScore(Gene region) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isSignificant(double score) {
		throw new UnsupportedOperationException();
	}

	/**
	 * True if all scores are significant in the same direction
	 */
	@Override
	public boolean isSignificant(Gene region) {
		boolean firstExpUp2 = scores.iterator().next().experiment2IsUp(region);
		for(DifferentialRegionScore<Gene> score : scores) {
			if(!score.isSignificant(region)) {
				return false;
			}
			if(score.experiment2IsUp(region) != firstExpUp2) {
				return false;
			}
		}
		return true;
	}

	@Override
	public String getExperimentID() {
		return expID;
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

	@Override
	public String getExperimentID1() {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getExperimentID2() {
		throw new UnsupportedOperationException();
	}

	/**
	 * True if experiment 2 is up in all scores
	 */
	@Override
	public boolean experiment2IsUp(Gene region) {
		boolean firstExpUp2 = scores.iterator().next().experiment2IsUp(region);
		for(DifferentialRegionScore<Gene> score : scores) {
			if(score.experiment2IsUp(region) != firstExpUp2) {
				throw new IllegalArgumentException("Not same direction of change for all scores");
			}
		}
		return firstExpUp2;
	}

}
