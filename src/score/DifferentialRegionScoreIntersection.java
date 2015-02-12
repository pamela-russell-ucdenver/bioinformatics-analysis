package score;

import java.util.Collection;

import org.apache.log4j.Level;
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
		//logger.setLevel(Level.DEBUG);
	}
	
	@Override
	public double getScore(Gene region) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isSignificant(double score, SignificanceType significanceType) {
		throw new UnsupportedOperationException();
	}

	/**
	 * True if all scores are significant in the same direction
	 */
	@Override
	public boolean isSignificant(Gene region, SignificanceType significanceType) {
		boolean firstExpUp2 = scores.iterator().next().experiment2IsUp(region);
		logger.debug("");
		logger.debug("Gene " + region.getName() + " exp " + scores.iterator().next().getExperimentID() + " is up: " + firstExpUp2);
		for(DifferentialRegionScore<Gene> score : scores) {
			if(!score.isSignificant(region, significanceType)) {
				logger.debug(score.getExperimentID() + "\tscore not significant\t" + score.getScore(region));
				return false;
			} else {
				logger.debug(score.getExperimentID() + "\tscore is significant\t" + score.getScore(region));
			}
			if(score.experiment2IsUp(region) != firstExpUp2) {
				logger.debug(score.getExperimentID() + "\tsignificant but wrong direction\t" + score.getScore(region));
				return false;
			} else {
				logger.debug(score.getExperimentID() + "\tsignificant in right direction\t" + score.getScore(region));
			}
		}
		logger.debug("Gene " + region.getName() + " is candidate for diff region score intersection");
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
