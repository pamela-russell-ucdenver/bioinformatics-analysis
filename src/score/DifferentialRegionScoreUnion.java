package score;

import java.util.Collection;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;

/**
 * A union of differential scores
 * One score must be significant in order to be called significant
 * @author prussell
 *
 */
public class DifferentialRegionScoreUnion implements DifferentialRegionScore<Gene> {

	private Collection<DifferentialRegionScore<Gene>> scores;
	private String expID;
	private static Logger logger = Logger.getLogger(DifferentialRegionScoreUnion.class.getName());
	
	/**
	 * @param regionScores Score objects
	 * @param experimentID Name of this experiment/test
	 */
	public DifferentialRegionScoreUnion(Collection<DifferentialRegionScore<Gene>> regionScores, String experimentID) {
		logger.info("");
		logger.info("Instantiating differential region score union object from " + regionScores.size() + " scores for experiment " + experimentID + ".");
		scores = regionScores;
		expID = experimentID;
		logger.info("");
		logger.info("Done instantiating differential score union for " + experimentID + ".");
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
	 * True if at least one score is significant
	 */
	@Override
	public boolean isSignificant(Gene region, SignificanceType significanceType) {
		for(DifferentialRegionScore<Gene> score : scores) {
			if(score.isSignificant(region, significanceType)) {
				return true;
			}
		}
		return false;
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
	 * True if experiment 2 is up in at least one score
	 */
	@Override
	public boolean experiment2IsUp(Gene region) {
		for(DifferentialRegionScore<Gene> score : scores) {
			try {
				if(score.experiment2IsUp(region)) {
					return true;
				}
			} catch(Exception e) {
				logger.warn(score.getClass().getSimpleName() + " can't asses which sample is up for gene " + region.getName() + ". Throwing exception.");
				throw(e);
			}
		}
		return false;
	}

}
