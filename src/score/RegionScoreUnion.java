package score;

import java.util.Collection;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;

/**
 * A union of scores
 * A region is significant if it is significant in one of the scores
 * @author prussell
 *
 */
public class RegionScoreUnion implements RegionScore<Gene> {
	
	private Collection<RegionScore<Gene>> scores;
	private String expID;
	private static Logger logger = Logger.getLogger(RegionScoreUnion.class.getName());
		
	/**
	 * @param regionScores Score objects
	 * @param experimentID Name of this experiment/test
	 */
	public RegionScoreUnion(Collection<RegionScore<Gene>> regionScores, String experimentID) {
		logger.info("");
		logger.info("Instantiating region score union object from " + regionScores.size() + " scores for experiment " + experimentID + ".");
		scores = regionScores;
		expID = experimentID;
		logger.info("");
		logger.info("Done instantiating score union for " + experimentID + ".");
	}
	
	@Override
	public double getScore(Gene region) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isSignificant(double score, SignificanceType significanceType) {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getExperimentID() {
		return expID;
	}

	@Override
	public boolean isSignificant(Gene region, SignificanceType significanceType) {
		for(RegionScore<Gene> score : scores) {
			if(score.isSignificant(region, significanceType)) {
				return true;
			}
		}
		return false;
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
