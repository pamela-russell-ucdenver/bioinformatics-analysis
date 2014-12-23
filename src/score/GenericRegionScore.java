package score;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.util.StringParser;

/**
 * A score that is read in from a table
 * Could be any type of score
 * @author prussell
 *
 */
public class GenericRegionScore extends AbstractRegionScore<Gene> {
	
	private Map<String, Double> scoresByGeneName;
	private double cutoff;
	private boolean scoreCutoffIsMax;
	private String expID;
	private static Logger logger = Logger.getLogger(GenericRegionScore.class.getName());
	
	public GenericRegionScore() {}
	
	/**
	 * @param scoreTable Score table where each line has gene name in column 0 and score in specified column (zero based)
	 * @param scoreColNum Zero based column number containing the score (gene name must be in column 0)
	 * @param scoreCutoff Score significance cutoff
	 * @param cutoffIsMax True if the cutoff is a maximum, false otherwise
	 * @param experimentID Experiment ID
	 * @throws IOException
	 */
	public GenericRegionScore(String scoreTable, int scoreColNum, double scoreCutoff, boolean cutoffIsMax, String experimentID) throws IOException {
		logger.info("");
		logger.info("Instantiating generic region score from column " + scoreColNum + " of table " + scoreTable + ". Cutoff is " + scoreCutoff + ".");
		cutoff = scoreCutoff;
		scoreCutoffIsMax = cutoffIsMax;
		expID = experimentID;
		initializeScores(scoreTable, scoreColNum);
		logger.info("");
		logger.info("Done instantiating generic region score.");
	}
	
	protected void initializeScores(String tableFile, int scoreColNum) throws IOException {
		scoresByGeneName = new HashMap<String, Double>();
		BufferedReader b = new BufferedReader(new FileReader(tableFile));
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			try {
				String name = s.asString(0);
				if(scoresByGeneName.containsKey(name)) {
					throw new IllegalStateException("Score map already contains key " + name + ".");
				}
				scoresByGeneName.put(name, Double.valueOf(s.asDouble(scoreColNum)));
			} catch(Exception e) {
				b.close();
				logger.error("Exception on line: " + line);
				throw e;
			}
		}
		b.close();
	}
		
	@Override
	public double getScore(Gene region) {
		String name = region.getName();
		Double rtrn = scoresByGeneName.get(name);
		if(rtrn == null) {
			throw new IllegalArgumentException("Score map does not contain key " + name + ".");
		}
		return rtrn.doubleValue();
	}

	@Override
	public boolean isSignificant(double score) {
		if(scoreCutoffIsMax) return score <= cutoff;
		return score >= cutoff;
	}

	@Override
	public String getExperimentID() {
		return expID;
	}

	@Override
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		validateConfigFileLine(line);
		StringParser s = new StringParser();
		String table = s.asString(0);
		int col = s.asInt(1);
		double cutoff = s.asDouble(2);
		boolean max = s.asBoolean(3);
		String id = s.asString(4);
		try {
			return new GenericRegionScore(table, col, cutoff, max, id);
		} catch (IOException e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	@Override
	public String getConfigFileLineFormat() {
		return "score_table[string]\tscore_col_num[int]\tscore_cutoff[double]\tcutoff_is_max[boolean]\texperiment_ID[string]";
	}
	
	@SuppressWarnings("unused")
	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 5) {
			crashWithHelpMessage(line, logger);
		}
		try {
			String table = s.asString(0);
			int col = s.asInt(1);
			double cutoff = s.asDouble(2);
			boolean max = s.asBoolean(3);
			String id = s.asString(4);
		} catch(Exception e) {
			crashWithHelpMessage(line, logger);
		}
	}


}
