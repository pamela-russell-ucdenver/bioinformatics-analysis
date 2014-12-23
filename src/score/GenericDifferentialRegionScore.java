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
 * A differential score that is read in from a table
 * E.g. differential expression
 * The table must contain numerical scores for 2 samples and a differential score, e.g. a P-value for differential expression
 * @author prussell
 *
 */
public class GenericDifferentialRegionScore extends GenericRegionScore implements DifferentialRegionScore<Gene> {
	
	private Map<String, Boolean> experiment2IsUp;
	private static Logger logger = Logger.getLogger(GenericDifferentialRegionScore.class.getName());
	private String experiment1ID;
	private String experiment2ID;
	
	public GenericDifferentialRegionScore() {}
	
	/**
	 * @param scoreTable Table containing gene name in column 0 and differential score (e.g. change P value), experiment 1 score, and experiment 2 score in specified columns
	 * @param scoreCutoff Significance cutoff for differential score
	 * @param cutoffIsMax True if the cutoff is a maximum, false otherwise
	 * @param experimentID Experiment ID
	 * @param differentialScoreColNum Zero based column number containing the differential score (e.g. change P value) (gene name must be in column 0)
	 * @param exp1scoreColNum Zero based column number containing the score for experiment 1
	 * @param exp2scoreColNum Zero based column number containing the score for experiment 2
	 * @throws IOException
	 */
	public GenericDifferentialRegionScore(String scoreTable, double scoreCutoff, boolean cutoffIsMax, String experimentID, int differentialScoreColNum, String exp1ID, String exp2ID, int exp1scoreColNum, int exp2scoreColNum) throws IOException {
		super(scoreTable, differentialScoreColNum, scoreCutoff, cutoffIsMax, experimentID);
		initializeComparisons(scoreTable, exp1scoreColNum, exp2scoreColNum);
		experiment1ID = exp1ID;
		experiment2ID = exp2ID;
	}

	/**
	 * Store the comparisons between the two samples, i.e. which sample is up
	 * @param tableFile Score table
	 * @param exp1scoreColNum Column number containing score for experiment 1
	 * @param exp2scoreColNum Column number containing score for experiment 2
	 * @throws IOException
	 */
	protected void initializeComparisons(String tableFile, int exp1scoreColNum, int exp2scoreColNum) throws IOException {
		experiment2IsUp = new HashMap<String, Boolean>();
		BufferedReader b = new BufferedReader(new FileReader(tableFile));
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			try {
				String name = s.asString(0);
				if(experiment2IsUp.containsKey(name)) {
					throw new IllegalStateException("Comparison map already contains key " + name + ".");
				}
				double score1 = s.asDouble(exp1scoreColNum);
				double score2 = s.asDouble(exp2scoreColNum);
				experiment2IsUp.put(name, Boolean.valueOf(score2 > score1));
			} catch(Exception e) {
				b.close();
				logger.error("Exception on line: " + line);
				throw e;
			}
		}
		b.close();
	}

	
	@Override
	public String getExperimentID1() {
		return experiment1ID;
	}

	@Override
	public String getExperimentID2() {
		return experiment2ID;
	}

	@Override
	public boolean experiment2IsUp(Gene region) {
		String name = region.getName();
		Boolean rtrn = experiment2IsUp.get(name);
		if(rtrn == null) {
			throw new IllegalArgumentException("Comparison map does not contain key " + name + ".");
		}
		return rtrn.booleanValue();
	}
	
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		validateConfigFileLine(line);
		StringParser s = new StringParser();
		s.parse(line);
		String st = s.asString(0);
		double cutoff = s.asDouble(1);
		boolean max = s.asBoolean(2);
		String id = s.asString(3);
		int dc = s.asInt(4);
		String id1 = s.asString(5);
		String id2 = s.asString(6);
		int c1 = s.asInt(7);
		int c2 = s.asInt(8);
		try {
			return new GenericDifferentialRegionScore(st, cutoff, max, id, dc, id1, id2, c1, c2);
		} catch (IOException e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	public String getConfigFileLineFormat() {
		return "scoreTable\tscoreCutoff\tcutoffIsMax\texperimentIDtifferentialScoreColNum\texp1ID\texp2ID\texp1scoreColNum\texp2scoreColNum";
	}
	
	
	@SuppressWarnings("unused")
	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 9) {
			crashWithHelpMessage(line, logger);
		}
		try {
			String st = s.asString(0);
			double cutoff = s.asDouble(1);
			boolean max = s.asBoolean(2);
			String id = s.asString(3);
			int dc = s.asInt(4);
			String id1 = s.asString(5);
			String id2 = s.asString(6);
			int c1 = s.asInt(7);
			int c2 = s.asInt(8);
		} catch(Exception e) {
			crashWithHelpMessage(line, logger);
		}
	}

	
}
