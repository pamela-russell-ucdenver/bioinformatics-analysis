package score;

import translation.DifferentialTranslationalEfficiency;
import translation.TranslationalEfficiencyFromBam;
import expression.DifferentialExpressionCuffdiff;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.util.StringParser;

/**
 * Static factory methods to create region score objects from config file lines
 * @author prussell
 *
 */
public class RegionScoreFactory {
	
	public static final String DIFF_EXP_CUFFDIFF = "diff_exp_cuffdiff";
	public static final String GENERIC_REGION_SCORE = "generic";
	public static final String GENERIC_DIFF_REGION_SCORE = "generic_diff";
	public static final String TRANSLATIONAL_EFFICIENCY = "translational_efficiency";
	public static final String DIFF_TRANSLATIONAL_EFFICIENCY = "diff_translational_efficiency";
	
	/**
	 * Create a score object specified by a config file line
	 * The name of the score is in the specified field
	 * The information needed to create the object is in the fields after the name
	 * There can be other fields before the name that are ignored
	 * @param line Config file line
	 * @param scoreNameField The zero-based field containing the name of the score
	 * @return The score object described by the information after the score name
	 */
	public static RegionScore<? extends Annotation> createScoreFromConfigFileLine(String line, int scoreNameField) {
		
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() < scoreNameField - 1) {
			throw new IllegalArgumentException("Line does not have enough fields to get score name from field " + scoreNameField + ": " + line);
		}
		String scoreName = s.asString(scoreNameField);
		String lineSuffix = s.removeFirstTokens(scoreNameField + 1);
		
		if(scoreName.equals(DIFF_EXP_CUFFDIFF)) {
			return new DifferentialExpressionCuffdiff().createFromConfigFileLine(lineSuffix);
		}
		
		if(scoreName.equals(DIFF_TRANSLATIONAL_EFFICIENCY)) {
			return new DifferentialTranslationalEfficiency().createFromConfigFileLine(lineSuffix);
		}
		
		if(scoreName.equals(GENERIC_DIFF_REGION_SCORE)) {
			return new GenericDifferentialRegionScore().createFromConfigFileLine(lineSuffix);
		}
		
		if(scoreName.equals(GENERIC_REGION_SCORE)) {
			return new GenericRegionScore().createFromConfigFileLine(lineSuffix);
		}
		
		if(scoreName.equals(TRANSLATIONAL_EFFICIENCY)) {
			return new TranslationalEfficiencyFromBam().createFromConfigFileLine(lineSuffix);
		}
		
		throw new IllegalArgumentException("Score " + scoreName + " not supported.");
		
	}

	/**
	 * Create a differential score object specified by a config file line
	 * The name of the score is in the specified field
	 * The information needed to create the object is in the fields after the name
	 * There can be other fields before the name that are ignored
	 * @param line Config file line
	 * @param scoreNameField The zero-based field containing the name of the score
	 * @return The score object described by the information after the score name
	 */
	public static DifferentialRegionScore<? extends Annotation> createDiffScoreFromConfigFileLine(String line, int scoreNameField) {
		
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() < scoreNameField - 1) {
			throw new IllegalArgumentException("Line does not have enough fields to get score name from field " + scoreNameField + ": " + line);
		}
		String scoreName = s.asString(scoreNameField);
		String lineSuffix = s.removeFirstTokens(scoreNameField + 1);
		
		if(scoreName.equals(DIFF_EXP_CUFFDIFF)) {
			return (DifferentialRegionScore<Gene>) new DifferentialExpressionCuffdiff().createFromConfigFileLine(lineSuffix);
		}
		
		if(scoreName.equals(DIFF_TRANSLATIONAL_EFFICIENCY)) {
			return (DifferentialRegionScore<Gene>) new DifferentialTranslationalEfficiency().createFromConfigFileLine(lineSuffix);
		}
		
		if(scoreName.equals(GENERIC_DIFF_REGION_SCORE)) {
			return (DifferentialRegionScore<Gene>) new GenericDifferentialRegionScore().createFromConfigFileLine(lineSuffix);
		}
		
		throw new IllegalArgumentException("Score " + scoreName + " not supported or is not a differential score.");
		
	}


	
}
	
