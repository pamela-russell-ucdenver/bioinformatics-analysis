package score;

/**
 * Ways in which a gene can be significant for a particular score
 * @author prussell
 *
 */
public enum SignificanceType {
	
	/**
	 * There is only one sample so a gene is just significant or not for the score
	 */
	SINGLE_SAMPLE_SIGNIFICANT,
	
	/**
	 * There is only one sample so a gene is just significant or not for the score
	 */
	SINGLE_SAMPLE_NOT_SIGNIFICANT,
	
	/**
	 * There are two samples and sample 1 is up with respect to sample 2
	 */
	SAMPLE_1_UP,
	
	/**
	 * There are two samples and sample 2 is up with respect to sample 1
	 */	
	SAMPLE_2_UP,
	
	/**
	 * There are two samples and either sample is up
	 */
	EITHER_SAMPLE_UP,
	
	/**
	 * Not significant
	 */
	TWO_SAMPLE_NOT_SIGNIFICANT;
	
	public String toString() {
		switch(this) {
		case EITHER_SAMPLE_UP:
			return "either_sample_up";
		case SAMPLE_1_UP:
			return "sample_1_up";
		case SAMPLE_2_UP:
			return "sample_2_up";
		case SINGLE_SAMPLE_SIGNIFICANT:
			return "single_sample_significant";
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			return "single_sample_not_significant";
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			return "two_sample_not_significant";
		default:
			throw new UnsupportedOperationException("Case not covered");
		}
	}

	/**
	 * Convert a two-sample significance type to its single sample equivalent
	 * @param significanceType Two-sample significance type
	 * @return Equivalent for a single sample score
	 */
	public static SignificanceType twoSampleToSingleSampleAnalog(SignificanceType significanceType) {
		switch(significanceType) {
		case EITHER_SAMPLE_UP:
			return SINGLE_SAMPLE_SIGNIFICANT;
		case SAMPLE_1_UP:
			throw new IllegalArgumentException(significanceType.toString() + " has no single sample equivalent");
		case SAMPLE_2_UP:
			throw new IllegalArgumentException(significanceType.toString() + " has no single sample equivalent");
		case SINGLE_SAMPLE_NOT_SIGNIFICANT:
			throw new IllegalArgumentException(significanceType.toString() + " is already a single sample significance type");
		case SINGLE_SAMPLE_SIGNIFICANT:
			throw new IllegalArgumentException(significanceType.toString() + " is already a single sample significance type");
		case TWO_SAMPLE_NOT_SIGNIFICANT:
			return SINGLE_SAMPLE_NOT_SIGNIFICANT;
		default:
			throw new IllegalArgumentException("Significance type " + significanceType.toString() + " not implemented.");
		}
	}
	
	/**
	 * Create from a string name
	 * @param name Name
	 * @return Corresponding significance type object
	 */
	public static SignificanceType fromString(String name) {
		for(SignificanceType st : SignificanceType.values()) {
			if(name.equals(st.toString())) {
				return st;
			}
		}
		throw new IllegalArgumentException("Significance type " + name + " not recognized.");
	}
	
	/**
	 * @return Comma separated list of the names
	 */
	public static String commaSeparatedList() {
		String rtrn = SignificanceType.values()[0].toString();
		for(int i = 1; i < SignificanceType.values().length; i++) {
			rtrn += "," + SignificanceType.values()[i].toString();
		}
		return rtrn;
	}

	
}
