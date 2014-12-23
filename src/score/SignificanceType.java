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
	SINGLE_SAMPLE,
	
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
	EITHER_SAMPLE_UP;
	
	public String toString() {
		switch(this) {
		case EITHER_SAMPLE_UP:
			return "either_sample_up";
		case SAMPLE_1_UP:
			return "sample_1_up";
		case SAMPLE_2_UP:
			return "sample_2_up";
		case SINGLE_SAMPLE:
			return "single_sample";
		default:
			throw new UnsupportedOperationException("Case not covered");
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
