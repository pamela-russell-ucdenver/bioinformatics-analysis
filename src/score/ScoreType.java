package score;

/**
 * Types of region scores
 * @author prussell
 *
 */
public enum ScoreType {
	
	/**
	 * One score that is not a comparison of two samples
	 */
	SINGLE_REGULAR,
	
	/**
	 * One score that is a comparison of two samples
	 */
	SINGLE_DIFFERENTIAL,
	
	/**
	 * An intersection of scores that are not comparisons of two samples
	 */
	INTERSECTION_REGULAR,
	
	/**
	 * An intersection of scores that are comparisons of two samples
	 */
	INTERSECTION_DIFFERENTIAL,
	
	/**
	 * A union of scores that are not comparisons of two samples
	 */
	UNION_REGULAR,
	
	/**
	 * A union of scores that are comparisons of two samples
	 */
	UNION_DIFFERENTIAL;
	
	public String toString() {
		switch(this) {
		case INTERSECTION_DIFFERENTIAL:
			return "intersection_differential";
		case INTERSECTION_REGULAR:
			return "intersection_regular";
		case SINGLE_DIFFERENTIAL:
			return "single_differential";
		case SINGLE_REGULAR:
			return "single_regular";
		case UNION_DIFFERENTIAL:
			return "union_differential";
		case UNION_REGULAR:
			return "union_regular";
		default:
			throw new UnsupportedOperationException("Case not covered");
		}
	}
	
	/**
	 * Instantiate from the string name
	 * @param name Name of score type
	 * @return Corresponding score type object
	 */
	public static ScoreType fromString(String name) {
		for(ScoreType st : ScoreType.values()) {
			if(name.equals(st.toString())) {
				return st;
			}
		}
		throw new IllegalArgumentException("Score type " + name + " not recognized.");
	}
	
	/**
	 * @return Comma separated list of score type names
	 */
	public static String commaSeparatedList() {
		String rtrn = ScoreType.values()[0].toString();
		for(int i = 1; i < ScoreType.values().length; i++) {
			rtrn += ", " + ScoreType.values()[i].toString();
		}
		return rtrn;
	}
	
}
