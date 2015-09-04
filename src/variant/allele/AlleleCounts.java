package variant.allele;

import htsjdk.variant.variantcontext.Allele;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class AlleleCounts {
	
	private Map<Allele, Integer> counts;
	private int totalCount;
	
	/**
	 * No initial alleles or counts
	 */
	public AlleleCounts() {
		counts = new HashMap<Allele, Integer>();
		totalCount = 0;
	}
	
	/**
	 * Initialize the possible alleles without any counts yet
	 * @param alleles Possible alleles
	 */
	public AlleleCounts(Set<Allele> alleles) {
		this();
		for(Allele allele : alleles) {
			counts.put(allele, Integer.valueOf(0));
		}
	}
	
	/**
	 * Initialize the possible alleles and initial counts
	 * @param alleles Alleles with multiplicities
	 */
	public AlleleCounts(List<Allele> alleles) {
		this();
		for(Allele allele : alleles) {
			add(allele);
		}
	}
	
	/**
	 * Initialize the possible alleles and initial counts
	 * @param countsByAllele Map of allele to initial count
	 */
	public AlleleCounts(Map<Allele, Integer> countsByAllele) {
		totalCount = 0;
		counts = countsByAllele;
	}
	
	/**
	 * Add another allele to the counts
	 * @param allele New allele
	 */
	public void add(Allele allele) {
		if(!counts.containsKey(allele)) {
			counts.put(allele, Integer.valueOf(0));
		}
		increment(allele);
	}
	
	/**
	 * Increment count for the allele
	 * @param allele Allele
	 */
	private void increment(Allele allele) {
		counts.put(allele, Integer.valueOf(counts.get(allele).intValue()+1));
		totalCount++;
	}
	
	/**
	 * Get the set of possible alleles
	 * @return Possible alleles
	 */
	public Set<Allele> getAlleles() {
		return counts.keySet();
	}
	
	/**
	 * Get total count of all alleles
	 * @return Total count
	 */
	public int getTotalCount() {
		return totalCount;
	}
	
	/**
	 * Get count for one allele
	 * @param allele Allele
	 * @return Allele count
	 */
	public int getCount(Allele allele) {
		return counts.get(allele).intValue();
	}
	
	/**
	 * Get counts for all alleles
	 * @return Map of allele to count
	 */
	public Map<Allele, Integer> getCounts() {
		return counts;
	}
	
	/**
	 * Get proportions for all alleles
	 * @return Map of allele to proportion
	 */
	public Map<Allele, Float> getProportions() {
		Map<Allele, Float> rtrn = new HashMap<Allele,Float>();
		for(Allele allele : getAlleles()) {
			rtrn.put(allele, Float.valueOf(getProportion(allele)));
		}
		return rtrn;
	}
	
	/**
	 * Get string representation of counts
	 * @param counts Map of allele to count or other statistic
	 * @return String representation
	 */
	public static String asString(Map<Allele, ? extends Number> counts) {
		Iterator<Allele> iter = counts.keySet().iterator();
		Allele first = iter.next();
		String rtrn = first.getBaseString() + ":" + counts.get(first).toString();
		while(iter.hasNext()) {
			Allele allele = iter.next();
			rtrn += ";" + allele.getBaseString() + ":" + counts.get(allele).toString();
		}
		return rtrn;
	}
	
	public boolean isEmpty() {
		return counts.isEmpty();
	}
	
	/**
	 * Get proportion for one allele
	 * @param allele Allele
	 * @return Allele proportion
	 */
	public float getProportion(Allele allele) {
		return (float) getCount(allele) / (float) totalCount;
	}
	
	/**
	 * Determine if the variant is heterozygous
	 * @return True iff variant is called heterozygous
	 */
	public boolean isHet() {
		throw new UnsupportedOperationException("Not implemented");
	}
	
}
