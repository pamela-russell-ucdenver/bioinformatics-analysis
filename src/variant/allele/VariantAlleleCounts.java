package variant.allele;

/**
 * Get allele counts for the variant in a dataset
 * @author prussell
 *
 * @param <T> Class encoding a variant, e.g. VariantContext
 * @param <S> Class encoding a dataset, e.g. SAMFileReader
 */
public interface VariantAlleleCounts<T,S> {
	
	/**
	 * Get allele counts for the variant in the dataset
	 * @param variant Variant
	 * @param dataset Dataset
	 * @return Allele counts
	 */
	public AlleleCounts getAlleleCounts(T variant, S dataset);
	
}
