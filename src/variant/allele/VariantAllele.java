package variant.allele;

import htsjdk.variant.variantcontext.Allele;

/**
 * Identify the allele present in a record for a specific variant
 * @author prussell
 *
 * @param <T> Class encoding a variant, e.g. VariantContext
 * @param <S> Class encoding a record, e.g. SAMRecord
 */
public interface VariantAllele<T,S> {
	
	/**
	 * Get the allele for the variant in the record, or null if not possible
	 * @param variant Variant
	 * @param record Record
	 * @return Allele or null
	 */
	public Allele getAllele(T variant, S record);
	
}
