package variant.genotype;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;

@SuppressWarnings("serial")
public class ConstructibleGenotypesContext extends GenotypesContext {

    /**
     * Create a fully resolved GenotypeContext containing genotypes, sample lookup table,
     * and sorted sample names
     *
     * @param genotypes our genotypes in arbitrary
     * @param sampleNameToOffset map optimized for efficient lookup.  Each genotype in genotypes must have its
     * sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that
     * genotype in the vector of genotypes
     * @param sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical
     * order.
     */
    protected ConstructibleGenotypesContext(final ArrayList<Genotype> genotypes,
                             final Map<String, Integer> sampleNameToOffset,
                             final List<String> sampleNamesInOrder) {
        super(genotypes, sampleNameToOffset, sampleNamesInOrder);
    }


}
