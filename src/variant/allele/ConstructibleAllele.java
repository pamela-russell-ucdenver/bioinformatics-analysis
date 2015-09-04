package variant.allele;

import htsjdk.variant.variantcontext.Allele;

@SuppressWarnings("serial")
public class ConstructibleAllele extends Allele {
	
    public ConstructibleAllele(final String bases, final boolean isRef) {
        super(bases.getBytes(), isRef);
    }
	
}
