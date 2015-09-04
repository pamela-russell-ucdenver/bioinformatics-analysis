package variant.allele;


import org.apache.log4j.Logger;

import sam.SAMRecordUtils;
import variant.VCFUtils;
import net.sf.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VcfSamAllele implements VariantAllele<VariantContext, SAMRecord> {
	
	public static Logger logger = Logger.getLogger(VcfSamAllele.class.getName());
	
	public VcfSamAllele() {}
	
	private static boolean isSNP(VariantContext variant) {
		if(variant.getReference().length() != 1) {
			return false;
		}
		for(Allele altAllele : variant.getAlternateAlleles()) {
			if(altAllele.length() != 1) {
				return false;
			}
		}
		return true;
	}
	
	@Override
	public Allele getAllele(VariantContext variant, SAMRecord record) {
		if(!isSNP(variant)) {
			return null;
		}
		byte[] readBases = record.getReadBases(); // Always returns reference strand regardless of mapped strand
		byte[] alleleBases = new byte[1];
		int zeroBasedVariantStart = VCFUtils.vcfToZeroBased(variant.getStart());
		int readPos = SAMRecordUtils.getReadPositionAtReferencePosition(record, zeroBasedVariantStart);
		if(readPos < 0) return null;
		byte base = readBases[readPos];
		alleleBases[0] = base;
		// TODO boolean is "isRef". We don't know. Can't even get this info from cigar (match/mismatch have same element)
		return new ConstructibleAllele(new String(alleleBases), true);
	}

}
