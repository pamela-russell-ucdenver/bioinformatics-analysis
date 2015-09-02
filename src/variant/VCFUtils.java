package variant;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.variant.ReferenceConfidenceVariantContextMerger;

public class VCFUtils {
	
	public static Logger logger = Logger.getLogger(VCFUtils.class.getName());

	/**
	 * Merge variants with same genomic location
	 * @param variants Variants including some with same location
	 * @param dict SAM sequence dictionary for reference genome
	 * @return Merged variant contexts
	 */
	public static Collection<VariantContext> merge(Collection<VariantContext> variants, SAMSequenceDictionary dict) {
		Map<GenomeLoc, List<VariantContext>> duplicateLocs = new TreeMap<GenomeLoc, List<VariantContext>>();
		Map<GenomeLoc, Byte> refBase = new TreeMap<GenomeLoc, Byte>(); 
		for(VariantContext vc : variants) {
			GenomeLocParser glp = new GenomeLocParser(dict);
			GenomeLoc loc = glp.createGenomeLoc(vc.getContig(), vc.getStart(), vc.getEnd(), true);
			if(!duplicateLocs.containsKey(loc)) {
				duplicateLocs.put(loc, new ArrayList<VariantContext>());
			}
			duplicateLocs.get(loc).add(vc);
			Byte base = null; //TODO do we need to set this?
			refBase.put(loc, base);
		}
		Collection<VariantContext> rtrn = new ArrayList<VariantContext>();
		for(GenomeLoc loc : duplicateLocs.keySet()) {
			for(VariantContext v : duplicateLocs.get(loc)) {
				logger.debug("MERGING\t" + v.toStringDecodeGenotypes());
			}
			VariantContext merged = ReferenceConfidenceVariantContextMerger.merge(duplicateLocs.get(loc), loc, refBase.get(loc), true);
			logger.debug("MERGED\t" + duplicateLocs.get(loc).size() + "\t" + loc.toString());
			logger.debug("MERGED_VARIANT_CONTEXT\t" + merged.toStringDecodeGenotypes());
			rtrn.add(merged);
		}
		return rtrn;
	}
	
	public static int vcfToZeroBased(int pos) {
		return pos - 1;
	}
	
	public static int zeroBasedToVcf(int pos) {
		return pos + 1;
	}

}
