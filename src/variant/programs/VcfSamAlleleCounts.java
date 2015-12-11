package variant.programs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import annotation.OverlapUtils;
import variant.allele.AlleleCounts;
import variant.allele.VariantAlleleCounts;
import variant.allele.VcfSamAllele;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class VcfSamAlleleCounts implements VariantAlleleCounts<VariantContext, SAMFileReader> {
	
	private static Logger logger = Logger.getLogger(VcfSamAlleleCounts.class.getName());
	private OverlapUtils overlap;
	
	/**
	 * @param annotationBed Bed file of transcripts of interest
	 * @param referenceSizes Reference genome chromosome sizes
	 * @throws IOException
	 */
	public VcfSamAlleleCounts(String annotationBed, String referenceSizes) throws IOException {
		overlap = new OverlapUtils(annotationBed, referenceSizes);
	}
	
	@Override
	public AlleleCounts getAlleleCounts(VariantContext variant, SAMFileReader dataset) {
		VcfSamAllele vsa = new VcfSamAllele();
		AlleleCounts rtrn = new AlleleCounts();
		SAMRecordIterator iter = dataset.queryOverlapping(variant.getContig(), variant.getStart(), variant.getEnd());
		while(iter.hasNext()) {
			try {
				SAMRecord record = iter.next();
				Allele allele = vsa.getAllele(variant, record);
				if(allele == null) continue;
				rtrn.add(allele);
			} catch(SAMFormatException e) {
				logger.warn("Caught exception, skipping read: " + e.getMessage());
			}
		}
		iter.close();
		return rtrn;
	}
	
	private static String getTableHeader() {
		String rtrn = "";
		rtrn += "chr\t";
		rtrn += "variant_start\t";
		rtrn += "variant_alleles\t";
		rtrn += "overlapping_features\t";
		rtrn += "overlapping_feature_strands\t";
		rtrn += "overlapping_feature_transcript_coords_from_5prime\t";
		rtrn += "overlapping_feature_exon_numbers_from_5prime\t";
		rtrn += "allele_counts\t";
		rtrn += "allele_proportions\t";
		return rtrn;
	}
	
	private static String getTableLine(VariantContext vc, AlleleCounts ac, FeatureCollection<Gene> overlappers) {
		Map<Allele, Integer> counts = ac.getCounts();
		if(counts.isEmpty()) {
			throw new IllegalArgumentException("Counts map is empty");
		}
		String line = vc.getContig() + "\t";
		line += vc.getStart() + "\t";
		Iterator<Allele> ai = vc.getAlleles().iterator();
		line += ai.next().getBaseString();
		while(ai.hasNext()) {
			line += "," + ai.next().getBaseString();
		}
		line += "\t";

		String overlapperIDs = "";
		String overlapperTranscriptStrands = "";
		String overlapperTranscriptCoords = "";
		String overlapperExonNumbers = "";
		net.sf.samtools.util.CloseableIterator<Gene> overlapIter = overlappers.sortedIterator();
		if(overlappers.getNumAnnotations() == 0) {
			overlapperIDs += "-";
			overlapperTranscriptStrands += "-";
			overlapperTranscriptCoords += "-";
			overlapperExonNumbers += "-";
		} else {
			while(overlapIter.hasNext()) {
				Gene gene = overlapIter.next();
				overlapperIDs += gene.getName() + ";";
				overlapperTranscriptStrands += gene.getOrientation().toString() + ";";
				overlapperTranscriptCoords += gene.getRelativePositionFrom5PrimeOfFeature(vc.getStart()) + ";";
				overlapperExonNumbers += OverlapUtils.exonNumberFrom5Prime(gene, vc.getContig(), vc.getStart()) + ";";
			}
		overlapIter.close();
		}
		line += overlapperIDs + "\t";
		line += overlapperTranscriptStrands + "\t";
		line += overlapperTranscriptCoords + "\t";
		line += overlapperExonNumbers + "\t";
		
		line += "\t" + AlleleCounts.asString(counts) + "\t";
		line += AlleleCounts.asString(ac.getProportions()) + "\t";
		
		return line;
	}
	
	private void writeCountTable(String vcf, String bam, String out) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(out));
		writer.write(getTableHeader() + "\n");
		VCFFileReader vcfReader = new VCFFileReader(new File(vcf));
		SAMFileReader samReader = new SAMFileReader(new File(bam));
		CloseableIterator<VariantContext> vcIter = vcfReader.iterator();
		int numDone = 0;
		while(vcIter.hasNext()) {
			numDone++;
			if(numDone % 100000 == 0) {
				logger.info("Finished " + numDone + " variants.");
			}
			VariantContext vc = vcIter.next();
			AlleleCounts ac = getAlleleCounts(vc, samReader);
			if(!ac.isEmpty()) {
				writer.write(getTableLine(vc, ac, overlap.getOverlappers(vc.getContig(), vc.getStart(), vc.getEnd() + 1)) + "\n");
			}
		}
		vcfReader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-v", "VCF file", true);
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-o", "Output file", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addStringArg("-t", "Bed file of transcripts", true);
		p.addStringArg("-s", "Reference size file", true);
		p.parse(args);
		String vcf = p.getStringArg("-v");
		String bam = p.getStringArg("-b");
		String out = p.getStringArg("-o");
		String bed = p.getStringArg("-t");
		String sizes = p.getStringArg("-s");
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
			VcfSamAllele.logger.setLevel(Level.DEBUG);
		}
		
		VcfSamAlleleCounts vsac = new VcfSamAlleleCounts(bed, sizes);
		vsac.writeCountTable(vcf, bam, out);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
