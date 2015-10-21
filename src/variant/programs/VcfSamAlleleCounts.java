package variant.programs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import variant.allele.AlleleCounts;
import variant.allele.VariantAlleleCounts;
import variant.allele.VcfSamAllele;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import guttmanlab.core.util.CommandLineParser;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class VcfSamAlleleCounts implements VariantAlleleCounts<VariantContext, SAMFileReader> {
	
	private static Logger logger = Logger.getLogger(VcfSamAlleleCounts.class.getName());
	
	public VcfSamAlleleCounts() {}
	
	@Override
	public AlleleCounts getAlleleCounts(VariantContext variant, SAMFileReader dataset) {
		VcfSamAllele vsa = new VcfSamAllele();
		AlleleCounts rtrn = new AlleleCounts();
		SAMRecordIterator iter = dataset.queryOverlapping(variant.getContig(), variant.getStart(), variant.getEnd());
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			Allele allele = vsa.getAllele(variant, record);
			if(allele == null) continue;
			rtrn.add(allele);
		}
		iter.close();
		return rtrn;
	}
	
	private static String getTableLine(VariantContext vc, AlleleCounts ac) {
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
		line += "\t" + AlleleCounts.asString(counts) + "\t";
		line += AlleleCounts.asString(ac.getProportions());
		return line;
	}
	
	private void writeCountTable(String vcf, String bam, String out) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(out));
		
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
				writer.write(getTableLine(vc, ac) + "\n");
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
		p.parse(args);
		String vcf = p.getStringArg("-v");
		String bam = p.getStringArg("-b");
		String out = p.getStringArg("-o");
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
			VcfSamAllele.logger.setLevel(Level.DEBUG);
		}
		
		VcfSamAlleleCounts vsac = new VcfSamAlleleCounts();
		vsac.writeCountTable(vcf, bam, out);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
