package variant.programs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

import guttmanlab.core.util.CommandLineParser;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import variant.programs.RecombinantInbredHaplotypeWriter.Parent;

public class RecombinantInbredRecombinationEvents {
	
	private static Logger logger = Logger.getLogger(RecombinantInbredRecombinationEvents.class.getName());
	
	public static String getTableHeader() {
		String rtrn = "ri_strain\t";
		rtrn += "chr\t";
		rtrn += "position1\t";
		rtrn += "position2\t";
		rtrn += "id1\t";
		rtrn += "id2\t";
		rtrn += "distance\t";
		rtrn += "position1_parent\t";
		rtrn += "position2_parent";
		return rtrn;
	}

	private class ConsecutiveRiSnps {
		
		private String chr;
		private int pos1;
		private int pos2;
		private String id1;
		private String id2;
		private Parent pos1parent;
		private Parent pos2parent;
		private String riName;
		private RecombinantInbredHaplotypeWriter rihw;
		
		public ConsecutiveRiSnps(RecombinantInbredHaplotypeWriter rw, VariantContext vc1, VariantContext vc2, String riSampleName) {
			rihw = rw;
			rihw.checkRIvariantsConsecutive(vc1, vc2);
			chr = vc1.getContig();
			pos1 = vc1.getStart();
			pos2 = vc2.getStart();
			id1 = vc1.getID();
			id2 = vc2.getID();
			riName = riSampleName;
			pos1parent = rw.riVariantOrigin(vc1.getGenotypes(), riName);
			pos2parent = rw.riVariantOrigin(vc2.getGenotypes(), riName);
		}
		
		public boolean containsTransition() {
			if(pos1parent == Parent.NEITHER) return false;
			if(pos2parent == Parent.NEITHER) return false;
			return pos1parent != pos2parent;
		}
				
		public String getTableLine() {
			String rtrn = riName + "\t";
			rtrn += chr + "\t";
			rtrn += pos1 + "\t";
			rtrn += pos2 + "\t";
			rtrn += id1 + "\t";
			rtrn += id2 + "\t";
			rtrn += Integer.valueOf(pos2 - pos1).toString() + "\t";
			rtrn += rihw.parentName(pos1parent) + "\t";
			rtrn += rihw.parentName(pos2parent) + "\t";
			return rtrn;
		}
		
	}
	
	private void writeTransitions(RecombinantInbredHaplotypeWriter rw, String outFile) throws IOException {
		logger.info("Writing transition regions to " + outFile + "...");
		BufferedWriter w = new BufferedWriter(new FileWriter(outFile));
		w.write(getTableHeader() + "\n");
		for(String riSample : rw.sampleNames) {
			logger.info("");
			logger.info(riSample);
			CloseableIterator<VariantContext> iter = rw.riVcfReader.iterator();
			VariantContext first = null;
			VariantContext second = iter.next();
			while(iter.hasNext()) {
				first = second;
				second = iter.next();
				if(!first.getContig().equals(second.getContig())) {
					logger.info("Finished " + first.getContig());
					continue;
				}
				ConsecutiveRiSnps crs = new ConsecutiveRiSnps(rw, first, second, riSample);
				if(crs.containsTransition()) {
					w.write(crs.getTableLine() + "\n");
				}
			}
			iter.close();
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-rv", "VCF file including parental strains and recombinant inbred strains", true);
		p.addStringArg("-p1", "Parental sequence 1 name", true);
		p.addStringArg("-p2", "Parental sequence 2 name", true);
		p.addStringArg("-o", "Output table", true);
		p.parse(args);
		String riVcf = p.getStringArg("-rv");
		String parent1 = p.getStringArg("-p1");
		String parent2 = p.getStringArg("-p2");
		String outTable = p.getStringArg("-o");
		
		RecombinantInbredHaplotypeWriter r = new RecombinantInbredHaplotypeWriter(riVcf, null, null, parent1, parent2, null, null);
		new RecombinantInbredRecombinationEvents().writeTransitions(r, outTable);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
