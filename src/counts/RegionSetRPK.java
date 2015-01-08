package counts;

import java.io.IOException;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;
import net.sf.samtools.util.CloseableIterator;

import org.apache.log4j.Logger;

/**
 * Get reads per kb over a set of regions
 * @author prussell
 *
 */
public class RegionSetRPK {

	private AnnotationCollection<? extends Annotation> data;
	private String chrSizeFile;
	private static Logger logger = Logger.getLogger(RegionSetRPK.class.getName());
	
	private RegionSetRPK(String bamFile, String chrSizes) {
		data = BAMFragmentCollectionFactory.createFromBam(bamFile);
		chrSizeFile = chrSizes;
	}
	
	private double getTotalCount(String regionBed) throws IOException {
		double rtrn = 0;
		AnnotationCollection<Gene> regions = BEDFileIO.loadFromFile(regionBed, chrSizeFile);
		int numGenes = regions.getNumAnnotations();
		CloseableIterator<Gene> iter = regions.sortedIterator();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		logger.info("");
		logger.info("Adding counts for " + numGenes + " genes...");
		while(iter.hasNext()) {
			Gene region = iter.next();
			rtrn += data.numOverlappers(region, true);
			countLogger.advance();
		}
		iter.close();
		logger.info("There are " + rtrn + " reads total.");
		return rtrn;
	}
	
	private double getRPK(String regionBed) throws IOException {
		AnnotationCollection<Gene> regions = BEDFileIO.loadFromFile(regionBed, chrSizeFile);
		double totalSize = 0;
		CloseableIterator<Gene> iter = regions.sortedIterator();
		while(iter.hasNext()) {
			Gene region = iter.next();
			totalSize += region.size();
		}
		iter.close();
		logger.info("Total size of regions is " + totalSize + ".");
		double totalCount = getTotalCount(regionBed);
		double kb = totalSize / 1000;
		double rtrn = totalCount / kb;
		logger.info("Reads per kilobase = " + rtrn + ".");
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		logger.warn("PROGRAM DOES NOT COLLAPSE OVERLAPPERS");
		logger.warn("MAKE SURE YOUR REGIONS DO NOT OVERLAP");
		logger.warn("");
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-c", "Chr size file", true);
		p.addStringArg("-r", "Region bed file", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-r");
		String chrSizes = p.getStringArg("-c");
		
		RegionSetRPK re = new RegionSetRPK(bamFile, chrSizes);
		
		@SuppressWarnings("unused")
		double rpk = re.getRPK(bedFile);
		
		logger.info("");
		logger.info("All done.");
		
	}
	

	
}
