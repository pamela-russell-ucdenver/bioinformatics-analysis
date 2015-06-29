package counts;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Random;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.util.CommandLineParser;

public class RegionEnrichmentOverCoordSpace {
	
	private AnnotationCollection<? extends Annotation> data;
	private CoordinateSpace coordSpace;
	private Random rand;
	private static Logger logger = Logger.getLogger(RegionEnrichmentOverCoordSpace.class.getName());
	
	private RegionEnrichmentOverCoordSpace(String bamFile, String chrSizes) {
		rand = new Random();
		data = BAMFragmentCollectionFactory.createFromBam(bamFile);
		coordSpace = new CoordinateSpace(chrSizes);
	}
	
	private Annotation getCoordSpacePosition(long overallTotalPosition) {
		Map<String, Integer> refLengths = coordSpace.getRefSeqLengths();
		String chr = null;
		long currTotalLength = 0;
		long previousTotalLength = 0;
		Iterator<String> chrIter = refLengths.keySet().iterator();
		while(currTotalLength < overallTotalPosition) {
			try {
				chr = chrIter.next();
			} catch(NoSuchElementException e) {
				logger.error("Ran out of chromosomes for random genome position " + overallTotalPosition);
				throw(e);
			}
			previousTotalLength = currTotalLength;
			currTotalLength += refLengths.get(chr);
		}
		int rtrnPos = (int) (overallTotalPosition - previousTotalLength);
		return new SingleInterval(chr, rtrnPos, rtrnPos);
	}
	
	private SingleInterval getRandomRegion(int regionLength) {
		long totalRegionLength = coordSpace.getTotalReferenceLength();
		int numTries = 0;
		while(numTries < 1000) {
			long overallStartPos = Math.abs(rand.nextLong()) % totalRegionLength;
			Annotation randStart = getCoordSpacePosition(overallStartPos);
			String chr = randStart.getReferenceName();
			int start = randStart.getReferenceStartPosition();
			int end = start + regionLength;
			try {
				if(end > coordSpace.getRefSeqLengths().get(chr).intValue()) {
					numTries++;
					continue;
				}
			} catch(NullPointerException e) {
				logger.error("Null pointer exception on chromosome " + chr);
				throw(e);
			}
			return new SingleInterval(chr, start, end);
		}
		throw new IllegalStateException("Tried 1000 random regions and couldn't fit one inside a chromosome. Check requested region length: " + regionLength);
	}
	
	private List<SingleInterval> getRandomRegions(int regionLength, int numRegions) {
		List<SingleInterval> rtrn = new ArrayList<SingleInterval>();
		while(rtrn.size() < numRegions) {
			rtrn.add(getRandomRegion(regionLength));
		}
		return rtrn;
	}
	
	private double getAverageCount(List<SingleInterval> regions) {
		double totalCount = 0;
		for(SingleInterval region : regions) {
			totalCount += data.numOverlappers(region, false);
		}
		return totalCount / regions.size();
	}
	
	private double getAverageCountRandomRegionsSameSize(Annotation region, int numRandomRegions) {
		return getAverageCount(getRandomRegions(region.size(), numRandomRegions));
	}
	
	private double getEnrichment(Annotation region, int numRandomRegions) {
		double rtrn = data.numOverlappers(region, false) / getAverageCountRandomRegionsSameSize(region, numRandomRegions);
		logger.info("Enrichment for " + region.toUCSC() + " is " + rtrn);
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-c", "Chr size file", true);
		p.addStringArg("-r", "Region bed file", true);
		p.addStringArg("-o", "Output table", true);
		p.addIntArg("-n", "Number of random regions per region", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-r");
		String outTable = p.getStringArg("-o");
		int numRand = p.getIntArg("-n");
		String chrSizes = p.getStringArg("-c");
		
		RegionEnrichmentOverCoordSpace re = new RegionEnrichmentOverCoordSpace(bamFile, chrSizes);
		AnnotationCollection<Gene> regions = BEDFileIO.loadFromFile(bedFile, chrSizes);
		FileWriter writer = new FileWriter(outTable);
		
		CloseableIterator<Gene> iter = regions.sortedIterator();
		while(iter.hasNext()) {
			Gene region = iter.next();
			double enrichment = re.getEnrichment(region, numRand);
			writer.write(region.toUCSC() + "\t" + enrichment + "\n");
		}
		iter.close();
		
		writer.close();
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
