package counts;

import guttmanlab.core.util.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.normalize.TranscriptAverageNormalization;
import nextgen.core.utils.AnnotationUtils;

/**
 * Map regions of interest to parent genes in an annotation
 * If there are several sub-regions within a parent, combine them into one annotation
 * Calculate enrichment of regions over parents
 * Enrichment is the ratio of average read depth
 * @author prussell
 *
 */
public class SingleSampleRegionEnrichmentOverParent {
	
	private TranscriptAverageNormalization data;
	private Map<Annotation, Annotation> childToParent;
	private static Logger logger = Logger.getLogger(SingleSampleRegionEnrichmentOverParent.class.getName());
	private static boolean mergeOverlappingParents = true;
	private static boolean subtractRegionsFromParent = false;
	private static boolean pairedEnd = true;
	private static String outDir = null;
	
	/**
	 * @param bamFile Alignment data
	 * @param parentAnnotationBed Bed file of parent annotations, also used to make transcriptome space
	 * @param regionBed Bed file of sub-regions of interest within parent annotations
	 * @param mergeOverlappingParents Instead of mapping child to largest parent, merge overlapping parents so child is mapped to merged annotation
	 * @param subtractRegionsFromParent For enrichment calculation, subtract sub-regions from parents, i.e., do not count sub-region as part of parent
	 * @throws IOException
	 */
	private SingleSampleRegionEnrichmentOverParent(String bamFile, String parentAnnotationBed, String regionBed) throws IOException {
		
		// Get the parent genes
		Map<String, Collection<Gene>> parentGenes = null;
		logger.info("");
		logger.info("Reading parent genes from " + parentAnnotationBed);
		Map<String, Collection<Gene>> fullParents = BEDFileParser.loadDataByChr(parentAnnotationBed);
		for(String chr : fullParents.keySet()) {
			logger.info(chr + "\t" + fullParents.get(chr).size());
		}
		if(mergeOverlappingParents) {
			// Merge overlapping parents
			logger.info("");
			logger.info("Merging overlapping parent genes");
			parentGenes = AnnotationUtils.collapseOverlappers(fullParents, false);
			for(String chr : fullParents.keySet()) {
				logger.info(chr + "\t" + parentGenes.get(chr).size());
			}
		} else {
			parentGenes = fullParents;
		}
		// Make transcriptome space
		logger.info("");
		logger.info("Making transcriptome space for parent genes");
		TranscriptomeSpace ts = new TranscriptomeSpace(parentGenes);
		// Make alignment model
		logger.info("");
		logger.info("Making alignment model");
		AlignmentModel alignmentData = new AlignmentModel(bamFile, ts, pairedEnd);

		// Get the regions of interest
		logger.info("");
		logger.info("Getting regions from " + regionBed);
		Map<String, Collection<Gene>> regions = BEDFileParser.loadDataByChr(regionBed);
		for(String chr : fullParents.keySet()) {
			logger.info(chr + "\t" + regions.get(chr).size());
		}
		
		// Map regions to parent genes
		logger.info("");
		logger.info("Mapping regions to parent genes");
		childToParent = new TreeMap<Annotation, Annotation>();
		Map<Gene, Gene> childToParentGene = AnnotationUtils.mapChildToLargestParent(regions, parentGenes);
		childToParent.putAll(childToParentGene);
		if(!subtractRegionsFromParent) {
			data = new TranscriptAverageNormalization(alignmentData, childToParent);
		} else {
			logger.info("Merging regions under single parent gene");
			Map<Annotation, Collection<Annotation>> parentToChildren = new TreeMap<Annotation, Collection<Annotation>>();
			for(Annotation child : childToParent.keySet()) {
				Annotation parent = childToParent.get(child);
				if(!parentToChildren.containsKey(parent)) {
					parentToChildren.put(parent, new TreeSet<Annotation>());
				}
				parentToChildren.get(parent).add(child);
			}
			Map<Gene, Gene> joinedChildToParentGene = new TreeMap<Gene, Gene>();
			for(Annotation parent : parentToChildren.keySet()) {
				Gene mergedChildren = new Gene(parentToChildren.get(parent), "regions_" + parent.getName());
				joinedChildToParentGene.put(mergedChildren, new Gene(parent));
			}
			logger.info("Subtracting merged regions from parent genes");
			childToParent.clear();		
			childToParent.putAll(joinedChildToParentGene);
			for(Annotation child : childToParent.keySet()) {
				Annotation parent = childToParent.get(child);
				Annotation parentMinusChild = parent.minus(child);
				parentMinusChild.setName(parent.getName() + "_minus_regions");
				childToParent.put(child, parentMinusChild);
			}
			data = new TranscriptAverageNormalization(alignmentData, childToParent);
		}

		// Write bed files of regions and parent genes used
		logger.info("");
		logger.info("There are " + childToParent.keySet().size() + " regions and " + childToParent.values().size() + " parent genes.");
		String outRegion = outDir + "/regions_used.bed";
		String outParent = outDir + "/parent_genes_used.bed";
		logger.info("Writing actual regions used to " + outRegion);
		BEDFileParser.writeBED(outRegion, childToParent.keySet());
		logger.info("Writing actual parent genes used to " + outParent);
		Collection<Annotation> parentsToWrite = new TreeSet<Annotation>();
		parentsToWrite.addAll(childToParent.values());
		BEDFileParser.writeBED(outParent, parentsToWrite);
		
	}

	/**
	 * Get enrichment of region over its parent gene
	 * @param region The region
	 * @return Enrichment of region over parent gene
	 */
	private double getEnrichment(Annotation region) {
		return data.getNormalizedCount(region);
	}
	
	private static double log2(double x) {
		return Math.log(x) / Math.log(2);
	}
	
	/**
	 * Write enrichments to a file
	 * @throws IOException
	 */
	private void writeEnrichmentsToFile() throws IOException {
		String outFile = outDir + "/region_enrichments_over_parent_gene.out";
		logger.info("");
		logger.info("Writing enrichments to " + outFile);
		FileWriter w = new FileWriter(outFile);
		w.write("Region\tParent\tlog2(enrichment)\n");
		for(Annotation region : childToParent.keySet()) {
			logger.info(region.getName());
			double enrichment = log2(getEnrichment(region));
			w.write(region.getName() + "\t" + childToParent.get(region).getName() + "\t" + enrichment + "\n");
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-rb", "Bed file of regions", true);
		p.addStringArg("-pb", "Bed file of parent genes", true);
		p.addBooleanArg("-mp", "Merge overlapping parent genes", false, mergeOverlappingParents);
		p.addBooleanArg("-sr", "Subtract regions from parent gene before computing ratio", false, subtractRegionsFromParent);
		p.addStringArg("-o", "Output directory", true);
		p.addBooleanArg("-p", "Treat paired reads as full fragments (fill in full fragment)", false, pairedEnd);
		p.parse(args);
		String bam = p.getStringArg("-b");
		String regionBed = p.getStringArg("-rb");
		String geneBed = p.getStringArg("-pb");
		mergeOverlappingParents = p.getBooleanArg("-mp");
		subtractRegionsFromParent = p.getBooleanArg("-sr");
		outDir = p.getStringArg("-o");
		pairedEnd = p.getBooleanArg("-p");
		
		SingleSampleRegionEnrichmentOverParent ss = new SingleSampleRegionEnrichmentOverParent(bam, geneBed, regionBed);
		ss.writeEnrichmentsToFile();
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
	
}
