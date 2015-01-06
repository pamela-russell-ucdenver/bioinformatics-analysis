package translation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import score.AbstractRegionScore;
import score.GenericRegionScore;
import score.RegionScore;
import score.SignificanceType;
import bam.BamCountRegionOverlappers;
import broad.core.math.ScanStatistics;
import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;
import guttmanlab.core.util.StringParser;

public class TranslationalEfficiency extends AbstractRegionScore<Gene> {
	
	private double ribosomeGlobalGenomeTotal;
	private double controlGlobalGenomeTotal;
	private double ribosomeGlobalExonTotal;
	private double controlGlobalExonTotal;
	private double normalizationFactor; // Ribosome total / control total
	private BAMSingleReadCollection ribosomeData;
	private BAMSingleReadCollection controlData;
	private String chrSizeFile;
	private String controlBamFile;
	private String ribosomeBamFile;
	private long totalChrSize;
	private static final double EXPRESSION_SCAN_PVAL_CUTOFF = 0.01;
	public static int TE_MIN_RAW_READS = 10; // Minimum number of reads in ribosome and control fraction to compute the TE of a CDS
	private Map<String, Double> expressionScanPvals;
	private Map<String, Double> ribosomeCdsScanPvals;
	private Map<String, Double> ribosomeCounts;
	private Map<String, Double> controlCounts;
	private boolean strandSpecific;
	private double controlGlobalGenomeLambda;
	private double ribosomeGlobalGenomeLambda;
	private String ribosomeName;
	private String controlName;
	private String experimentID;
	
	public static Logger logger = Logger.getLogger(TranslationalEfficiency.class.getName());
	
	public TranslationalEfficiency() {}
	
	/**
	 * Automatically set experiment name to ribosome sample name
	 * @param ribosomeBam Bam file of ribosome profiling sample
	 * @param controlBam Bam file of control sample
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal Optional total number of ribosome reads mapped to genome (instead of computing from data)
	 * @param controlGenomeTotal Optional total number of control reads mapped to genome (instead of computing from data)
	 * @param ribosomeExonTotal Optional total number of ribosome reads mapped to exons (instead of computing from data)
	 * @param controlExonTotal Optional total number of control reads mapped to exons (instead of computing from data)
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @throws IOException
	 */
	public TranslationalEfficiency(String ribosomeBam, String controlBam, String geneBed, String chrSizes, double ribosomeGenomeTotal, double controlGenomeTotal, double ribosomeExonTotal, double controlExonTotal, boolean isStrandSpecific) throws IOException {
		this(ribosomeBam, controlBam, geneBed, chrSizes, ribosomeGenomeTotal, controlGenomeTotal, ribosomeExonTotal, controlExonTotal, isStrandSpecific, null);
	}
	
	/**
	 * @param ribosomeBam Bam file of ribosome profiling sample
	 * @param controlBam Bam file of control sample
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal Optional total number of ribosome reads mapped to genome (instead of computing from data)
	 * @param controlGenomeTotal Optional total number of control reads mapped to genome (instead of computing from data)
	 * @param ribosomeExonTotal Optional total number of ribosome reads mapped to exons (instead of computing from data)
	 * @param controlExonTotal Optional total number of control reads mapped to exons (instead of computing from data)
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @param experimentId Experiment ID
	 * @throws IOException
	 */
	public TranslationalEfficiency(String ribosomeBam, String controlBam, String geneBed, String chrSizes, double ribosomeGenomeTotal, double controlGenomeTotal, double ribosomeExonTotal, double controlExonTotal, boolean isStrandSpecific, String experimentId) throws IOException {
		logger.info("");
		
		logger.info("Creating translational efficiency object...");
		
		// Make control name from control bam file
		StringParser t = new StringParser();
		t.parse(controlBam, "/");
		String withExtension2 = t.asString(t.getFieldCount() - 1);
		String s1 = withExtension2.replaceAll(".sorted.bam", "");
		controlName = s1.replaceAll(".bam", "");
		logger.info("Control name is " + controlName + ".");
		
		// Make sample name from ribosome bam file
		StringParser s = new StringParser();
		s.parse(ribosomeBam, "/");
		String withExtension = s.asString(s.getFieldCount() - 1);
		String s2 = withExtension.replaceAll(".sorted.bam", "");
		ribosomeName = s2.replaceAll(".bam", "");
		logger.info("Ribosome sample name is " + ribosomeName + ".");
		
		// Experiment ID
		if(experimentId == null) {
			experimentID = ribosomeName;
		} else {
			experimentID = experimentId;
		}
		
		// Set whether the libraries are strand specific (affects counts over annotations)
		strandSpecific = isStrandSpecific;
		
		// Initialize caches of expression P values and region counts
		expressionScanPvals = new HashMap<String, Double>();
		ribosomeCdsScanPvals = new HashMap<String, Double>();
		ribosomeCounts = new HashMap<String, Double>();
		controlCounts = new HashMap<String, Double>();
		
		// Save chromosome size file to use when loading annotations
		chrSizeFile = chrSizes;
		
		// Load read mapping data
		controlBamFile = controlBam;
		ribosomeBamFile = ribosomeBam;
		ribosomeData = new BAMSingleReadCollection(new File(ribosomeBamFile));
		controlData = new BAMSingleReadCollection(new File(controlBamFile));

		// Compute global read counts
		controlGlobalGenomeTotal = controlGenomeTotal;
		if(controlGlobalGenomeTotal <= 0) {
			logger.info("Computing total genome read count for control sample " + controlName + "...");
			controlGlobalGenomeTotal = controlData.getNumAnnotations();
		}
		logger.info(controlGlobalGenomeTotal + " total reads in control fraction.");
		ribosomeGlobalGenomeTotal = ribosomeGenomeTotal;
		if(ribosomeGlobalGenomeTotal <= 0) {
			logger.info("Computing total genome read count for ribosome sample " + ribosomeName + "...");
			ribosomeGlobalGenomeTotal = ribosomeData.getNumAnnotations();
		}
		logger.info(ribosomeGlobalGenomeTotal + " total reads in ribosome fraction.");

		// Compute scan distribution parameters
		CoordinateSpace chrs = new CoordinateSpace(chrSizeFile);
		totalChrSize = chrs.getTotalReferenceLength();
		controlGlobalGenomeLambda = controlGlobalGenomeTotal / totalChrSize;
		logger.info("Total chromosome size is " + totalChrSize + ". Global mapped reads in control fraction is " + controlGlobalGenomeTotal + ". Control global lambda is " + controlGlobalGenomeLambda + ".");
		ribosomeGlobalGenomeLambda = ribosomeGlobalGenomeTotal / totalChrSize;
		
		// Calculate TE normalization factor
		if(controlExonTotal < 0) {
			logger.info("Computing total exon read count for control sample " + controlName + "...");
			controlExonTotal = exonTotal(geneBed, true);
		}
		controlGlobalExonTotal = controlExonTotal;
		logger.info(controlGlobalExonTotal + " total exon reads in control fraction.");
		if(ribosomeExonTotal < 0) {
			logger.info("Computing total exon read count for ribosome sample " + ribosomeName + "...");
			ribosomeExonTotal = exonTotal(geneBed, false);
		}
		ribosomeGlobalExonTotal = ribosomeExonTotal;
		logger.info(ribosomeGlobalExonTotal + " total exon reads in ribosome fraction.");
		normalizationFactor = ribosomeGlobalExonTotal / controlGlobalExonTotal;
		logger.info("Normalization factor is " + normalizationFactor + ".");
		
		logger.info("");
		logger.info("Done creating translational efficiency object for " + ribosomeName + " and " + controlName + ".");
		
	}
	
	/**
	 * @return Ribosome sample name
	 */
	public String getRibosomeName() {
		return ribosomeName;
	}
	
	/**
	 * @return Control sample name
	 */
	public String getControlName() {
		return controlName;
	}
	
	/**
	 * @param ribosomeBam Bam file of ribosome profiling sample
	 * @param controlBam Bam file of control sample
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal Optional total number of ribosome reads mapped to genome (instead of computing from data)
	 * @param controlGenomeTotal Optional total number of control reads mapped to genome (instead of computing from data)
	 * @param ribosomeExonTotal Optional total number of ribosome reads mapped to exons (instead of computing from data)
	 * @param controlExonTotal Optional total number of control reads mapped to exons (instead of computing from data)
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @param experimentId Experiment ID
	 * @throws IOException
	 */
	public static TranslationalEfficiency factory(String ribosomeBam, String controlBam, String geneBed, String chrSizes, double ribosomeGenomeTotal, double controlGenomeTotal, double ribosomeExonTotal, double controlExonTotal, boolean isStrandSpecific, String experimentId) throws IOException {
		return new TranslationalEfficiency(ribosomeBam, controlBam, geneBed, chrSizes, ribosomeGenomeTotal, controlGenomeTotal, ribosomeExonTotal, controlExonTotal, isStrandSpecific, experimentId);
	}
	
	/**
	 * @param gene Gene
	 * @return Whether the gene contains a significant number of control reads compared to genomic background
	 */
	public boolean isExpressed(Annotation gene) {
		double expressionPval = getExpressionScanPval(gene);
		return expressionPval < EXPRESSION_SCAN_PVAL_CUTOFF;
	}
	
	/**
	 * @param gene Gene
	 * @return Whether the CDS of the gene contains a significant number of ribosome footprints compared to genomic background
	 */
	public boolean cdsHasRibosomeFootprints(Gene gene) {
		double ribosomePval = getCdsRibosomeScanPval(gene);
		return ribosomePval < EXPRESSION_SCAN_PVAL_CUTOFF;
	}
	
	/**
	 * Get translational efficiency score for a region
	 * @param region The region
	 * @param parentGene Parent gene to check for expression
	 * @return The TE score, or NaN if parent gene is not expressed or either CDS does not meet min cutoff
	 */
	public double getTE(Annotation region, Annotation parentGene) {
		if(!isExpressed(parentGene)) {
			return Double.NaN;
		}
		double ribosomeReads = Math.max(getRibosomeCount(region),1);
		if(ribosomeReads < TE_MIN_RAW_READS) {
			return Double.NaN;
		}
		double controlReads = Math.max(getControlCount(region),1);
		if(ribosomeReads < TE_MIN_RAW_READS) {
			return Double.NaN;
		}
		return (ribosomeReads / controlReads) / normalizationFactor;
	}
	
	/**
	 * Get translational efficiency score for the CDS of a gene
	 * @param gene The gene
	 * @return The TE score of the CDS or NaN if gene is not expressed
	 */
	public double getCdsTE(Gene gene) {
		Annotation cds = gene.getCodingRegion();
		return getTE(cds, gene);
	}
	
	/**
	 * Write TE scores for the CDS of each gene to a bed file
	 * @param geneBed Bed file of genes to analyze
	 * @param outputBed Bed file to write
	 * @throws IOException 
	 */
	public void writeCdsTEsToBed(String geneBed, String outputBed) throws IOException {
		logger.info("");
		FileWriter w = new FileWriter(outputBed);
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		logger.info("Writing bed file with CDS TE scores for " + genes.getNumAnnotations() + " genes in " + geneBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		while(iter.hasNext()) {
			countLogger.advance();
			Gene gene = iter.next();
			double te = getCdsTE(gene);
			String bedLine = gene.toBED(te);
			w.write(bedLine + "\n");
		}
		w.close();
		iter.close();
		logger.info("Done writing to " + outputBed + ".");
	}
	
	
	/**
	 * Write TE scores for the CDS of each gene to a table
	 * @param geneBed Bed file of genes to analyze
	 * @param outputTable Table file to write
	 * @throws IOException 
	 */
	public void writeCdsTEsToTable(String geneBed, String outputTable) throws IOException {
		logger.info("");
		FileWriter w = new FileWriter(outputTable);
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		logger.info("Writing table with CDS TE scores for " + numGenes + " genes in " + geneBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		String header = "gene\t";
		header += "global_genome_count_control\t";
		header += "global_genome_count_ribosome\t";
		header += "global_exon_count_control\t";
		header += "global_exon_count_ribosome\t";
		header += "global_lambda\t";
		header += "gene_size\t";
		header += "gene_lambda\t";
		header += "expression_scan_pval\t";
		header += "gene_count_control\t";
		header += "gene_count_ribosome\t";
		header += "cds_count_control\t";
		header += "cds_count_ribosome\t";
		header += "TE_normalization_factor\t";
		header += "TE_score_CDS\t";
		w.write(header + "\n");
		while(iter.hasNext()) {
			Gene gene = iter.next();
			countLogger.advance();
			Annotation cds = gene.getCodingRegion();
			double expPval = getExpressionScanPval(gene);
			double rCountCDS = getRibosomeCount(cds);
			double mCountCDS = getControlCount(cds);
			double mCountGene = getControlCount(gene);
			double rCountGene = getRibosomeCount(gene);
			double te = getCdsTE(gene);
			int geneSize = gene.size();
			double geneLambda = mCountGene / geneSize;
			String line = gene.getName() + "\t";
			line += controlGlobalGenomeTotal + "\t";
			line += ribosomeGlobalGenomeTotal + "\t";
			line += controlGlobalExonTotal + "\t";
			line += ribosomeGlobalExonTotal + "\t";
			line += controlGlobalGenomeLambda + "\t";
			line += geneSize + "\t";
			line += geneLambda + "\t";
			line += expPval + "\t";
			line += mCountGene + "\t";
			line += rCountGene + "\t";
			line += mCountCDS + "\t";
			line += rCountCDS + "\t";
			line += normalizationFactor + "\t";
			line += te + "\t";
			w.write(line + "\n");
		}
		w.close();
		iter.close();
		logger.info("Done writing to " + outputTable + ".");
	}
	
	
	/**
	 * Get total number of reads mapping to exons in an annotation
	 * @param geneBed Bed file of genes
	 * @param control If true, use control sample; if false, use ribosome sample
	 * @return Number of reads mapping to exons
	 * @throws IOException
	 */
	private double exonTotal(String geneBed, boolean control) throws IOException {
		BamCountRegionOverlappers b = control ? new BamCountRegionOverlappers(controlBamFile, geneBed, chrSizeFile) : new BamCountRegionOverlappers(ribosomeBamFile, geneBed, chrSizeFile);
		return b.getTotalOverlappers();
	}
	
	/**
	 * @param gene Gene
	 * @return Scan P value for read count over gene in control sample
	 */
	public double getExpressionScanPval(Annotation gene) {
		if(expressionScanPvals.containsKey(gene.toBED())) {
			return expressionScanPvals.get(gene.toBED()).doubleValue();
		}
		int controlCount = (int)getControlCount(gene);
		int geneSize = gene.size();
		double rtrn = ScanStatistics.calculatePVal(controlCount, controlGlobalGenomeLambda, geneSize, totalChrSize);
		expressionScanPvals.put(gene.toBED(), Double.valueOf(rtrn));
		return rtrn;
	}
	
	
	/**
	 * @param gene Gene
	 * @return Scan P value for ribosome footprint count over CDS of gene
	 */
	public double getCdsRibosomeScanPval(Gene gene) {
		if(ribosomeCdsScanPvals.containsKey(gene.toBED())) {
			return ribosomeCdsScanPvals.get(gene.toBED()).doubleValue();
		}
		Annotation cds = gene.getCodingRegion();
		int ribosomeCdsCount = (int)getRibosomeCount(cds);
		int cdsSize = cds.size();
		double rtrn = ScanStatistics.calculatePVal(ribosomeCdsCount, ribosomeGlobalGenomeLambda, cdsSize, totalChrSize);
		ribosomeCdsScanPvals.put(gene.toBED(), Double.valueOf(rtrn));
		return rtrn;
	}
	
	
	/**
	 * @param gene Gene
	 * @return Ribosome read count over gene
	 */
	public double getRibosomeCount(Annotation gene) {
		Gene geneToUse = new Gene(gene);
		if(!strandSpecific) {
			geneToUse.setOrientation(Strand.BOTH);
		}
		if(ribosomeCounts.containsKey(geneToUse.toBED())) {
			return ribosomeCounts.get(geneToUse.toBED()).doubleValue();
		}
		int ribosome = ribosomeData.numOverlappers(geneToUse, false);
		ribosomeCounts.put(geneToUse.toBED(), Double.valueOf(ribosome));
		return ribosome;
	}
	
	/**
	 * @param gene Gene
	 * @return control read count over gene
	 */
	public double getControlCount(Annotation gene) {
		Gene geneToUse = new Gene(gene);
		if(!strandSpecific) {
			geneToUse.setOrientation(Strand.BOTH);
		}
		if(controlCounts.containsKey(geneToUse.toBED())) {
			return controlCounts.get(geneToUse.toBED()).doubleValue();
		}
		int rtrn = controlData.numOverlappers(geneToUse, false);
		controlCounts.put(geneToUse.toBED(), Double.valueOf(rtrn));
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-r", "Ribosome bam file", true);
		p.addStringArg("-m", "Control bam file", true);
		p.addStringArg("-g", "Gene annotation bed file for computing totals", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-gc", "Bed file of genes for computing TE of CDS", true);
		p.addStringArg("-ob", "Output bed file with gene score set to CDS TE", false, null);
		p.addStringArg("-ot", "Output table of CDS TE", false, null);
		p.addDoubleArg("-mtg", "Control global genome total (instead of computing from data)", false, -1);
		p.addDoubleArg("-rtg", "Ribosome global genome total (instead of computing from data)", false, -1);
		p.addDoubleArg("-mte", "Control global exon total (instead of computing from data)", false, -1);
		p.addDoubleArg("-rte", "Ribosome global exon total (instead of computing from data)", false, -1);
		p.addBooleanArg("-ss", "Libraries are strand specific", false, true);
		p.addStringArg("-e", "Experiment ID", true);
		p.parse(args);
		String ribosomeBam = p.getStringArg("-r");
		String controlBam = p.getStringArg("-m");
		String geneAnnotationBed = p.getStringArg("-g");
		String chrSizes = p.getStringArg("-c");
		String geneBed = p.getStringArg("-gc");
		String outputBed = p.getStringArg("-ob");
		String outputTable = p.getStringArg("-ot");
		double ribosomeGenomeTotal = p.getDoubleArg("-rtg");
		double controlGenomeTotal = p.getDoubleArg("-mtg");
		double ribosomeExonTotal = p.getDoubleArg("-rte");
		double controlExonTotal = p.getDoubleArg("-mte");
		boolean strandSpecific = p.getBooleanArg("-ss");
		String experimentId = p.getStringArg("-e");
		
		TranslationalEfficiency te = new TranslationalEfficiency(ribosomeBam, controlBam, geneAnnotationBed, chrSizes, ribosomeGenomeTotal, controlGenomeTotal, ribosomeExonTotal, controlExonTotal, strandSpecific, experimentId);
		
		if(outputBed != null) te.writeCdsTEsToBed(geneBed, outputBed);
		
		if(outputTable != null) te.writeCdsTEsToTable(geneBed, outputTable);
		
		logger.info("");
		logger.info("All done.");
		
	}

	@Override
	public double getScore(Gene region) {
		return getCdsTE(region);
	}

	@Override
	public boolean isSignificant(double score, SignificanceType significanceType) {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getExperimentID() {
		return experimentID;
	}
	
	@Override
	public RegionScore<Gene> createFromConfigFileLine(String line) {
		validateConfigFileLine(line);
		StringParser s = new StringParser();
		s.parse(line);
		String rb = s.asString(0);
		String cb = s.asString(1);
		String geneBed = s.asString(2);
		String chrSizes = s.asString(3);
		double rgt = s.asDouble(4);
		double cgt = s.asDouble(5);
		double ret = s.asDouble(6);
		double cet = s.asDouble(7);
		boolean ss = s.asBoolean(8);
		try {
			return new TranslationalEfficiency(rb, cb, geneBed, chrSizes, rgt, cgt, ret, cet, ss);
		} catch (IOException e) {
			logger.error("Caugh exception:");
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	@Override
	public String getConfigFileLineFormat() {
		return TranslationalEfficiency.class.getSimpleName() + ":\tribosomeBam\tcontrolBam\tgeneAnnotationBed\tchrSizes\tribosomeGenomeTotal\tcontrolGenomeTotal\tribosomeExonTotal\tcontrolExonTotal\tstrandSpecific\texperimentId";
	}

	@SuppressWarnings("unused")
	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 9) {
			logger.error("Field count is not 9: " + line);
			crashWithHelpMessage(line, logger);
		}
		try {
			String rb = s.asString(0);
			String cb = s.asString(1);
			String geneBed = s.asString(2);
			String chrSizes = s.asString(3);
			double rgt = s.asDouble(4);
			double cgt = s.asDouble(5);
			double ret = s.asDouble(6);
			double cet = s.asDouble(7);
			boolean ss = s.asBoolean(8);
		} catch(Exception e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			crashWithHelpMessage(line, logger);
		}
	}


}
