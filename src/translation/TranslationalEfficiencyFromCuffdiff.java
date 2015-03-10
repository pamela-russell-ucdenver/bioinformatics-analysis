package translation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import score.AbstractRegionScore;
import score.RegionScore;
import score.SignificanceType;
import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

public class TranslationalEfficiencyFromCuffdiff extends AbstractRegionScore<Gene> {
	
	private double ribosomeGlobalExonTotal;
	private double controlGlobalExonTotal;
	private double normalizationFactor; // Ribosome total / control total
	public static int TE_MIN_RAW_READS = 10; // Minimum number of reads in ribosome and control fraction to compute the TE of a CDS
	private Map<String, Double> ribosomeCounts;
	private Map<String, Double> controlCounts;
	private String ribosomeName;
	private String controlName;
	private String experimentID;
	
	public static Logger logger = Logger.getLogger(TranslationalEfficiencyFromCuffdiff.class.getName());
	
	public TranslationalEfficiencyFromCuffdiff() {}
	
	public TranslationalEfficiencyFromCuffdiff(String cuffdiffCountTrackingFile, String ribosomeSampleName, String controlSampleName, String experimentId) throws IOException {
		ribosomeName = ribosomeSampleName;
		controlName = controlSampleName;
		experimentID = experimentId;
		initializeCountsAndTotals(cuffdiffCountTrackingFile);
		normalizationFactor = ribosomeGlobalExonTotal / controlGlobalExonTotal;
	}
	
	private void initializeCountsAndTotals(String cuffdiffCountTrackingFile) throws IOException {
		int ribosomeCountCol = getCountColNum(cuffdiffCountTrackingFile, ribosomeName);
		int controlCountCol = getCountColNum(cuffdiffCountTrackingFile, controlName);
		BufferedReader r = new BufferedReader(new FileReader(cuffdiffCountTrackingFile));
		r.readLine();
		ribosomeGlobalExonTotal = 0;
		controlGlobalExonTotal = 0;
		ribosomeCounts = new HashMap<String, Double>();
		controlCounts = new HashMap<String, Double>();
		StringParser s = new StringParser();
		while(r.ready()) {
			s.parse(r.readLine());
			double ribosomeCount = s.asDouble(ribosomeCountCol);
			double controlCount = s.asDouble(controlCountCol);
			ribosomeGlobalExonTotal += ribosomeCount;
			controlGlobalExonTotal += controlCount;
			String name = s.asString(0);
			ribosomeCounts.put(name, Double.valueOf(ribosomeCount));
			controlCounts.put(name, Double.valueOf(controlCount));
		}
		r.close();
	}
	
	private static int getCountColNum(String cuffdiffCountTrackingFile, String sampleName) throws IOException {
		String colHeader = sampleName + "_count";
		BufferedReader r = new BufferedReader(new FileReader(cuffdiffCountTrackingFile));
		StringParser s = new StringParser();
		String line = r.readLine();
		s.parse(line);
		for(int i = 0; i < s.getFieldCount(); i++) {
			if(s.asString(i).equals(colHeader)) {
				r.close();
				return i;
			}
		}
		r.close();
		throw new IllegalArgumentException("Column header " + colHeader + " not found in line " + line);
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
	 * Get translational efficiency score for a region
	 * @param region The region
	 * @return The TE score, or NaN if either region does not meet min cutoff
	 */
	public double getTE(Annotation region) {
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
	 * Write TE scores for the CDS of each gene to a table
	 * @param geneBed Bed file of genes to analyze
	 * @param outputTable Table file to write
	 * @throws IOException 
	 */
	public void writeTEsToTable(String geneBed, String chrSizeFile, String outputTable) throws IOException {
		logger.info("");
		FileWriter w = new FileWriter(outputTable);
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		int numGenes = genes.getNumAnnotations();
		logger.info("Writing table with TE scores for " + numGenes + " genes in " + geneBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		String header = "gene\t";
		header += "global_exon_count_control\t";
		header += "global_exon_count_ribosome\t";
		header += "gene_size\t";
		header += "gene_count_control\t";
		header += "gene_count_ribosome\t";
		header += "TE_normalization_factor\t";
		header += "TE_score\t";
		w.write(header + "\n");
		while(iter.hasNext()) {
			Gene gene = iter.next();
			try {
				double mCountGene = getControlCount(gene);
				double rCountGene = getRibosomeCount(gene);
				double te = getTE(gene);
				int geneSize = gene.size();
				String line = gene.getName() + "\t";
				line += controlGlobalExonTotal + "\t";
				line += ribosomeGlobalExonTotal + "\t";
				line += geneSize + "\t";
				line += mCountGene + "\t";
				line += rCountGene + "\t";
				line += normalizationFactor + "\t";
				line += te + "\t";
				w.write(line + "\n");
			} catch (IllegalArgumentException e) {
				logger.warn("Gene not found in cuffdiff table: " + gene.getName() + ". Skipping.");
			}
		}
		w.close();
		iter.close();
		logger.info("Done writing to " + outputTable + ".");
	}
	
		
	/**
	 * @param gene Gene
	 * @return Ribosome read count over gene
	 */
	public double getRibosomeCount(Annotation gene) {
		if(ribosomeCounts.containsKey(gene.getName())) {
			return ribosomeCounts.get(gene.getName()).doubleValue();
		}
		throw new IllegalArgumentException("Gene not found in cuffdiff table: " + gene.getName());
	}
	
	/**
	 * @param gene Gene
	 * @return control read count over gene
	 */
	public double getControlCount(Annotation gene) {
		if(controlCounts.containsKey(gene.getName())) {
			return controlCounts.get(gene.getName()).doubleValue();
		}
		throw new IllegalArgumentException("Gene not found in cuffdiff table: " + gene.getName());
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-cc", "Cuffdiff count tracking file", true);
		p.addStringArg("-rn", "Ribosome sample name", true);
		p.addStringArg("-cn", "Control sample name", true);
		p.addStringArg("-gc", "Bed file of genes with names matching Cuffdiff table", true);
		p.addStringArg("-ot", "Output table of TE scores", false, null);
		p.addStringArg("-e", "Experiment ID", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.parse(args);
		String geneBed = p.getStringArg("-gc");
		String outputTable = p.getStringArg("-ot");
		String experimentId = p.getStringArg("-e");
		String chrSizeFile = p.getStringArg("-c");
		String cuffdiffCountTrackingFile = p.getStringArg("-cc");
		String ribosomeSampleName = p.getStringArg("-rn");
		String controlSampleName = p.getStringArg("-cn");
		
		TranslationalEfficiencyFromCuffdiff te = new TranslationalEfficiencyFromCuffdiff(cuffdiffCountTrackingFile, ribosomeSampleName, controlSampleName, experimentId);
		
		if(outputTable != null) te.writeTEsToTable(geneBed, chrSizeFile, outputTable);
		
		logger.info("");
		logger.info("All done.");
		
	}

	@Override
	public double getScore(Gene region) {
		return getTE(region);
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
		String cuffdiff = s.asString(0);
		String rName = s.asString(1);
		String cName = s.asString(2);
		String expId = s.asString(3);
		try {
			return new TranslationalEfficiencyFromCuffdiff(cuffdiff, rName, cName, expId);
		} catch(IOException e) {
			logger.info("Caught exception: " + e.getMessage());
			System.exit(-1);
			return null;
		}
	}

	@Override
	public String getConfigFileLineFormat() {
		return TranslationalEfficiencyFromCuffdiff.class.getSimpleName() + ":\tcuffdiffCountTrackingFile\tribosomeSampleName\tcontrolSampleName\texperimentId";
	}

	@SuppressWarnings("unused")
	@Override
	public void validateConfigFileLine(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		if(s.getFieldCount() != 4) {
			logger.error("Field count is not 4: " + line);
			crashWithHelpMessage(line, logger);
		}
		try {
			String cuffdiff = s.asString(0);
			String rName = s.asString(1);
			String cName = s.asString(2);
			String expId = s.asString(3);
		} catch(Exception e) {
			logger.error("Caught exception:");
			e.printStackTrace();
			crashWithHelpMessage(line, logger);
		}
	}


}
