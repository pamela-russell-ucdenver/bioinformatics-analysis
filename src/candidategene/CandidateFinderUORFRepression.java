package candidategene;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import score.SignificanceType;
import translation.DifferentialTranslationalEfficiency;
import translation.ORFFinder;
import translation.TranslationalEfficiency;
import translation.UpstreamORF;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;

/**
 * Find genes that have a uORF such that
 * uORF and CDS are differentially expressed
 * in opposite directions in the two samples
 * @author prussell
 *
 */
public class CandidateFinderUORFRepression implements CandidateFinder<Gene> {
	
	private DifferentialTranslationalEfficiency diffTE;
	private ORFFinder orfFinder;
	Map<String, FeatureCollection<Gene>> genes;
	private static Logger logger = Logger.getLogger(CandidateFinderUORFRepression.class.getName());

	/**
	 * @param differentialTE Differential translational efficiency object
	 * @param orffinder ORF finder object
	 * @param Bed file of genome annotation
	 * @param chrSizes Chromosome size file
	 * @throws IOException 
	 */
	public CandidateFinderUORFRepression(DifferentialTranslationalEfficiency differentialTE, ORFFinder orffinder, String geneBed, String chrSizes) throws IOException {
		diffTE = differentialTE;
		orfFinder = orffinder;
		genes = BEDFileIO.loadFromFileByReferenceName(geneBed, chrSizes);
	}
	
	/**
	 * @param genomeFasta Fasta file of chromosomes
	 * @param ribosomeBam1 Bam file of ribosome profiling sample 1
	 * @param ribosomeBam 2Bam file of ribosome profiling sample 2
	 * @param controlBam1 Bam file of control sample 1
	 * @param controlBam2 Bam file of control sample 2
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal1 Optional total number of ribosome reads mapped to genome (instead of computing from data), sample 1
	 * @param ribosomeGenomeTotal2 Optional total number of ribosome reads mapped to genome (instead of computing from data), sample 2
	 * @param controlGenomeTotal1 Optional total number of control reads mapped to genome (instead of computing from data), sample 1
	 * @param controlGenomeTotal2 Optional total number of control reads mapped to genome (instead of computing from data), sample 2
	 * @param ribosomeExonTotal1 Optional total number of ribosome reads mapped to exons (instead of computing from data), sample 1
	 * @param ribosomeExonTotal2 Optional total number of ribosome reads mapped to exons (instead of computing from data), sample 2
	 * @param controlExonTotal1 Optional total number of control reads mapped to exons (instead of computing from data), sample 1
	 * @param controlExonTotal2 Optional total number of control reads mapped to exons (instead of computing from data), sample 2
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @param cutoffLog2ratio
	 * @throws IOException
	 */
	public CandidateFinderUORFRepression(String genomeFasta, String ribosomeBam1, String ribosomeBam2, String controlBam1, String controlBam2, String geneBed, String chrSizes, 
			double ribosomeGenomeTotal1, double ribosomeGenomeTotal2, double controlGenomeTotal1, double controlGenomeTotal2, double ribosomeExonTotal1, double ribosomeExonTotal2, 
			double controlExonTotal1, double controlExonTotal2, boolean isStrandSpecific, double cutoffLog2ratio) throws IOException {
		this(DifferentialTranslationalEfficiency.factory(ribosomeBam1, ribosomeBam2, controlBam1, controlBam2, geneBed, chrSizes, ribosomeGenomeTotal1, 
				ribosomeGenomeTotal2, controlGenomeTotal1, controlGenomeTotal2, ribosomeExonTotal1, ribosomeExonTotal2, controlExonTotal1, controlExonTotal2, isStrandSpecific, cutoffLog2ratio),
				new ORFFinder(genomeFasta, chrSizes), geneBed, chrSizes);
	}
	
	/**
	 * Get uORFs that have significantly changed TE in opposite direction of CDS
	 * @param gene Gene
	 * @return Candidate uORFs if they exist and the CDS also has changed TE, or empty collection if none
	 */
	public AnnotationCollection<UpstreamORF> getCandidateUORFs(Gene gene) {
		FeatureCollection<UpstreamORF> rtrn = new FeatureCollection<UpstreamORF>(orfFinder.getCoordSpace());
		double geneDiffTE = diffTE.getScore(gene);
		if(!diffTE.isSignificant(geneDiffTE, SignificanceType.EITHER_SAMPLE_UP)) {
			// Gene TE does not change, therefore there can be no candidate uORFs
			return rtrn;
		}
		boolean cdsIsUp = geneDiffTE > 0;
		AnnotationCollection<UpstreamORF> uorfs = UpstreamORF.findAllUpstreamORFs(orfFinder, gene);
		CloseableIterator<UpstreamORF> iter = uorfs.sortedIterator();
		while(iter.hasNext()) {
			UpstreamORF uorf = iter.next();
			double uorfDiffTE = diffTE.getScore(uorf);
			boolean uorfIsUp = uorfDiffTE > 0;
			if(diffTE.isSignificant(uorfDiffTE, SignificanceType.EITHER_SAMPLE_UP)) {
				if((cdsIsUp && !uorfIsUp) || (!cdsIsUp && uorfIsUp)) {
					logger.debug("Candidate uORF: " + uorf.getCodingRegion().toUCSC());
					rtrn.addAnnotation(uorf);
				}
			}
		}
		iter.close();
		return rtrn;
	}
	
	@Override
	public boolean isCandidate(Gene region) {
		return(getCandidateUORFs(region).getNumAnnotations() > 0);
	}

	@Override
	public String getOutputBedLine(Gene region) {
		AnnotationCollection<UpstreamORF> candidateUORFs = getCandidateUORFs(region);
		String rtrn = "";
		CloseableIterator<UpstreamORF> iter = candidateUORFs.sortedIterator();
		while(iter.hasNext()) {
			UpstreamORF uorf = iter.next();
			rtrn += uorf.getCodingRegion().toBED() + "\n";
		}
		iter.close();
		if(rtrn.equals("")) {
			return null;
		}
		return rtrn;
	}

	
	@Override
	public String getOutputTableLine(Gene region) {
		AnnotationCollection<UpstreamORF> candidateUORFs = getCandidateUORFs(region);
		String rtrn = "";
		CloseableIterator<UpstreamORF> iter = candidateUORFs.sortedIterator();
		while(iter.hasNext()) {
			UpstreamORF uorf = iter.next();
			String line = diffTE.getControl1Name() + "\t";
			line += diffTE.getRibosome1Name() + "\t";
			line += diffTE.getControl2Name() + "\t";
			line += diffTE.getRibosome2Name()	+ "\t";
			line += region.getName() + "\t";
			line += uorf.getCodingRegion().toUCSC() + "\t";
			line += diffTE.getTE1().getControlCount(region.getCodingRegion()) + "\t";
			line += diffTE.getTE1().getRibosomeCount(region.getCodingRegion()) + "\t";
			line += diffTE.getTE1().getScore(region) + "\t";
			line += diffTE.getTE2().getControlCount(region.getCodingRegion()) + "\t";
			line += diffTE.getTE2().getRibosomeCount(region.getCodingRegion()) + "\t";
			line += diffTE.getTE2().getScore(region) + "\t";
			line += diffTE.getScore(region) + "\t";
			line += diffTE.getTE1().getControlCount(uorf.getCodingRegion()) + "\t";
			line += diffTE.getTE1().getRibosomeCount(uorf.getCodingRegion()) + "\t";
			line += diffTE.getTE1().getScore(uorf) + "\t";
			line += diffTE.getTE2().getControlCount(uorf.getCodingRegion()) + "\t";
			line += diffTE.getTE2().getRibosomeCount(uorf.getCodingRegion()) + "\t";
			line += diffTE.getTE2().getScore(uorf) + "\t";
			line += diffTE.getScore(uorf) + "\t";
			line += diffTE.getTE1().getExpressionScanPval(region) + "\t";
			line += diffTE.getTE2().getExpressionScanPval(region) + "\t";
			rtrn += line + "\n";
		}
		iter.close();
		if(rtrn.equals("")) {
			return null;
		}
		return rtrn;
	}

	@Override
	public String getOutputTableHeader() {
		String rtrn = "control_1\t";
		rtrn += "ribosome_1\t";
		rtrn += "control_2\t";
		rtrn += "ribosome_2\t";
		rtrn += "gene\t";
		rtrn += "uORF\t";
		rtrn += "CDS_control_count_1\t";
		rtrn += "CDS_ribosome_count_1\t";
		rtrn += "CDS_TE_1\t";
		rtrn += "CDS_control_count_2\t";
		rtrn += "CDS_ribosome_count_2\t";
		rtrn += "CDS_TE_2\t";
		rtrn += "CDS_log2_TE_ratio\t";
		rtrn += "uORF_control_count_1\t";
		rtrn += "uORF_ribosome_count_1\t";
		rtrn += "uORF_TE_1\t";
		rtrn += "uORF_control_count_2\t";
		rtrn += "uORF_ribosome_count_2\t";
		rtrn += "uORF_TE_2\t";
		rtrn += "uORF_log2_TE_ratio\t";
		rtrn += "expression_scan_pval_control_1\t";
		rtrn += "expression_scan_pval_control_2\t";
		return rtrn;
	}

	/**
	 * Write candidates to a table and bed file
	 * @param outFilePrefix Output file prefix
	 * @throws IOException
	 */
	private void writeResults(String outFilePrefix) throws IOException {
		String outTable = outFilePrefix + ".out";
		String outBed = outFilePrefix + ".bed";
		logger.info("");
		logger.info("Writing candidate uORFs to table " + outTable + " and bed file " + outBed + "...");
		FileWriter wt = new FileWriter(outTable);
		FileWriter wb = new FileWriter(outBed);
		wt.write(getOutputTableHeader() + "\n");
		for(String chr : genes.keySet()) {
			logger.info(chr);
			int numGenes = genes.get(chr).getNumAnnotations();
			CountLogger cl = new CountLogger(numGenes, 10);
			CloseableIterator<Gene> iter = genes.get(chr).sortedIterator();
			while(iter.hasNext()) {
				cl.advance();
				Gene gene = iter.next();
				String line = getOutputTableLine(gene);
				if(line != null) {
					wt.write(line);
					wt.flush();
					wb.write(getOutputBedLine(gene));
					wb.flush();
				}
			}
			iter.close();
		}
		wt.close();
		wb.close();
		logger.info("Done writing file.");
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-gf", "Genome fasta", true);
		p.addStringArg("-rb1", "Ribosome bam file 1", true);
		p.addStringArg("-rb2", "Ribosome bam file 2", true);
		p.addStringArg("-cb1", "Control bam file 1", true);
		p.addStringArg("-cb2", "Control bam file 2", true);
		p.addStringArg("-gb", "Gene bed file", true);
		p.addStringArg("-cs", "Chromosome size file", true);
		p.addDoubleArg("-rgt1", "Ribosome genome total 1", false, -1);
		p.addDoubleArg("-rgt2", "Ribosome genome total 2", false, -1);
		p.addDoubleArg("-cgt1", "Control genome total 1", false, -1);
		p.addDoubleArg("-cgt2", "Control genome total 2", false, -1);
		p.addDoubleArg("-ret1", "Ribosome exon total 1", false, -1);
		p.addDoubleArg("-ret2", "Ribosome exon total 2", false, -1);
		p.addDoubleArg("-cet1", "Control exon total 1", false, -1);
		p.addDoubleArg("-cet2", "Control exon total 2", false, -1);
		p.addBooleanArg("-ss", "Data are strand specific", true);
		p.addDoubleArg("-l2", "Cutoff for absolute value of log2(ratio) for translational efficiency", true);
		p.addStringArg("-ot", "Output table", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addIntArg("-min", "Minimum number of reads mapping to CDS in ribosome and control fraction to compute TE", false, TranslationalEfficiency.TE_MIN_RAW_READS);
		p.parse(args);
		String genomeFasta = p.getStringArg("-gf");
		String ribosomeBam1 = p.getStringArg("-rb1");
		String ribosomeBam2 = p.getStringArg("-rb2");
		String controlBam1 = p.getStringArg("-cb1");
		String controlBam2 = p.getStringArg("-cb2");
		String geneBed = p.getStringArg("-gb");
		String chrSizes = p.getStringArg("-cs");
		double ribosomeGenomeTotal1 = p.getDoubleArg("-rgt1");
		double ribosomeGenomeTotal2 = p.getDoubleArg("-rgt2");
		double controlGenomeTotal1 = p.getDoubleArg("-cgt1");
		double controlGenomeTotal2 = p.getDoubleArg("-cgt2");
		double ribosomeExonTotal1 = p.getDoubleArg("-ret1");
		double ribosomeExonTotal2 = p.getDoubleArg("-ret2");
		double controlExonTotal1 = p.getDoubleArg("-cet1");
		double controlExonTotal2 = p.getDoubleArg("-cet2");
		boolean isStrandSpecific = p.getBooleanArg("-ss");
		double cutoffLog2ratio = p.getDoubleArg("-l2");
		String outFile = p.getStringArg("-ot");
		
		TranslationalEfficiency.TE_MIN_RAW_READS = p.getIntArg("-min");
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
			TranslationalEfficiency.logger.setLevel(Level.DEBUG);
			DifferentialTranslationalEfficiency.logger.setLevel(Level.DEBUG);
		}
		
		logger.info("");
		logger.info("Creating candidate uORF finder object...");
		
		CandidateFinderUORFRepression c = new CandidateFinderUORFRepression(genomeFasta, 
				ribosomeBam1, ribosomeBam2, controlBam1, controlBam2, 
				geneBed, chrSizes, ribosomeGenomeTotal1, ribosomeGenomeTotal2, 
				controlGenomeTotal1, controlGenomeTotal2, ribosomeExonTotal1, 
				ribosomeExonTotal2, controlExonTotal1, controlExonTotal2, 
				isStrandSpecific, cutoffLog2ratio);
		
		logger.info("");
		logger.info("Done creating candidate uORF finder.");
		
		c.writeResults(outFile);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
