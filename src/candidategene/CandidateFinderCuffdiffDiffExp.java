package candidategene;

import java.io.FileWriter;
import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;

import org.apache.log4j.Logger;

import score.SignificanceType;
import expression.CuffdiffRecord;
import expression.DifferentialExpressionCuffdiff;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.util.CommandLineParser;

/**
 * Find genes that are called differentially expressed by Cuffdiff
 * @author prussell
 */
public class CandidateFinderCuffdiffDiffExp implements CandidateFinder<Gene> {
	
	private DifferentialExpressionCuffdiff diffexp;
	private static Logger logger = Logger.getLogger(CandidateFinderCuffdiffDiffExp.class.getName());
	private AnnotationCollection<Gene> genes;
	
	/**
	 * @param cuffdiffOutputIsoformExpDiff Cuffdiff isoform_exp.diff file
	 * @param genesBed Bed file of genes
	 * @param chrSizes Chromosome size file
	 * @throws IOException
	 */
	public CandidateFinderCuffdiffDiffExp(String cuffdiffOutputIsoformExpDiff, String genesBed, String chrSizes) throws IOException {
		logger.info("");
		logger.info("Instantiating candidate finder with cuffdiff file " + cuffdiffOutputIsoformExpDiff + " and genes in " + genesBed + ".");
		logger.info("Cuffdiff Q value cutoff is " + DifferentialExpressionCuffdiff.QVAL_CUTOFF + ".");
		diffexp = new DifferentialExpressionCuffdiff(cuffdiffOutputIsoformExpDiff);
		genes = BEDFileIO.loadFromFile(genesBed, chrSizes);
	}
	
	/**
	 * Write table of cuffdiff information for the genes
	 * @param outTable Output table
	 * @throws IOException
	 */
	private void writeTable(String outTable) throws IOException {
		logger.info("");
		logger.info("Writing table to " + outTable + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		FileWriter w = new FileWriter(outTable);
		w.write(getOutputTableHeader());
		while(iter.hasNext()) {
			Gene gene = iter.next();
			try {
				w.write(getOutputTableLine(gene));
			} catch(IllegalArgumentException e) {
				logger.warn("Caught exception; skipping gene " + gene.getName() + ". " + e.getMessage());
			}
		}
		w.close();
		iter.close();
		logger.info("Done writing table.");
	}
	
	/**
	 * Write bed file of significant genes with cuffdiff scores
	 * @param outBed Output bed file
	 * @throws IOException
	 */
	private void writeBedFileOfSignificantGenes(String outBed) throws IOException {
		logger.info("");
		logger.info("Writing bed file of significant genes with scores to file " + outBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		FileWriter w = new FileWriter(outBed);
		while(iter.hasNext()) {
			Gene gene = iter.next();
			try {
				if(isCandidate(gene)) {
					w.write(getOutputBedLine(gene));
				}
			} catch(IllegalArgumentException e) {
				logger.warn("Caught exception; skipping gene " + gene.getName() + ". " + e.getMessage());
			}
		}
		w.close();
		iter.close();
		logger.info("Done writing bed file.");
	}
	
	@Override
	public boolean isCandidate(Gene region) {
		return diffexp.isSignificant(diffexp.getScore(region), SignificanceType.EITHER_SAMPLE_UP);
	}

	@Override
	public String getOutputTableLine(Gene region) {
		CuffdiffRecord record = diffexp.getRecord(region);
		String rtrn = record.getSample1() + "\t";
		rtrn += record.getSample2() + "\t";
		rtrn += record.getGene() + "\t";
		rtrn += record.getLocus() + "\t";
		rtrn += record.getFpkm1() + "\t";
		rtrn += record.getFpkm2() + "\t";
		rtrn += record.getLog2fpkmRatio() + "\t";
		rtrn += record.getPval() + "\t";
		rtrn += record.getQval() + "\t";
		rtrn += isCandidate(region) + "\t";
		rtrn += "\n";
		return rtrn;
	}

	@Override
	public String getOutputBedLine(Gene region) {
		double score = diffexp.getScore(region);
		return region.toBED(score) + "\n";
	}

	@Override
	public String getOutputTableHeader() {
		String rtrn = "sample1\t";
		rtrn += "sample2\t";
		rtrn += "gene_name\t";
		rtrn += "locus\t";
		rtrn += "fpkm1\t";
		rtrn += "fpkm2\t";
		rtrn += "log2_fpkm_ratio\t";
		rtrn += "pval\t";
		rtrn += "qval\t";
		rtrn += "is_candidate\t";
		rtrn += "\n";
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-ic", "Cuffdiff isoform_exp.diff file", true);
		p.addDoubleArg("-q", "Q value cutoff", false, DifferentialExpressionCuffdiff.QVAL_CUTOFF);
		p.addStringArg("-ot", "Output table of cuffdiff info", false, null);
		p.addStringArg("-ob", "Output bed file of significant genes", false, null);
		p.addStringArg("-g", "Genes bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.parse(args);
		String cuffdiffFile = p.getStringArg("-ic");
		double qvalCutoff = p.getDoubleArg("-q");
		String outTable = p.getStringArg("-ot");
		String outBed = p.getStringArg("-ob");
		String geneBed = p.getStringArg("-g");
		String chrSizes = p.getStringArg("-c");
		
		DifferentialExpressionCuffdiff.QVAL_CUTOFF = qvalCutoff;
		
		CandidateFinderCuffdiffDiffExp cf = new CandidateFinderCuffdiffDiffExp(cuffdiffFile, geneBed, chrSizes);
		
		if(outTable != null) {
			cf.writeTable(outTable);
		}
		
		if(outBed != null) {
			cf.writeBedFileOfSignificantGenes(outBed);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
