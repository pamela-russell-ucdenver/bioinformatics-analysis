package expression;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import guttmanlab.core.util.StringParser;

/**
 * A record from a Cuffdiff .diff output file
 * @author prussell
 *
 */
public class CuffdiffRecord {
	
	private String id;
	private String geneID;
	private String gene;
	private String locus;
	private String sample1;
	private String sample2;
	private String testStatus;
	private double fpkm1;
	private double fpkm2;
	private double log2fpkmRatio;
	private double testStatistic;
	private double pVal;
	private double qVal;
	private boolean isSignificant;
	
	/**
	 * @param cuffdiffOutputLine Line from .diff output file
	 */
	public CuffdiffRecord(String cuffdiffOutputLine) {
		StringParser s = new StringParser();
		s.parse(cuffdiffOutputLine);
		id = s.asString(0);
		geneID = s.asString(1);
		gene = s.asString(2);
		locus = s.asString(3);
		sample1 = s.asString(4);
		sample2 = s.asString(5);
		testStatus = s.asString(6);
		fpkm1 = s.asDouble(7);
		fpkm2 = s.asDouble(8);
		log2fpkmRatio = s.asDouble(9);
		testStatistic = s.asDouble(10);
		pVal = s.asDouble(11);
		qVal = s.asDouble(12);
		String sig = s.asString(13);
		if(sig.equals("yes")) isSignificant = true;
		if(sig.equals("no")) isSignificant = false;
	}
	
	/**
	 * Load a Cuffdiff output file and get map of ID to full record
	 * @param cuffdiffOutputFile Cuffdiff .diff differential expression output file e.g. isoform_exp.diff or gene_exp.diff
	 * @return Map of ID to full record
	 * @throws IOException
	 */
	public static Map<String, CuffdiffRecord> loadRecordsById(String cuffdiffOutputFile) throws IOException {
		BufferedReader b = new BufferedReader(new FileReader(cuffdiffOutputFile));
		b.readLine(); // discard header line
		Map<String, CuffdiffRecord> rtrn = new HashMap<String, CuffdiffRecord>();
		while(b.ready()) {
			String line = b.readLine();
			CuffdiffRecord record = new CuffdiffRecord(line);
			String recordId = record.getId();
			if(rtrn.containsKey(recordId)) {
				b.close();
				throw new IllegalStateException("Map already has record for ID " + recordId);
			}
			rtrn.put(recordId, record);
		}
		b.close();
		return rtrn;
	}

	public String getId() {
		return id;
	}

	public String getGene() {
		return gene;
	}
	
	public String getGeneID() {
		return geneID;
	}

	public String getLocus() {
		return locus;
	}

	public String getSample1() {
		return sample1;
	}

	public String getSample2() {
		return sample2;
	}

	public String getTestStatus() {
		return testStatus;
	}

	public double getFpkm1() {
		return fpkm1;
	}

	public double getFpkm2() {
		return fpkm2;
	}

	public double getLog2fpkmRatio() {
		return log2fpkmRatio;
	}

	public double getTestStatistic() {
		return testStatistic;
	}

	public double getPval() {
		return pVal;
	}

	public double getQval() {
		return qVal;
	}

	public boolean isSignificant() {
		return isSignificant;
	}


}
