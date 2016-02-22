package search;

import guttmanlab.core.pipeline.util.FastqParser;
import guttmanlab.core.pipeline.util.FastqSequence;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import samtools.util.SamtoolsUtils;

/**
 * A query is considered to match a target if they share a perfect kmer match of the specified length
 * Case is ignored
 * Ns are treated as wildcards that everything matches
 * Reverse complement matches are NOT included
 * @author prussell
 *
 */
public class PerfectKmerSearch {
	
	/**
	 * A sequence and position on the sequence
	 * E.g. the start position of a match on a target sequence
	 * @author prussell
	 *
	 */
	private class SequencePos implements Comparable<SequencePos> {
		
		private Sequence seq;
		private int pos;
		
		/**
		 * @param seq The sequence
		 * @param pos The position on the sequence
		 */
		public SequencePos(Sequence seq, int pos) {
			this.seq = seq;
			this.pos = pos;
		}
		
		public boolean equals(Object o) {
			if(!o.getClass().equals(SequencePos.class)) return false;
			SequencePos s = (SequencePos)o;
			return seq.equals(s.getSequence()) && pos == s.getPos();
		}
		
		public String toString() {
			return seq.getName() + ":" + pos;
		}
		
		public int hashCode() {
			String s = toString() + ";" + seq.getSequenceBases();
			return s.hashCode();
		}
		
		public Sequence getSequence() {return seq;}
		public int getPos() {return pos;}

		@Override
		public int compareTo(SequencePos o) {
			int n = seq.getName().compareTo(o.getSequence().getName());
			if(n != 0) return n;
			return pos - o.getPos();
		}
		
	}
	
	/**
	 * A match between kmers of query and target sequences
	 * @author prussell
	 *
	 */
	private class IndividualKmerMatch {
		
		private Sequence query; // Query sequence
		private KmerSubsequence queryKmer; // Kmer sequence and start position of kmer on query sequence
		private SequencePos target; // Start position of match on target sequence
		
		/**
		 * @param queryName Query sequence name
		 * @param queryKmer Kmer from query sequence that matches the target, which stores the start position on the query
		 * @param target Target sequence and start position
		 */
		public IndividualKmerMatch(Sequence query, KmerSubsequence queryKmer, SequencePos target) {
			this.query = query;
			this.queryKmer = queryKmer;
			this.target = target;
		}
		
		/**
		 * @return Query/target pair object for this query and target
		 */
		public QueryTargetPair getQueryTargetPair() {
			return new QueryTargetPair(query, target.getSequence());
		}
		
		public String toString() {
			return query.getName() + ":" + queryKmer.getOrigSeqPos() + "->" + target.toString();
		}
		
		public boolean equals(Object o) {
			if(o.getClass().equals(IndividualKmerMatch.class)) return false;
			IndividualKmerMatch m = (IndividualKmerMatch)o;
			return m.toString().equals(toString());
		}
		
		public int hashCode() {return toString().hashCode();}
		
		public Sequence getQuery() {return query;}
		public Sequence getTarget() {return target.getSequence();}
		public int getQueryStartPos() {return queryKmer.getOrigSeqPos();}		
		public String getQueryName() {return query.getName();}
		public int getTargetStartPos() {return target.getPos();}
		public String getTargetName() {return target.getSequence().getName();}
		
	}
	
	/**
	 * A query and target pair
	 * @author prussell
	 *
	 */
	private class QueryTargetPair {
		
		private Sequence query;
		private Sequence target;
		
		public QueryTargetPair(Sequence query, Sequence target) {
			this.query = query;
			this.target = target;
		}
		
		public Sequence getQuery() {return query;}
		public Sequence getTarget() {return target;}
		
		public String toString() {
			return query.getName() + "->" + target.getName();
		}
		
		public boolean equals(Object o) {
			if(!o.getClass().equals(getClass())) return false;
			QueryTargetPair q = (QueryTargetPair)o;
			return query.equals(q.getQuery()) && target.equals(q.getTarget());
		}
		
		public int hashCode() {
			HashCodeBuilder h = new HashCodeBuilder();
			h.append(query);
			h.append(target);
			return h.toHashCode();
		}
		
	}
	
	/**
	 * A match between a query and target sequence
	 * No gaps are accounted for
	 * @author prussell
	 *
	 */
	private class QueryTargetMatch {
		
		private SequencePos queryMatchStart;
		private SequencePos targetMatchStart;
		private int matchLength;
		
		/**
		 * @param queryMatchStart Query sequence and start position of match
		 * @param targetMatchStart Target sequence and start position of match
		 * @param matchLength Match length
		 */
		public QueryTargetMatch(SequencePos queryMatchStart, SequencePos targetMatchStart, int matchLength) {
			this.queryMatchStart = queryMatchStart;
			this.targetMatchStart = targetMatchStart;
			this.matchLength = matchLength;
		}
		
		/**
		 * Get this match as a SAM record
		 * @return SAM record
		 */
		public SAMRecord toSAMRecord() {
			// Populate cigar
			Cigar cigar = new Cigar();
			int queryStart = queryMatchStart.getPos();
			if(queryStart > 0) {
				// Soft clip beginning of read
				CigarElement s = new CigarElement(queryStart, CigarOperator.S);
				cigar.add(s);
			}
			CigarElement m = new CigarElement(matchLength, CigarOperator.M);
			cigar.add(m);
			int queryLength = queryMatchStart.getSequence().getLength();
			int endSoftClip = queryLength - matchLength - queryStart;
			if(endSoftClip > 0) {
				// Soft clip end of read
				CigarElement s = new CigarElement(endSoftClip, CigarOperator.S);
				cigar.add(s);
			}

			// Set record fields
			SAMRecord rtrn = new SAMRecord(samHeader);
			rtrn.setCigar(cigar);
			rtrn.setAlignmentStart(targetMatchStart.getPos() + 1);
			rtrn.setReadName(queryMatchStart.getSequence().getName());
			rtrn.setReadPairedFlag(false);
			rtrn.setReferenceName(targetMatchStart.getSequence().getName());
			rtrn.setMappingQuality(255); // mapping quality unknown
			
			return rtrn;
		}
		
	}
	
	/**
	 * Organize kmer matches by query and target, then get the "first" match for each query/target pair
	 * @param kmerMatches Collection of kmer matches with different queries and targets allowed
	 * @return Map of query/target pair to "first" kmer match from the set wrt query and target coordinates
	 */
	private Map<QueryTargetPair, QueryTargetMatch> firstKmerMatchEachQueryTargetPair(Collection<IndividualKmerMatch> kmerMatches) {
		
		Map<QueryTargetPair, Collection<IndividualKmerMatch>> matchesByPair = new HashMap<QueryTargetPair, Collection<IndividualKmerMatch>>();
		for(IndividualKmerMatch match : kmerMatches) {
			QueryTargetPair qtp = match.getQueryTargetPair();
			if(!matchesByPair.containsKey(qtp)) {
				matchesByPair.put(qtp, new ArrayList<IndividualKmerMatch>());
			}
			matchesByPair.get(qtp).add(match);
		}
		
		Map<QueryTargetPair, QueryTargetMatch> rtrn = new HashMap<QueryTargetPair, QueryTargetMatch>();
		for(QueryTargetPair qtp : matchesByPair.keySet()) {
			rtrn.put(qtp, firstKmerMatch(matchesByPair.get(qtp)));
		}
		
		return rtrn;
		
	}
	
	/**
	 * Get the first kmer match of this query to each of its targets
	 * @param query Query sequence
	 * @return The first match to each target
	 */
	private Collection<QueryTargetMatch> firstKmerMatchEachTarget(Sequence query) {
		Map<QueryTargetPair, QueryTargetMatch> matches = firstKmerMatchEachQueryTargetPair(getIndividualKmerMatches(query));
		return matches.values();
	}
	
	/**
	 * Get the first kmer match of this query to each of its targets as SAM records
	 * @param record Query read
	 * @return The first match to each target as SAM records
	 */	
	private Collection<SAMRecord> samRecordFirstKmerMatchEachTarget(FastqSequence record) {
		
		// Check how many N's are in the read
		// TODO this adds a lot of time
		String querySeq = record.getSequence();
		int numNs = 0;
		for(int i = 0; i < querySeq.length(); i++) {
			if(Character.toUpperCase(querySeq.charAt(i)) == 'N') {numNs++;}
		}
		if((double) numNs / (double) querySeq.length() > MAX_PCT_N) {
			throw new TooManyNsException("Read has >" + MAX_PCT_N + " Ns:\t" + record.getName() + "\t" + querySeq);
		}
		
		return samRecordFirstKmerMatchEachTarget(new Sequence(record.getName(), record.getSequence()));
	}
	
	/**
	 * Get the first kmer match of this query to each of its targets as SAM records
	 * @param query Query sequence
	 * @return The first match to each target as SAM records
	 */
	private Collection<SAMRecord> samRecordFirstKmerMatchEachTarget(Sequence query) {
		Collection<QueryTargetMatch> matches = firstKmerMatchEachTarget(query);
		Collection<SAMRecord> rtrn = new ArrayList<SAMRecord>();
		for(QueryTargetMatch match : matches) {
			rtrn.add(match.toSAMRecord());
		}
		return rtrn;
	}
	
	/**
	 * Iterate through fastq file and for each query and target, write first kmer match to a bam file
	 * @param queryFastq Query fastq file
	 * @param outputBam Bam file to write
	 * @throws IOException
	 */
	private void writeFirstKmerMatchEachTarget(String queryFastq, String outputBam) throws IOException {
		
		BAMFileWriter writer = new BAMFileWriter(new File(outputBam));
		writer.setSortOrder(SAMFileHeader.SortOrder.unsorted, false);
		writer.setHeader(samHeader);
		FastqParser reader = new FastqParser();
		reader.start(new File(queryFastq));
		int numTooShort = 0;
		int numIllegalChar = 0;
		int numUniquelyMapped = 0;
		int numMultiMapped = 0;
		int numUnmapped = 0;
		int numDone = 0;
		int numTooManyNs = 0;
		while(reader.hasNext()) {
			FastqSequence query = reader.next();
			numDone++;
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " reads");
			}
			try {
				Collection<SAMRecord> alignments = samRecordFirstKmerMatchEachTarget(query);
				if(alignments.size() == 0) numUnmapped++;
				if(alignments.size() == 1) numUniquelyMapped++;
				if(alignments.size() > 1) numMultiMapped++;
				for(SAMRecord alignment : alignments) {
					writer.addAlignment(alignment);
				}
			} catch(SequenceTooShortException e) {
				numTooShort++;
			} catch(IllegalCharacterException e) {
				numIllegalChar++;
			} catch(TooManyNsException e) {
				numTooManyNs++;
			}
		}
		logger.info("");
		logger.info("RESULTS");
		logger.info("Reads mapped uniquely:\t" + numUniquelyMapped);
		logger.info("Reads mapped to multiple targets:\t" + numMultiMapped);
		logger.info("Reads unmapped:\t" + numUnmapped);
		if(numTooShort > 0) {
			logger.warn("Reads skipped because they were shorter than " + k + ":\t" + numTooShort);
		}
		if(numIllegalChar > 0) {
			logger.warn("Reads skipped because they contain an illegal character " + k + ":\t" + numIllegalChar);
		}
		if(numTooManyNs > 0) {
			logger.warn("Reads skipped because they contain > " + MAX_PCT_N + " N's:\t" + numTooManyNs);
		}
		logger.info("");
		reader.close();
		writer.close();
		
	}
	
	/**
	 * Get the "first" kmer match between a query and a target, out of a set of multiple kmer matches between these sequences
	 * @param kmerMatches Set of kmer matches all with same query and target
	 * @return The "first" match, i.e. smallest position on query and target represented in the set of kmer matches
	 */
	private QueryTargetMatch firstKmerMatch(Collection<IndividualKmerMatch> kmerMatches) {
		Iterator<IndividualKmerMatch> iter = kmerMatches.iterator();
		if(!iter.hasNext()) {
			throw new IllegalArgumentException("Iterator empty");
		}
		IndividualKmerMatch match = iter.next();
		Sequence query = match.getQuery();
		Sequence target = match.getTarget();
		String queryName = match.getQueryName();
		String targetName = match.getTargetName();
		int queryStart = match.getQueryStartPos();
		int targetStart = match.getTargetStartPos();
		while(iter.hasNext()) {
			IndividualKmerMatch next = iter.next();
			String nextQueryName = next.getQueryName();
			if(!nextQueryName.equals(queryName)) {
				throw new IllegalArgumentException("Must have only one query: " + queryName + ", " + nextQueryName);
			}
			String nextTargetName = next.getTargetName();
			if(!nextTargetName.equals(targetName)) {
				throw new IllegalArgumentException("Must have only one target: " + targetName + ", " + nextTargetName);
			}
			int nextQueryStart = next.getQueryStartPos();
			if(nextQueryStart < queryStart) {queryStart = nextQueryStart;}
			int nextTargetStart = next.getTargetStartPos();
			if(nextTargetStart < targetStart) {targetStart = nextTargetStart;}
		}
		SequencePos queryMatchPos = new SequencePos(query, queryStart);
		SequencePos targetMatchPos = new SequencePos(target, targetStart);
		return new QueryTargetMatch(queryMatchPos, targetMatchPos, k);
	}
	
	/**
	 * A kmer sequence and the start position of the original sequence it came from
	 * @author prussell
	 *
	 */
	private class KmerSubsequence {
		
		private String seq; // The kmer sequence
		private int origSeqPos; // Start position of the kmer on the original sequence
		
		/**
		 * @param seq Kmer sequence
		 * @param origSeqPos Start position on original sequence
		 */
		public KmerSubsequence(String seq, int origSeqPos) {
			if(seq == null) throw new IllegalArgumentException("Sequence is null");
			this.seq = seq;
			this.origSeqPos = origSeqPos;
		}
		
		public boolean equals(Object o) {
			if(!o.getClass().equals(KmerSubsequence.class)) return false;
			KmerSubsequence k = (KmerSubsequence)o;
			return k.toString().equals(toString());
		}
		
		public int hashCode() {
			return toString().hashCode();
		}
		
		public String toString() {
			return origSeqPos + ":" + seq;
		}
		
		public String getSeq() {return seq;}
		public int getOrigSeqPos() {return origSeqPos;}
		
	}
	
	private int k; // Kmer length
	private Map<String, Collection<SequencePos>> targetKmers; // Key is kmer; value is collection of sequences with kmer and the match position
	private static Logger logger = Logger.getLogger(PerfectKmerSearch.class.getName());
	private SAMFileHeader samHeader; // SAM header for target sequences
	private static double MAX_PCT_N = 0.05;
	
	/**
	 * The legal characters converted to upper case, not including N
	 */
	public static final char[] alphabet = {'A', 'C', 'G', 'T'};
	
	/**
	 * @param k Length of kmers to match
	 * @param fasta Fasta file of target sequences
	 */
	public PerfectKmerSearch(int k, String fasta) {
		this.k = k;
		createIndex(fasta);
		samHeader = SamtoolsUtils.createSamHeader(fasta);
	}

	
	/**
	 * Check if the char (converted to upper case) is in the alphabet or is N
	 * @param c Char to check
	 * @return True iff upper case version of char is in alphabet or is N
	 */
	private static boolean charIsLegal(char c) {
		char cu = Character.toUpperCase(c);
		if(cu == 'N') return true;
		for(char a : alphabet) {
			if(cu == a) return true;
		}
		return false;
	}
	
	@SuppressWarnings("serial")
	private class SequenceTooShortException extends RuntimeException {
		public SequenceTooShortException(String message) {
			super(message);
		}
	}
	
	@SuppressWarnings("serial")
	private class TooManyNsException extends RuntimeException {
		public TooManyNsException(String message) {
			super(message);
		}
	}
	
	@SuppressWarnings("serial")
	private class IllegalCharacterException extends RuntimeException {
		public IllegalCharacterException(String message) {
			super(message);
		}
	}
	
	/**
	 * Check that a sequence is valid
	 * @param seq Sequence
	 */
	private void validateSequence(Sequence seq) {
		String bases = seq.getSequenceBases();
		int len = bases.length();
		if(len < k) {
			throw new SequenceTooShortException("Sequence length <" + k + ": " + seq.getName());
		}
		for(int i = 0; i < len; i++) {
			char c = bases.charAt(i);
			if(!charIsLegal(c)) {
				throw new IllegalCharacterException("Illegal char in sequence " + seq.getName() + ": " + c);
			}
		}
	}
	
	/**
	 * Store kmers and their matches to target sequences
	 * @param fasta Fasta file of target sequences
	 */
	private void createIndex(String fasta) {
		logger.info("");
		logger.info("Creating index for target fasta " + fasta + "...");
		targetKmers = new HashMap<String, Collection<SequencePos>>();
		Collection<Sequence> targets = new FastaFileIOImpl().readFromFile(fasta);
		int numSkipped = 0;
		for(Sequence target : targets) {
			logger.debug("");
			logger.debug("TARGET\t" + target.getName());
			logger.debug("TARGET_SEQ\t" + target.getSequenceBases());
			try {
				validateSequence(target);
			} catch(SequenceTooShortException e) {
				logger.warn("Caught exception, skipping target sequence:\t" + e.getMessage());
				numSkipped++;
				continue;
			}
			for(KmerSubsequence kmer : getKmers(target.getSequenceBases())) {
				String kmerSeq = kmer.getSeq();
				if(!targetKmers.containsKey(kmerSeq)) {
					targetKmers.put(kmerSeq, new TreeSet<SequencePos>());
				}
				targetKmers.get(kmerSeq).add(new SequencePos(target, kmer.getOrigSeqPos()));
				logger.debug("ADDED\t" + kmer + "\t" + target.getName());
			}
		}
		if(numSkipped > 0) {
			logger.warn("");
			logger.warn("Skipped " + numSkipped + " target sequences that did not validate");
			logger.warn("");
		}
	}
	
	/**
	 * Create multiple versions of the sequence for every possible value of N's
	 * @param sequence Sequence to expand
	 * @return Collection of versions with all possible values of N's taken from the alphabet
	 */
	private Collection<String> expandNs(String sequence) {
		for(int i = 0; i < sequence.length(); i++) {
			if(Character.toUpperCase(sequence.charAt(i)) == 'N') {
				Collection<String> expanded = new HashSet<String>();
				StringBuilder sb = new StringBuilder(sequence);
				for(char a : alphabet) {
					sb.replace(i, i+1, Character.toString(a));
					expanded.add(sb.toString());
				}
				return expandNs(expanded);
			}
		}
		Collection<String> rtrn = new HashSet<String>();
		rtrn.add(sequence);
		logger.debug("EXPANDED_Ns\t" + sequence);
		return rtrn;
	}
		
	/**
	 * Create multiple versions of each sequence for every possible value of N's
	 * @param sequences Sequences to expand
	 * @return Collection of versions with all possible values of N's taken from the alphabet
	 */
	private Collection<String> expandNs(Collection<String> sequences) {
		Collection<String> rtrn = new HashSet<String>();
		for(String seq : sequences) {
			rtrn.addAll(expandNs(seq));
		}
		return rtrn;
	}
	
	/**
	 * Get all kmer substrings of the sequence, converted to upper case
	 * N's are expanded so the method returns multiple versions of each kmer where there is an N
	 * Clients should call validateSequence() before calling this method to get meaningful error messages
	 * @param sequence Target sequence to break into kmers
	 * @return Set of kmers converted to upper case with Ns expanded to all possible values
	 * @throws IllegalArgumentException if sequence is shorter than kmer length
	 */
	private Collection<KmerSubsequence> getKmers(String sequence) {
		int len = sequence.length();
		Collection<KmerSubsequence> rtrn = new ArrayList<KmerSubsequence>();
		StringBuilder builder = new StringBuilder(sequence.substring(0, k));
		logger.debug("");
		logger.debug("EXPANDING_Ns for kmer\t" + builder.toString());
		Collection<String> expandedNs = expandNs(builder.toString().toUpperCase());
		for(String s : expandedNs) {
			rtrn.add(new KmerSubsequence(s, 0));
		}
		for(int i = k; i < len; i++) {
			builder.deleteCharAt(0);
			builder.append(sequence.charAt(i));
			logger.debug("");
			logger.debug("EXPANDING_Ns for kmer\t" + builder.toString());
			Collection<String> expandedNs2 = expandNs(expandNs(builder.toString().toUpperCase()));
			for(String s : expandedNs2) {
				rtrn.add(new KmerSubsequence(s, i - k + 1));
			}
		}
		if(logger.getLevel().equals(Level.DEBUG)) {
			StringBuilder s = new StringBuilder();
			for(KmerSubsequence kmer : rtrn) s.append(kmer.getSeq() + ";");
			logger.debug("KMERS\t" + sequence + "\t" + s.toString());
		}
		return rtrn;
	}
	
	/**
	 * Get all kmer matches of this query to the stored targets, based on kmer matches
	 * There can be multiple matches to a given target
	 * @param query Query sequence
	 * @return Set of perfect kmer matches
	 */
	private Collection<IndividualKmerMatch> getIndividualKmerMatches(Sequence query) {
		validateSequence(query);
		Collection<KmerSubsequence> queryKmers = getKmers(query.getSequenceBases());
		Collection<IndividualKmerMatch> rtrn = new HashSet<IndividualKmerMatch>();
		for(KmerSubsequence queryKmer : queryKmers) {
			String kmerSeq = queryKmer.getSeq();
			if(!targetKmers.containsKey(kmerSeq)) {continue;}
			for(SequencePos sp : targetKmers.get(kmerSeq)) {
				IndividualKmerMatch match = new IndividualKmerMatch(query, queryKmer, sp);
				rtrn.add(match);
			}
		}
		return rtrn;
	}
	
	/**
	 * Write kmer index out to a file
	 * @param outFile File to write
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeKmerIndex(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(String kmer : targetKmers.keySet()) {
			Iterator<SequencePos> iter = targetKmers.get(kmer).iterator();
			String targets = iter.next().toString();
			while(iter.hasNext()) {
				targets += "," + iter.next().toString();
			}
			w.write(kmer + "\t" + targets + "\n");
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-fa", "Reference fasta", true);
		p.addStringArg("-fq", "Query fastq", true);
		p.addStringArg("-b", "Output bam", true);
		p.addIntArg("-k", "Kmer length", true);
		p.addDoubleArg("-mn", "Max proportion of N's in query sequence", false, MAX_PCT_N);
		p.parse(args);
		String fasta = p.getStringArg("-fa");
		String fastq = p.getStringArg("-fq");
		String bam = p.getStringArg("-b");
		int k = p.getIntArg("-k");
		MAX_PCT_N = p.getDoubleArg("-mn");
		
		logger.setLevel(Level.INFO);
		PerfectKmerSearch pks = new PerfectKmerSearch(k, fasta);
		pks.writeFirstKmerMatchEachTarget(fastq, bam);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
	

}
