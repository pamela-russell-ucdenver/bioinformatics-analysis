package sam;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

public class SAMRecordUtils {
	
	/**
	 * Convert reference position to position within a mapped read
	 * Returns position on the read sequence itself regardless of mapping orientation
	 * @param record Mapped read
	 * @param zeroBasedRefPos Zero-based reference position
	 * @return Zero-based position on read, or -1 if no position
	 */
	public static int getReadPositionAtReferencePosition(SAMRecord record, int zeroBasedRefPos) {
		int oneBasedRtrn = -1;
		int oneBasedRefPos = zeroBasedRefPos + 1;
		/*
		 *  We can only identify the reference position if it is in one of the aligned blocks,
		 *  so ignore parts of the read that did not align
		 */
		for(AlignmentBlock block : record.getAlignmentBlocks()) {
			int refStart = block.getReferenceStart();
			int refEnd = refStart + block.getLength() - 1;
			if(!(oneBasedRefPos >= refStart && oneBasedRefPos <= refEnd)) {
				continue;
			}
			int readStart = block.getReadStart();
			oneBasedRtrn = readStart + oneBasedRefPos - refStart;
		}
		if(oneBasedRtrn < 0) {
			// We didn't find the position in a mapped block
			return oneBasedRtrn;
		}
		// Convert to read coordinates if mapped on minus strand
		if(record.getReadNegativeStrandFlag()) {
			oneBasedRtrn = record.getReadLength() - oneBasedRtrn + 1;
		}
		return oneBasedRtrn - 1;
	}
	
}
