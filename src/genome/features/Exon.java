package genome.features;

import genome.features.GenomicFeature;
import genome.features.GenomicFeatureI;

// Loosy based on IGV Broad Feature Exon

public class Exon extends GenomicFeature {

	/**
     * The index of the exon relative to the start codon.  The exon with the start
     * codon is number "1".
     */
    private int number;
    private int readingFrame = -1;

    /**
     * Coding start position.  This is the leftmost position of the coding region, not neccessarily the 5'utr end
     */
    private int codingStart;
    private int codingEnd;
    private boolean isUTR = false;
    
    private String type;
    
	public Exon(String chromosome, int start, int end, char strand, String type, String attributes) {
		super(chromosome, start, end, strand, attributes);
		this.type = type;
        // By default the entire exon is a coding region
		this.readingFrame = 0;
        this.codingStart = start;
        this.codingEnd = end;
	}

	// should this constructor be allowed since there is no 'type' in a GenomicFeature?
	public Exon(GenomicFeatureI feature) {
		super(feature.getChromosome(), feature.getStart(), feature.getEnd(), feature.getStrand(), feature.getAttributes());
		this.readingFrame = 0;
		this.codingStart = feature.getStart();
		this.codingEnd = feature.getEnd();
	}

	public boolean contains(GenomicFeatureI cds) {
        if (cds == null) {
            return false;
        }
        if (!this.getChromosome().equals(cds.getChromosome()) ||
                this.getStrand() != cds.getStrand()) {
            return false;
        }
        if ((cds.getStart() >= this.getStart()) && (cds.getEnd() <= this.getEnd())) {
            return true;
        } else {
            return false;
        }
	}

	public void setCodingStart(int start) {
		codingStart = start;
	}

	public void setCodingEnd(int end) {
		codingEnd = end;
	}

    public int getReadingFrame() {
        return readingFrame;
    }
    
    public void setReadingFrame(int offset) {
        this.readingFrame = offset;
    }

	public void setUTR(boolean b) {
		isUTR = b;
	}

}
