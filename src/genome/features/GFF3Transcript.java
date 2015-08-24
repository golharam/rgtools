package genome.features;

import genome.features.GenomicFeature;
import genome.features.GenomicFeatureI;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class GFF3Transcript extends GenomicFeature {

	private String id;
    private List<Exon> exons = new ArrayList<Exon>();
    private List<Exon> cdss = new ArrayList<Exon>();
    private Exon fivePrimeUTR;
    private Exon threePrimeUTR;

    /* This is used by GeneModelBuilder when building transcripts from a set of features */
    public GFF3Transcript(String id) {
    	super(null, Integer.MAX_VALUE, Integer.MIN_VALUE, '.', "gene_id \""+ id + "\";");
    	this.id = id;
    }
    
	public GFF3Transcript(String chromosome, int start, int end, char strand, String attributes) {
		super(chromosome, start, end, strand, attributes);
		this.id = (String)getAttribute("gene_id");
	}

	public void addExon(Exon exon) {
        exons.add(exon);
        setStart(Math.min(exon.getStart(), this.getStart()));
        setEnd(Math.max(exon.getEnd(), this.getEnd()));
    }

    public void addCDS(Exon cds) {
        cdss.add(cds);
        setStart(Math.min(cds.getStart(), this.getStart()));
        setEnd(Math.max(cds.getEnd(), this.getEnd()));
    }

    public void setFivePrimeUTR(Exon exon) {
        fivePrimeUTR = exon;
        setStart(Math.min(exon.getStart(), this.getStart()));
        setEnd(Math.max(exon.getEnd(), this.getEnd()));
    }

    public void setThreePrimeUTR(Exon exon) {
        threePrimeUTR = exon;
        setStart(Math.min(exon.getStart(), this.getStart()));
        setEnd(Math.max(exon.getEnd(), this.getEnd()));
    }

    public int getExonCount() {
    	return exons.size();
    }
    
    /**
     * Create a transcript from its constituent parts.
     *
     * @return
     */
    public GFF3Transcript createTranscript() {
        char strand = '+';
        String name = null;
        String chr;
        
        //FeatureUtils.sortFeatureList(exons);

        // Combine CDS and exons into exons
        for (Exon cds : cdss) {
            //insertCDS(cds);
            Exon exon = findMatchingExon(cds);
            if (exon == null) {
            	// not sure these next two lines are necessary,
            	// since the coding start is automatically set to the size of the exon
                cds.setCodingStart(cds.getStart());
                cds.setCodingEnd(cds.getEnd());
                exons.add(cds);
            } else {
                exon.setCodingStart(cds.getStart());
                exon.setCodingEnd(cds.getEnd());
                exon.setReadingFrame(cds.getReadingFrame());
            }
        }

		for (Exon exon : exons) {
            chr = exon.getChromosome();
            strand = exon.getStrand();
            setStart(Math.min(exon.getStart(), this.getStart()));
            setEnd(Math.max(exon.getEnd(), this.getEnd()));
        }

        sortExons();

        // If 5'UTR is represented by an exon, adjust its start, else add an exon to represent 5'utr
        if (fivePrimeUTR != null) {
            adjustBoundariesByUTR(fivePrimeUTR);
        }

        if (threePrimeUTR != null) {
            adjustBoundariesByUTR(threePrimeUTR);
        }

        return this;
    }
    
    private Exon findMatchingExon(Exon cds) {
        for (Exon exon : exons) {
            if (exon.contains(cds)) {
                return exon;
            }
        }
        return null;
    }

    /**
     * Sort the exon collection, if any, by start position.
     */
    public void sortExons() {
        if (exons != null) {
            Collections.sort(exons, new Comparator<Exon>() {
                public int compare(Exon arg0, Exon arg1) {
                    return arg0.getStart() - arg1.getStart();
                }
            });
        }
    }

    private void adjustBoundariesByUTR(Exon UTR) {
        UTR.setUTR(true);
        addExon(UTR);
        Exon exon = findMatchingExon(UTR);
        if (exon != null) {
            if (exon.getStrand() == '+') {
                exon.setStart(UTR.getEnd());
            } else {
                exon.setEnd(UTR.getStart());
            }
        }
    }
}
