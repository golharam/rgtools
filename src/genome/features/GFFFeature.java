package genome.features;

import java.util.Comparator;

/* A GFF/GTF Feature is simply a GenomicFeature with attributes */
public class GFFFeature extends GenomicFeature implements Comparator<GFFFeature> {

	String source;
	double score;
	char frame;
	String type;
	
	public GFFFeature(String chromosome, int start, int end, char strand, String type, String source, double score, char frame, String attributes) {
		super(chromosome, start, end, strand, attributes);
		this.source = source;
		this.score = score;
		this.frame = frame;
		this.type = type;
	}
	
	public GFFFeature(Exon exon) {
		super(exon.getChromosome(), exon.getStart(), exon.getEnd(), exon.getStrand(), exon.getAttributes());
		this.source = "";
		this.score = 0.0;
		this.frame = 0;
		this.type = "exon";
	}
	
	public String getType() {
		return type;
	}
	
	public String getSource() {
		return source;
	}
	
	public double getScore() {
		return score;
	}
	
	public String toString() {
		return getChromosome() + "\t" + source + "\t" + type + "\t" + getStart() + "\t" + getEnd() + "\t" + score + "\t" + getStrand() + "\t" + frame + "\t" + getAttributes();
	}
	

	@Override
	public int compare(GFFFeature a, GFFFeature b) {
		char a_chr = a.getChromosome().charAt(3);
		char b_chr = b.getChromosome().charAt(3);
		if (a_chr < b_chr) {
			return -1;
		}
		if (a_chr > b_chr) {
			return 1;
		}
		if (a_chr == b_chr) {
			return a.getStart() - b.getStart();
		} else {
			throw new IndexOutOfBoundsException(
					"Something went wrong comparing two features:\n"
							+ a.toString() + "\n" + b.toString());
		}
	}
}
