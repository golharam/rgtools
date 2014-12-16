package genome.features;

public interface GenomicFeatureI {
	public String getChromosome();
	public int getStart();
	public void setStart(int start);
	public int getEnd();
	public void setEnd(int end);
	public char getStrand();
	public int getLength();
	public Object getAttribute(String attributeName);
	public String getAttributes();
	public String toString();
}
