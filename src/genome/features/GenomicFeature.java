package genome.features;

import java.util.HashMap;

public class GenomicFeature implements GenomicFeatureI {
	private String chromosome;
	private int start;
	private int end;
	private char strand;
	private String attributes;
	private HashMap<String,Object> attributeHashMap = new HashMap<String,Object>();

	public GenomicFeature(String chromosome, int start, int end, char strand, String attributes) {
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.attributes = attributes;
		initAttributeHashMap();
	}
	
	@Override
	public String getChromosome() {
		return chromosome;
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public void setStart(int start) {
		this.start = start;
	}

	@Override
	public int getEnd() {
		return end;
	}

	@Override
	public void setEnd(int end) {
		this.end = end;
	}

	@Override
	public char getStrand() {
		return strand;
	}
	
	@Override
	public Object getAttribute(String attributeName) {
        return attributeHashMap.get(attributeName);
	}

	@Override
	public String getAttributes() {
		return attributes;
	}

	@Override
	public String toString() {
		return chromosome + "\t" + start + "\t" + end  + "\t" + strand + "\t" + attributes;
	}
	
    private void initAttributeHashMap() {
        String[] values = attributes.split(";");
        for (String attribute : values) {
            attribute = attribute.trim();
            int equalindex = attribute.indexOf("=");
            String splitData = "=";
            if (equalindex == -1) //gtf uses space and gff3 uses =
                splitData = " ";
            String[] data = attribute.split(splitData);
            String value = "";
            if (data.length >= 2 && data[1].indexOf('"') != -1){ // an attribute field could be empty
                value = data[1].replaceAll('"' + "","").trim();
            } else if (data.length >= 2) {
                value = data[1].trim();
            }
            attributeHashMap.put(data[0].trim(), value);
        }
     }
}
