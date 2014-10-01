package reader;


import genome.features.GenomicFeature;
import java.io.*;
import java.util.*;
import java.util.logging.Logger;

import org.biojava3.genome.parsers.gff.Feature;
import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.Location;

/**
 * http://www.bioperl.org/wiki/GTF
 * Read and write FeatureLists as GFF/GTF formatted files.
 *<br><br>
 * The GFF moniker is applied to a variety of tab-delimited formats
 * that mock the notion of a standard. This class should parse most files
 * bearing at least a passing resemblance to any of the formats. You will, however, need
 * to research the semantics of the files you encounter. Generally,
 * the format consists of 9 tab-delimited fields:
 * <br>
 * <pre>
 * seqname   source   featureType   start   end   score   strand   frame   attributes
 * </pre>
 * The 9th field consists of key-value pairs separated by semicolons, the first of which JavaGene interprets
 * as the group id (as used in GFF1). It is the precise meaning of this 9th field that
 * varies from week to week. The Feature and FeatureList objects provide various utility methods to
 * ease the task of accessing and using the attributes. The proper interpretation of any
 * particular attribute, however, is left to you.
 *
 * @author Hanno Hinsch
 */
public class BEDReader {

	private BufferedReader br = null;
	
    private static final Logger log = Logger.getLogger(BEDReader.class.getName());

    public BEDReader(String filename) {
        try {
			br = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			br = null;
		}
    }
    
    public GenomicFeature next() {
    	if (br == null)
    		return null;
    	
    	try {
	        String s = br.readLine();
	        if (s == null)
	        	return null;
	        
	        while (s.charAt(0) == '#')
	        	s = br.readLine();
	        
	        s = s.trim();
	
	        if (s.length() > 0) {
                GenomicFeature f = parseLine(s);
                return f;
	        }    	
    	} catch (IOException e) {
			e.printStackTrace();
    	}
		return null;
    }
    
    public static FeatureList read(String filename) throws IOException {
    	log.info("Reading " + filename);
        BufferedReader br = new BufferedReader(new FileReader(filename));
        return readAsFeatureList(br);
    }
    
    /**
     * Read a file into a FeatureList. Each line of the file becomes one Feature object.
     *
     * @param filename The path to the GFF file.
     * @return A FeatureList.
     * @throws IOException Something went wrong -- check exception detail message.
     */
    public static List<GenomicFeature> readOld(String filename) throws IOException {
        log.info("Reading " + filename);
        BufferedReader br = new BufferedReader(new FileReader(filename));
        return read(br);
    }

    public static List<GenomicFeature> read(File file) throws IOException {
        log.info("Reading " + file.getAbsolutePath());
        BufferedReader br = new BufferedReader(new FileReader(file));
        return read(br);
    }
    
    private static List<GenomicFeature> read(BufferedReader br) throws IOException {
    	List<GenomicFeature> features = new LinkedList<GenomicFeature>();

        String s;
        for (s = br.readLine(); null != s; s = br.readLine()) {
            s = s.trim();

            if (s.length() > 0) {
                if (s.charAt(0) == '#') {
                    //ignore comment lines
                    if(s.startsWith("##fasta"))
                        break;
                } else {

                    GenomicFeature f = parseLine(s);
                    if (f != null) {
                        features.add(f);
                    }
                }
            }

        }

        br.close();
        return features;
    }
    
    private static FeatureList readAsFeatureList(BufferedReader br) throws IOException {
        FeatureList features = new FeatureList();

        String s;
        for (s = br.readLine(); null != s; s = br.readLine()) {
            s = s.trim();

            if (s.length() > 0) {
                if (s.charAt(0) == '#') {
                    //ignore comment lines
                    if(s.startsWith("##fasta"))
                        break;
                } else {
                    FeatureI f = parseFeatureLine(s);
                    if (f != null) {
                        features.add(f);
                    }
                }
            }

        }

        br.close();
        return features;
    }

    private static Feature parseFeatureLine(String s) {
        //FIXME better errors on parse failures

    	String[] fields = s.split("\t");

    	// BED files have at least 3 tokens
        if (fields.length < 3) {
        	log.warning("Line contains less than 3 fields:\n" + s);
            return null;
        }

        String seqname = fields[0];
        String locStart = fields[1];
        String locEnd = fields[2];
        
        int locationStart = Integer.parseInt(locStart);
        int locationEnd = Integer.parseInt(locEnd);
        if(locationStart > locationEnd){
            int temp = locationStart;
            locationStart = locationEnd;
            locationEnd = temp;
        }
        //added by scooter willis to deal with glimmer predictions that
        //have the start after the end but is a negative strand
        char strand = '+';
        Location location = Location.fromBio(locationStart, locationEnd, strand);
        
        if (fields.length >= 4) {
        	return new Feature(seqname, "", "", location, 0.0, -1, "name="+fields[3]);
        } else {
        	return new Feature(seqname, "", "", location, 0.0, -1, "");
        }
    }

    /**
     * create Feature from line of GFF file
     */
    private static GenomicFeature parseLine(String s) {
        //FIXME better errors on parse failures

    	String[] fields = s.split("\t");

    	// BED files have at least 3 tokens
        if (fields.length < 3) {
        	log.warning("Line contains less than 3 fields:\n" + s);
            return null;
        }

        String seqname = fields[0];
        String locStart = fields[1];
        String locEnd = fields[2];
        
        int locationStart = Integer.parseInt(locStart);
        int locationEnd = Integer.parseInt(locEnd);
        if(locationStart > locationEnd){
            int temp = locationStart;
            locationStart = locationEnd;
            locationEnd = temp;
        }
        
        String attributes = "";
        char strand = '+';
        if (fields.length == 6) {
        	attributes = "name="+fields[3]+";";
        	strand = fields[5].charAt(0);
        }
        
        return new GenomicFeature(seqname, locationStart-1, locationEnd, strand, attributes);
    }
}
