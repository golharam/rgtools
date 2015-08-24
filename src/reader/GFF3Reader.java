package reader;

import genome.features.GFFFeature;
import genome.features.GenomicFeatureI;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * Based on BioPerl's GFFReader
 */
public class GFF3Reader {

	private String gffFile = null;
	private BufferedReader br = null;
	
    private static final Logger log = Logger.getLogger(GFF3Reader.class.getName());
    
    public GFF3Reader(String filename) {
    	gffFile = filename;
        try {
			br = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			br = null;
		}
    }
    
    public GFF3Reader(File gffFile) {
    	this.gffFile = gffFile.getAbsolutePath();
    	
        try {
			br = new BufferedReader(new FileReader(gffFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			br = null;
		}
    }

    public GFFFeature next() {
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
	        	GFFFeature f = parseLine(s);
                return f;
	        }    	
    	} catch (IOException e) {
			e.printStackTrace();
    	}
		return null;
    }
    
    public void close() {
    	try {
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	br = null;
    }
    
    /**
     * Read a file into a GenomicFeatureList. Each line of the file becomes one Feature object.
     *
     * @param filename The path to the GFF file.
     * @return A FeatureList.
     * @throws IOException Something went wrong -- check exception detail message.
     */
    public static List<GFFFeature> read(String filename) throws IOException {
        log.info("GFF3Reader.read(): Reading " + filename);
        BufferedReader br = new BufferedReader(new FileReader(filename));
        return read(br);
    }

    public static List<GFFFeature> read(File file) throws IOException {
        log.info("GFF3Reader.read(): Reading " + file.getAbsolutePath());
        BufferedReader br = new BufferedReader(new FileReader(file));
        return read(br);
    }
    
    private static List<GFFFeature> read(BufferedReader br) throws IOException {
    	List<GFFFeature> features = new LinkedList<GFFFeature>();

        String s;
        for (s = br.readLine(); null != s; s = br.readLine()) {
            s = s.trim();

            if (s.length() > 0) {
                if (s.charAt(0) == '#') {
                    //ignore comment lines
                    if(s.startsWith("##fasta"))
                        break;
                } else {
                	GFFFeature f = parseLine(s);
                    if (f != null) {
                        features.add(f);
                    }
                }
            }

        }

        br.close();
        return features;
    }
    
    /**
     * Read a file into a FeatureList. Each line of the file becomes one Feature object.
     *
     * @param filename The path to the GFF file.
     * @return A FeatureList.
     * @throws IOException Something went wrong -- check exception detail message.
     */
    public static List<GFFFeature> readByFeatureType(String filename, String featureType) throws IOException {
        log.info("GFF3Reader.readByFeatureType(): Reading " + filename + " for " + featureType);

        List<GFFFeature> features = new LinkedList<GFFFeature>();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        int line = 0;
        String s;
        for (s = br.readLine(); null != s; s = br.readLine()) {
            s = s.trim();
            line++;
            if (s.length() > 0) {
                if (s.charAt(0) == '#') {
                    //ignore comment lines
                    if(s.startsWith("##fasta"))
                        break;
                } else {

                	GFFFeature f = parseLine(s);
                    if (f == null) 
                    	log.warning("Incorrectly formatted feature on line " + line + ": " + s );
                    if ((f != null) && (f.getType().equals(featureType))) {
                        features.add(f);
                    }
                }
            }

        }

        br.close();
        return features;
    }
    
    /**
     * create Feature from line of GFF file
     */
    private static GFFFeature parseLine(String s) {
    	String[] fields = s.split("\t");

    	// GFF files have 9 tokens
        if (fields.length < 9) {
        	log.warning("Line contains less than 9 fields:\n" + s);
            return null;
        }

        String seqname = fields[0];
        String source = fields[1];
        String type = fields[2];
        String locStart = fields[3];
        String locEnd = fields[4];
        Double score;
        try {
        	score = Double.parseDouble(fields[5]);
        } catch (Exception e) {
        	score = 0.0;
        }
        char strand = fields[6].charAt(0);
        
        //added by scooter willis to deal with glimmer predictions that
        //have the start after the end but is a negative strand
        int locationStart = Integer.parseInt(locStart);
        int locationEnd = Integer.parseInt(locEnd);
        if (locationStart > locationEnd) {
            int temp = locationStart;
            locationStart = locationEnd;
            locationEnd = temp;
        }
        
        char frame = fields[7].charAt(0);
        String attributes = new String(fields[8]);
        
        return new GFFFeature(seqname, locationStart, locationEnd, strand, type, source, score, frame, attributes);
    }

    public static void main(String args[]) throws Exception {
    	// Two ways to read a GFF file.  
    	// 1.  Load all at once:
        List<GFFFeature> featureList = GFF3Reader.read("test_data/reference/hg19/refGene.hg19.gtf");
        System.out.println("Features");
        for (GenomicFeatureI feature : featureList) {
            System.out.println(feature);
        }

        // 2.  Open file and load one by one
        GFF3Reader gffReader = new GFF3Reader("test_data/reference/hg19/refGene.hg19.gtf");
        System.out.println("Features");
        GFFFeature feature = null;
        while ((feature = gffReader.next()) != null) {
            System.out.println(feature);
        }
        gffReader.close();
    }
}
