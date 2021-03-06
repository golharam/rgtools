package tools.vcfbrowser;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.table.AbstractTableModel;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Interval;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;
import ca.mcgill.mcb.pcingola.vcf.VcfInfo;
import htsjdk.samtools.util.Log;


public class VarTableModel extends AbstractTableModel {
    private final Log log = Log.getInstance(VarTableModel.class);

	protected static String[] defaultColumnNames = {"CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER"};	
	protected static String[] defaultColumnDescs = {"Chromosome",
        "Position",
        "dbSNP ID",
        "Reference Base",
        "Alternate Base(s)",
        "Quality of Call",
        "Filter Status"};	
	
	private ArrayList<String> columnNames = new ArrayList<String>();
	private ArrayList<String> columnDescs = new ArrayList<String>();
	private HashMap<String, Integer> keyToDataColumnIndex = new HashMap<String, Integer>();
    protected ArrayList<Object[]> data;
    protected ArrayList<VcfEntry> variants;

	private HashMap<String, VcfInfo> vcfInfoById = new HashMap<String, VcfInfo>();
	private HashMap<String, VcfInfo> vcfFormatById = new HashMap<String, VcfInfo>();
	private ArrayList<String> keyList = new ArrayList<String>();
	private ArrayList<String> perSampleKeyList = new ArrayList<String>();
	private ArrayList<String> sampleNames = new ArrayList<String>();
	private HashMap<String, String> bamFiles = new HashMap<String,String>();	// SAMPLE, URI
 
	// TODO: Redo this to make use to VCFFile
	protected void readHeader(VcfFileIterator vcfFile) {
        // Read the header and construct a table of all the possible info fields
        VcfHeader header = vcfFile.readHeader();
    	
    	String headerLines[] = header.toString().split("\n");
		for (String line : headerLines) {
			if (line.startsWith("##INFO=")) {
				VcfInfo vcfInfo = new VcfInfo(line);
				vcfInfoById.put(vcfInfo.getId(), vcfInfo);
				keyList.add(vcfInfo.getId());
			} else
			if (line.startsWith("##FORMAT=")) {
				VcfInfo vcfInfo = new VcfInfo(line);
				vcfFormatById.put(vcfInfo.getId(), vcfInfo);
				perSampleKeyList.add(vcfInfo.getId());
			} else
			if (line.startsWith("#CHROM")) {
				String fields[] = line.split("\t");
				for (int i = 9; i < fields.length; i++)
					sampleNames.add(fields[i]);
			} else
			if (line.startsWith("##bamfile")) {
				/* Format of this line in VCF should be:
				##bamfile=<SAMPLE=1793568725,URI=http://localhost:10000/1793568725/SID1793568725.dedup.indelrealigner.recal.bam>
				*/
				// this code taken from ca.mcgill.mcb.pcingola.vcf.VcfInfo
				int start = line.indexOf('<');
				int end = line.lastIndexOf('>');
				String params = line.substring(start + 1, end);
				// Find SAMPLE
				Pattern pattern = Pattern.compile("SAMPLE=([^,]+),");
				Matcher matcher = pattern.matcher(params);
				String sample = null;
				if (matcher.find()) sample = matcher.group(1);
				else throw new RuntimeException("Cannot find 'SAMPLE' in bamfile line: '" + line + "'");

				// Find SAMPLE
				pattern = Pattern.compile("URI=(.+)");
				matcher = pattern.matcher(params);
				String URI = null;
				if (matcher.find()) URI = matcher.group(1);
				else throw new RuntimeException("Cannot find 'URI' in bamfile line: '" + line + "'");
				
				log.info("Registering " + sample + " -> " + URI);
				bamFiles.put(sample, URI);
			}
		}
	}
	
	/* Input is VCF file and number of rows to load initially.  
	 * nRows is for when the VCF file is huge and we don't need to load everything, you can specify
	 * how many rows to load.  -1 indicates to load all rows.
	 */
	public VarTableModel(File VCF, int nRows) {
		VcfFileIterator vcfFile = new VcfFileIterator(VCF.getAbsolutePath());	
		readHeader(vcfFile);
		
		// Construct an array of columns
		int i = 0;
		for (String colName : defaultColumnNames) {
			columnNames.add(colName);
			columnDescs.add(defaultColumnDescs[i]);
			keyToDataColumnIndex.put(colName, i);
			i++;
		}		
		for (String key : keyList) {
			columnNames.add(key);
			VcfInfo info = vcfInfoById.get(key);
			columnDescs.add(info.getDescription());
			keyToDataColumnIndex.put(key, i);
			i++;
		}
		for (String sample : sampleNames) {
			for (String perSampleKey : perSampleKeyList) {
    			columnNames.add(sample + "-" + perSampleKey);
    			VcfInfo info = vcfFormatById.get(perSampleKey);
    			columnDescs.add(info.getDescription());
    			keyToDataColumnIndex.put(sample + "-" + perSampleKey, i);
    			i++;
			}
		}
		log.info("Found " + columnNames.size() + " columns: ");
		for (String colName : columnNames) {
			log.info("Column: " + colName + " (" + keyToDataColumnIndex.get(colName) + ")");
		}

		// read in data in the same order as columns
		log.info("Reading data...");
		long rows = 0; 
		long rowLimit = (nRows == -1) ? rows+1 : nRows;
		data = new ArrayList<Object[]>();
		variants = new ArrayList<VcfEntry>();
		
		while (rows < rowLimit) {
			VcfEntry ve = vcfFile.next();
			
			if (ve == null) 
				break;
			
    		Object rowData[] = new Object[columnNames.size()];
    		
    		rowData[0] = ve.getChromosomeNameOri();
    		rowData[1] = ve.getStart()+1;
    		rowData[2] = ve.getId();
    		rowData[3] = ve.getRef();
    		rowData[4] = ve.getAltsStr();
    		rowData[5] = ve.getQuality();
    		rowData[6] = ve.getFilterPass();
    		
    		int colIndex = 7;
    		for (String key : keyList) {
    			rowData[colIndex] = ve.getInfo(key);
    			colIndex++;
    		}
    		
    		List<VcfGenotype> genotypes = ve.getVcfGenotypes();
    		for (VcfGenotype genotype : genotypes) {
    			for (String perSampleKey : perSampleKeyList) {
    				if (perSampleKey.equals("GT")) 
    					rowData[colIndex] = genotype.getGenotypeStr();
    				else
    					rowData[colIndex] = genotype.get(perSampleKey);	

        			colIndex++;
    			}
    		}
    		data.add(rowData);
    		variants.add(ve);
    		
    		rows++;
    		if ((rows % 10000) == 0) {
    	        log.info("Finished reading " + rows + " rows");
    		}
    		if (nRows == -1)
    			rowLimit = rows+1;
		}
        log.info("Read " + data.size() + " rows");
	}

	@Override
	public int getColumnCount() {
		return columnNames.size();
	}

	public String getColumnDescription(int colIndex) {
		return columnDescs.get(colIndex);
	}
	
	@Override
	public int getRowCount() {
		return data.size();
	}

	@Override
	public Object getValueAt(int arg0, int arg1) {
		return data.get(arg0)[arg1];
	}

	public String getColumnName(int col) {
		return columnNames.get(col);
	}

	public int getColumnIndex(String columnName) {
		return columnNames.indexOf(columnName);
	}
	
	public String[] getSampleNames() {
		return sampleNames.toArray(new String[0]);
	}

	public String getBAMFileURI(String sample) {
		return bamFiles.get(sample);
	}
	
	// Return a list of variants from data[] that are in a specific range (inclusive of the ends)
	public ArrayList<VcfEntry> getVariantsInRange(String chr, int start, int end) {
		ArrayList<VcfEntry> variantsInRange = new ArrayList<VcfEntry>();
		Interval targetInterval = new Interval(null, start, end, false, null);
		
		for (VcfEntry v : variants) {
			if (v.getChromosomeNameOri().equals(chr) && v.intersects(targetInterval)) {
				variantsInRange.add(v);
			}
		}
		return variantsInRange;
	}
}
