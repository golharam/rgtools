package workflow;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.avalon.fortress.util.dag.CyclicDependencyException;

import edu.isi.pegasus.common.util.XMLWriter;

public class Job implements Comparable<Job> {
	private String jobName;
	private int cores;
	private String cmd;
	String[] inputFiles;
	String outputFile;
	private LinkedList<Job> parentJobs;
	private LinkedList<Job> childJobs;

	private String jobID = null;
	
	private int m_order;
    private boolean m_seen;
	
    protected Job(String jobName, int cores, String cmd) {
		this.jobName = jobName;
		this.cores = cores;
		this.cmd = cmd;
		childJobs = new LinkedList<Job>();
		parentJobs = new LinkedList<Job>();
    }
    
	protected Job(String jobName, int cores, String cmd, String inputFile, String outputFile) {
		this(jobName, cores, cmd);

		this.inputFiles = new String[1];
		this.inputFiles[0] = inputFile;
		this.outputFile = outputFile;
	}

	protected Job(String jobName, int cores, String cmd, String[] inputFiles, String outputFile) {
		this(jobName, cores, cmd);

		this.inputFiles = inputFiles;
		this.outputFile = outputFile;
	}

	public String getName() {
		return jobName;
	}
	
	public LinkedList<Job> getParentJobs() {
		return parentJobs;
	}
	
	public LinkedList<Job> getChildJobs() {
		return childJobs;
	}
	
	/* This is used internally to maintain tree structure */
	private void addChild(Job child) {
		childJobs.add(child);
	}
	
	public String getJobID() {
		return jobID;
	}
	
	public void setJobID(String id) {
		jobID = id;
	}

	public int getCores() {
		return cores;
	}
	
	/* A job is runnable if:
	 * 1.  The output file(s) do not exist, or
	 * 2.  the input files are newer than the output files
	 */
	public boolean isRunnable() {
		File outputFile = new File(this.outputFile);
		if (!outputFile.exists())
			return true;

		for (int i = 0; i < inputFiles.length; i++) {
			File inputFile = new File(inputFiles[i]);
			if (outputFile.lastModified() > inputFile.lastModified())
				return false;
		}
		
		return true;
	}
	
	/* Remove the parent of this */
	/*
	public void removeParent(Job parent) {
	
		for (Job p: parentJobs) {
			if (p.getName().equals(parent.getName())) {
				parentJobs.remove(p);
				
				p.childJobs.remove(this);
			}
		}
	}
	*/
	/* Use this method to add a job dependency */
	public void addParent(Job parent) {
		parentJobs.add(parent);
		parent.addChild(this);
	}
	
	public void toXML(XMLWriter writer, int indent) {
        writer.startElement("job", indent);
        writer.writeAttribute("name", jobName);
        writer.writeAttribute("cores", new Integer(cores).toString());
        
    	for (int i = 0; i < inputFiles.length; i++) {
    		writer.startElement("inputfile", indent+1);
    		writer.writeData(inputFiles[i]);
    		writer.endElement();
        }
    	writer.startElement("outputfile", indent+1);
    	writer.writeData(outputFile);
    	writer.endElement();
    	
    	writer.startElement("command", indent+1);
        writer.writeData(cmd);
        writer.endElement();
        
        writer.endElement(indent);
	}
	
	public String toString() {
		return String.format("Job %s", jobName);
	}

	public ArrayList<String> getCommandStrings() {
		ArrayList<String> cmds = new ArrayList<String>();
		cmds.add(cmd);
		return cmds;
	}

	/*
	 * The following functions are used for DAG work:
	 * reset()
	 * resolveOrder(...)
	 * 
	 */
	protected void reset() {
        m_order = 0;
        m_seen = false;
	}

	public int getOrder() {
		return m_order;
	}
    /**
     * Recurse through the tree from this vertex assigning an order to each
     *  and at the same time checking for any cyclic dependencies.
     *
     * @throws CyclicDependencyException If a cyclic dependency is discovered.
     */
    public void resolveOrder() throws CyclicDependencyException {
        resolveOrder(getName());
    }

    /**
     * Recursively searches for cycles by traveling down the dependency lists
     *  of this vertex, looking for the start vertex.
     *
     * @param path The path to the Vertex.  It is worth the load as it makes a
     *             descriptive error message possible.
     *
     * @return The highest order of any of the dependent vertices.
     *
     * @throws CyclicDependencyException If a cyclic dependency is discovered.
     */
    private int resolveOrder( String path ) throws CyclicDependencyException {
        m_seen = true;
        try
        {
            int highOrder = -1;
            for ( Iterator<Job> iter = childJobs.iterator(); iter.hasNext(); )
            {
                Job dv = iter.next();
                if ( dv.m_seen )
                {
                    throw new CyclicDependencyException( path + " -> " + dv.getName() );
                }
                else
                {
                    highOrder = Math.max(
                        highOrder, dv.resolveOrder( path + " -> " + dv.getName() ) );
                }
            }
            
            // Set this order so it is one higher than the highest dependency.
            m_order = highOrder + 1;
            return m_order;
        }
        finally
        {
            m_seen = false;
        }
    }

    @Override
	public int compareTo(Job other) {
        int orderInd;

        if ( m_order < other.m_order )
        {
            orderInd = -1;
        }
        else if ( m_order > other.m_order )
        {
            orderInd = 1;
        }
        else
        {
            orderInd = 0;
        }

        return orderInd;
	}

}
