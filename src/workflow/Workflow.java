package workflow;

import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.avalon.fortress.util.dag.CyclicDependencyException;

import edu.isi.pegasus.common.util.XMLWriter;


public class Workflow {
	private String name;
	
	private LinkedList<Job> jobs;
	
	private String SCHEMA_VERSION = "0.1";
	
	protected Workflow(String name) {
		this.name = name;
		jobs = new LinkedList<Job>();
	}
	
	public Job createJob(String jobName, int cores, String cmd, String inputFile, String outputFile) {
		Job job = new Job(jobName, cores, cmd, inputFile, outputFile);
		jobs.add(job);
		return job;
	}

	public Job createJob(String jobName, int cores, String cmd, String[] inputFiles, String outputFile) {
		Job job = new Job(jobName, cores, cmd, inputFiles, outputFile);
		jobs.add(job);
		return job;
	}

	public Job getJob(String jobName) {
		for (Job job: jobs) {
			if (job.getName().equals(jobName)) 
				return job;
		}
		return null;
	}	
		
    /**
     * Generate a DAX representation on STDOUT.
     */
    public void writeToSTDOUT() {
        XMLWriter mWriter = new XMLWriter(new BufferedWriter(new OutputStreamWriter(
                System.out)));
        toXML(mWriter);
        mWriter.close();
    }

    /**
     * Generate a DAX representation and pipe it into the Writer
     *
     * @param writer A Writer object
     * @param close Whether writer should be closed on return.
     */
    public void writeToWriter(Writer writer, boolean close) {
        XMLWriter mWriter = new XMLWriter(writer);
        toXML(mWriter);
        if (close) {
            mWriter.close();
        }
    }

    /**
     * Generates a DAX representation.
     *
     * @param writer @
     */
    public void toXML(XMLWriter writer) {
        int indent = 0;
        writer.startElement("workflow");
//        writer.writeAttribute("xmlns", SCHEMA_NAMESPACE);
//        writer.writeAttribute("xmlns:xsi", SCHEMA_NAMESPACE_XSI);
//        writer.writeAttribute("xsi:schemaLocation", SCHEMA_NAMESPACE + " " + SCHEMA_LOCATION);
        writer.writeAttribute("version", SCHEMA_VERSION);
        writer.writeAttribute("name", name);
//        writer.writeAttribute("index", Integer.toString(mIndex));
//        writer.writeAttribute("count", Integer.toString(mCount));
/*
        //print notification invokes
        writer.
                writeXMLComment(
                "Section 1: Invokes - Adds notifications for a workflow (can be empty)",
                true);
        for (Invoke i : mInvokes) {
            i.toXML(writer, indent + 1);
        }
        //print file
        writer.writeXMLComment(
                "Section 2: Files - Acts as a Replica Catalog (can be empty)",
                true);
        for (File f : mFiles) {
            f.toXML(writer, indent + 1);
        }

        //print executable
        writer.
                writeXMLComment(
                "Section 3: Executables - Acts as a Transformaton Catalog (can be empty)",
                true);
        for (Executable e : mExecutables) {
            e.toXML(writer, indent + 1);
        }

        //print transformation
        writer.
                writeXMLComment(
                "Section 4: Transformations - Aggregates executables and Files (can be empty)",
                true);
        for (Transformation t : mTransformations) {
            t.toXML(writer, indent + 1);
        }
*/
        // print jobs
        writer.writeXMLComment("Section 1: Define jobs", true);
		for (Job job : jobs) {
        	job.toXML(writer, indent+1);
        }

        // print dependencies
        writer.writeXMLComment("Section 2: Dependencies - Parent Child relationships (can be empty)", true);
        for (Job job: jobs) {
        	LinkedList<Job> parentJobs = job.getParentJobs();
        	if (parentJobs != null) {
    			writer.startElement("child", indent+1);
    			writer.writeAttribute("ref", job.getName());
        		for (Job parentJob: parentJobs) {
        			writer.startElement("parent", indent+2);
        			writer.writeAttribute("ref", parentJob.getName());
        			writer.endElement();
        		}
    			writer.endElement(indent+1);
        	}
        }
        
        //end workflow
        writer.endElement();
    }
    
    /**
     * Verify a set of jobs and all their dependencies have no cycles.  All
     *  jobs in the graph must exist in the list.
     *     *
     * @throws CyclicDependencyException  if there is a cycle.
     */
    private void verifyDAG() throws CyclicDependencyException {
        // Reset the orders of all the vertices.
    	for (Job job: jobs) {
    		job.reset();
    	}
        
        // Assert that all jobs are in the job list and resolve each of their orders.
        Iterator<Job> it = jobs.iterator();
        while ( it.hasNext() )
        {
            Job v = it.next();
            
            // Make sure that any dependencies are also in the jobs list.  This adds
            //  a little bit to the load, but we don't test this and the test would have
            //  failed, this would lead to some very hard to track down problems elsewhere.
            Iterator<Job> dit = v.getChildJobs().iterator();
            while (dit.hasNext())
            {
                Job dv = dit.next();
                if ( !jobs.contains( dv ) )
                {
                    throw new IllegalStateException( "A dependent job (" + dv.getName() + ") of "
                        + "job (" + v.getName() + ") was not included in the jobs list." );
                }
            }
            
            v.resolveOrder();
        }
    }

    /**
     * Sort a set of jobs so that no dependency is before its jobs.  If
     * we have a job named "Parent" and one named "Child" that is listed as
     * a dependency of "Parent", we want to ensure that "Child" always comes
     * after "Parent".  As long as there are no cycles in the list, we can sort
     * any number of job that may or may not be related.  Both "Parent"
     * and "Child" must exist in the job list, but "Child" will also be
     * referenced as a dependency of "Parent".
     *
     * <p>
     *   <b>Implementation Detail:</b> This particular algorithm is a more
     *   efficient variation of the typical Topological Sort algorithm.  It uses
     *   a Queue (Linked List) to ensure that each edge (connection between
     *   two vertices) or vertex is checked only once.  The efficiency is
     *   O = (|V| + |E|).
     * </p>
     *
     * @param vertices
     * @throws CyclicDependencyException
     */
	public Job[] getOrderedJobExecutionArray() {
        // Verify the graph and set the job orders in the process.
        try {
			verifyDAG();
		} catch (CyclicDependencyException e) {
			e.printStackTrace();
		}
        
        // We now that there are no cycles and that each of the vertices has an order
        //  that will allow them to be sorted.
        Collections.sort(jobs);
        Collections.reverse(jobs);
        
        return jobs.toArray(new Job[0]);
	}

}
