package workflow;

import java.io.IOException;
import java.util.LinkedList;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


public class Parser extends DefaultHandler {
	private Workflow workflow = null;
	
	boolean inJob = false;
	boolean inJobInputFile = false;
	boolean inJobOutputFile = false;
	boolean inJobCommand = false;
	String jobName = null;
	int jobCores = 0;
	LinkedList<String> jobInputFiles = null;
	LinkedList<String> jobOutputFiles = null;
	String jobCommand = null;

	boolean inJobDependency = false;
	String childJobName = null;
	
	private void addDependency(String parentJobName, String childJobName) throws Exception {
		Job childJob = workflow.getJob(childJobName);
		Job parentJob = workflow.getJob(parentJobName);
		childJob.addParent(parentJob);		
	}

	//Event Handlers		
	public void	startElement(String uri, String localName, String qName, Attributes attributes) {
		if (qName.equals("workflow")) {
			workflow = WorkflowFactory.createNewWorkflow(attributes.getValue("name"), attributes.getValue("version"));
		} else if (qName.equals("job")) {
			inJob = true;
			jobName = attributes.getValue("name");
			jobCores = Integer.parseInt(attributes.getValue("cores"));
		} else if (qName.equals("inputfile")) {
			inJobInputFile = true;
		} else if (qName.equals("outputfile")) {
			inJobOutputFile = true;
		} else if (qName.equals("command")) {
			inJobCommand = true;
			jobCommand = new String();
		} else if (qName.equals("child")) {
			inJobDependency = true;
			childJobName = attributes.getValue("ref");
		} else if (qName.equals("parent")) {
			String parentJobName = attributes.getValue("ref");
			try {
				addDependency(parentJobName, childJobName);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} /* else {
			// we should never get here.
			System.out.println("start element: " + qName);
			if (attributes.getLength() > 0) {
				for (int i = 0; i < attributes.getLength(); i++) {
					System.out.println("\tattribute: " + attributes.getQName(i) + " = " + attributes.getValue(i));
				}
			}			
			
		}
		*/
	}

	 public void characters(char[] ch, int start, int length) {
		 if (inJobInputFile) {	 
			 String data = new String(ch, start, length);
			 
			 if (jobInputFiles == null)
				 jobInputFiles = new LinkedList<String>();
			 
			 jobInputFiles.add(data);
		 } else if (inJobOutputFile) {
			 String data = new String(ch, start, length);

			 if (jobOutputFiles == null)
				 jobOutputFiles = new LinkedList<String>();
			 
			 jobOutputFiles.add(data);
		 } else if (inJobCommand) {
			 jobCommand = jobCommand.concat(new String(ch, start, length));
		 } 
	 }
	 
	 public void endElement(String uri, String localName, String qName) {
		 if ((inJob) && (qName.equals("job"))) {
			inJob = false;
			workflow.createJob(jobName, jobCores, jobCommand, jobInputFiles.toArray(new String[0]), jobOutputFiles.get(0));

			jobName = null;
			jobCores = 0;
			jobInputFiles = null;
			jobOutputFiles = null;
		 } else if ((inJobInputFile) && (qName.equals("inputfile"))) {
			 inJobInputFile = false;
		 } else if ((inJobOutputFile) && (qName.equals("outputfile"))) {
			 inJobOutputFile = false;
		 } else if ((inJobCommand) && (qName.equals("command"))) {
			 inJobCommand = false;
		 } else if ((inJobDependency) && (qName.equals("child"))) {
			 inJobDependency = false;
		 }
	 }
	 
	 /* Only this function should be called
	 * SAXParser calls the event handlers in this class to recreate the XML workflow
	 */
	public Workflow parse(String workflowFile) {
		//get a factory
		SAXParserFactory spf = SAXParserFactory.newInstance();
		try {
			//get a new instance of parser
			SAXParser sp = spf.newSAXParser();

			//parse the file and also register this class for call backs
			sp.parse(workflowFile, this);

		}catch(SAXException se) {
			se.printStackTrace();
		}catch(ParserConfigurationException pce) {
			pce.printStackTrace();
		}catch (IOException ie) {
			ie.printStackTrace();
		}		
		return workflow;
	}

}
