package workflow;

public class WorkflowFactory {
	public static Workflow createNewWorkflow(String name) {
		return new Workflow(name);
	}

	public static Workflow createNewWorkflow(String name, String version) {
		if (version.equals("0.1"))
			return new Workflow(name);
		else
			return null;
	}
}
