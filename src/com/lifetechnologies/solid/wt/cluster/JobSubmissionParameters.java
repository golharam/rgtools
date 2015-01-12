package com.lifetechnologies.solid.wt.cluster;

public class JobSubmissionParameters {
	
	private String environment;
	private String queueName;
	private long memoryRequirement;
	private boolean rerunnable;
	private String resourceString;
	private String additionalOptions;
	
	public static final String DEFAULT_PBS_RESOURCE_STRING = "nodes=1:ppn=1,mem=${pbs.mem},vmem=${pbs.vmem}";
	public static final String DEFAULT_LSF_RESOURCE_STRING = "select[mem>${lsf.select.mem} && swp>${lsf.select.swp}] rusage[mem=${lsf.rusage.mem}:swp=${lsf.rusage.swp}]";
	public static final String DEFAULT_SGE_RESOURCE_STRING = "mem_free=${sge.mem},virtual_free=${sge.vmem}";
	
	public String getEnvironment() {
		return environment;
	}
	public void setEnvironment(String environment) {
		this.environment = environment;
	}
	public String getQueueName() {
		return queueName;
	}
	public void setQueueName(String queueName) {
		this.queueName = queueName;
	}
	public long getMemoryRequirement() {
		return memoryRequirement;
	}
	public void setMemoryRequirement(long memoryRequirement) {
		this.memoryRequirement = memoryRequirement;
	}
	public boolean isRerunnable() {
		return rerunnable;
	}
	public void setRerunnable(boolean rerunnable) {
		this.rerunnable = rerunnable;
	}
	public String getResourceString() {
		return resourceString;
	}
	public void setResourceString(String resourceString) {
		this.resourceString = resourceString;
	}
	public String getAdditionalOptions() {
		return additionalOptions;
	}
	public void setAdditionalOptions(String additionalArguments) {
		this.additionalOptions = additionalArguments;
	}
	public JobSubmissionParameters clone() {
		JobSubmissionParameters j = new JobSubmissionParameters();
		j.setEnvironment(this.environment);
		j.setMemoryRequirement(this.memoryRequirement);
		j.setQueueName(this.queueName);
		j.setRerunnable(this.rerunnable);
		j.setResourceString(this.resourceString);
		j.setAdditionalOptions(this.additionalOptions);
		return j;
	}

}
