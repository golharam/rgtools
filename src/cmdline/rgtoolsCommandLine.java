package cmdline;

import htsjdk.samtools.util.Log;

/**
 * 
 * This is the main class of rgtools and is the way of executing individual command line programs.
 *
 * CommandLinePrograms are listed in a single command line interface based on the java package specified to instanceMain.
 *
 * If you want your own single command line program, extend this class and give instanceMain a new list of java packages in which to
 * search for classes that extend CommandLineProgram.
 *
 */
public class rgtoolsCommandLine {
    private static final Log log = Log.getInstance(rgtoolsCommandLine.class);

	public static void main(String[] args) {

	}

}
