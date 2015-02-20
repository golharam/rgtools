package tools;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Rectangle;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableRowSorter;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JSeparator;
import javax.swing.RowFilter;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JTable;
import javax.swing.JScrollPane;

import tools.vcfbrowser.AppSettings;
//import chromviewer.JDAP.JDAPSample;


import tools.vcfbrowser.AppSettings;
import tools.vcfbrowser.VarTableModel;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

// Version 0.01: initial version
// Version 0.01a: ???
// Version 0.01b: Allow specifying VCF File to load on command line
// Version 0.01c: Load IGV only upon double-clicking a row
// Version 0.01d: Added ability to click a column header to sort the table
// Version 0.01e: Added ability to filter on gene column only
//				  Added Go To menu item to go to specific row
// Version 0.01f: Show each sample's AD/DP/GT/GP fields instead of just the genotype basecall
// 				  Allow mouseover popups on column header to give description of what the column is
// Version ?????: Added Menu Item to Start IGV
//				  Added ability to allow user to specify BAM file location(s) through menu item
//				  Added command line option to not prompt to load BAM files
// 				  (change the command line option to the options parser)
//				  Changed behavior of StartIGV Menu Item to connect to running IGV.
//				  Added ability to specify column ordering from command-line when loading VCF file using --priortizeColumns
//				  Added ability to hide column from command-line when loading VCF file using --hideColumns
//				  Added status bar (panel)
// 				  Implemented Recent menu items
//				  Added ability to load only a specified # of rows

// TODO: Add ability to export VCF as text file using VCFToTab
// TODO: Implement row labels with row number
// TODO: Implement status bar indicating selected row number
// TODO: Implement advanced filtering (JEXL Support)
// TODO: Have VCFBrowser automatically update the VCF file with the information  about BAM location
// TODO: Add ability to annotate variants
// TODO: Filter variants that match the reference for all samples as the variants may have
//       been called with other samples NOT in the VCF file and the variants may not apply
//       to the samples in the VCF file.
// TODO: Allow split windows for browsing
// TODO: Add ability to rename columns

public class VCFBrowser extends JFrame {
	private static String APP_NAME = "eVCFerator - The VCF Browser";
	private static String VERSION = "0.01f";
	
	private static final Log log = Log.getInstance(VCFBrowser.class);

	private JPanel contentPane;

	private File workingDirectory = null;

	AppSettings m_appSettings = null;
	
	JMenu mnRecent = null;
	int mruListSize = 20; // allow 20 items in the Most Recently Used list
	//JMenuItem[] mruMenuItem = null;
	
	private JTable table;
	private VarTableModel varTableModel;
	
	private JLabel statusLabel;
	
	String prioritizeColumns = null;
	String hideColumns = null;
	boolean loadBAMs = true;
	int loadRows = -1;
	String[] bamFiles;
	TableRowSorter<VarTableModel> sorter;
	RowFilter<VarTableModel, Object> rowFilter = null;
	String filterText;
	
//	private Socket igv_socket = null;
//    private PrintWriter igv_out = null;
//    private BufferedReader igv_in = null;
	private boolean m_bUseIGV = false;
	
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		log.info(APP_NAME + " Version " + VERSION);
		
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					VCFBrowser frame = new VCFBrowser();
					if (args.length != 0) {
						frame.parseOptions(args);
					}
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public VCFBrowser() {
        loadSettings();

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 450, 300);
		
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);
		
		JMenuItem mntmOpenVcf = new JMenuItem("Open VCF");
		mntmOpenVcf.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onFileOpenVCFMenuItemSelected(e);
			}
		});
		mnFile.add(mntmOpenVcf);

		mnRecent = new JMenu("Open Recent");
		mnFile.add(mnRecent);
		setMRU();		

		
		JMenuItem mntmExit = new JMenuItem("Exit");
		mntmExit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				System.exit(0);
			}
		});
		mnFile.add(mntmExit);
		
		JMenu mnEdit = new JMenu("Edit");
		menuBar.add(mnEdit);
		
		JMenuItem mntmFind = new JMenuItem("Find");
		mntmFind.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onEditFind(e);
			}
		});
		mnEdit.add(mntmFind);
		
		JMenuItem mntmGoTo = new JMenuItem("Go to");
		mntmGoTo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onEditGoTo(e);
			}
		});
		mnEdit.add(mntmGoTo);

		JMenu mnIGV = new JMenu("IGV");
		menuBar.add(mnIGV);
		
		JMenuItem mntmStartIGV = new JMenuItem("Connect to IGV");
		mntmStartIGV.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onConnectToIGV(e);
			}
		});
		mnIGV.add(mntmStartIGV);
		
		JMenuItem mntmSpecifyBAMFiles = new JMenuItem("Specify BAM File(s)");
		mntmSpecifyBAMFiles.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onSpecifyBAMFiles(e);
			}
		});
		mnIGV.add(mntmSpecifyBAMFiles);
		
		JMenu mnHelp = new JMenu("Help");
		menuBar.add(mnHelp);
		
		JMenuItem mntmAbout = new JMenuItem("About");
		mnHelp.add(mntmAbout);
		mntmAbout.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onHelpAbout(e);
			}
		});
		
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		
		table = new JTable() {
			//Implement table header tool tips. 
            protected JTableHeader createDefaultTableHeader() {
                return new JTableHeader(columnModel) {
                    public String getToolTipText(MouseEvent e) {
                        String tip = null;
                        java.awt.Point p = e.getPoint();
                        int index = columnModel.getColumnIndexAtX(p.x);
                        int realIndex = columnModel.getColumn(index).getModelIndex();
                        return onGetToolTipText(realIndex);
                    }
                };
            }		
        };
		table.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				//if (e.getButton() == java.awt.event.MouseEvent.BUTTON3) {
				//	System.out.println("Right Click");
				//} else
				if (e.getClickCount() == 2) {
					int row = table.getSelectedRow();
					onRowDoubleClicked(row);
				}
			}
		});
		table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				updateStatus();
			}
		});
		table.setAutoCreateRowSorter(true);
		
		JScrollPane scrollPane = new JScrollPane(table, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		contentPane.add(scrollPane, BorderLayout.CENTER);

		table.setFillsViewportHeight(true);
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		scrollPane.setViewportView(table);
		
		JPanel statusPanel = new JPanel();
		contentPane.add(statusPanel, BorderLayout.SOUTH);
		statusPanel.setPreferredSize(new Dimension(contentPane.getWidth(), 16));
		statusPanel.setLayout(new BoxLayout(statusPanel, BoxLayout.X_AXIS));
		
		statusLabel = new JLabel("Ready.");
		statusLabel.setHorizontalAlignment(SwingConstants.LEFT);
		statusPanel.add(statusLabel);
		
		try {
			UIManager.setLookAndFeel(
			        UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (InstantiationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalAccessException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (UnsupportedLookAndFeelException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}

    private void loadSettings() {
        File userData = new File(System.getProperty("user.home") + File.separator + ".VCFBrowser.dat");
        FileInputStream in;
        ObjectInputStream s;

        m_appSettings = null;

        try {
            in = new FileInputStream(userData);
            s = new ObjectInputStream(in);
            m_appSettings = (AppSettings)s.readObject();
            s.close();
            in.close();
        } catch (FileNotFoundException e) {
            m_appSettings = null;
        } catch (IOException e) {
            m_appSettings = null;
        } catch (ClassNotFoundException e) {
            m_appSettings = null;
        }

        if (m_appSettings == null)
            m_appSettings = new AppSettings();
    }

    private void saveSettings() {
        File userData = new File(System.getProperty("user.home") + File.separator + ".VCFBrowser.dat");
        FileOutputStream f;
        ObjectOutputStream s;

        try {
            f = new FileOutputStream(userData);
            s = new ObjectOutputStream(f);

            s.writeObject(m_appSettings);
            s.flush();
            s.close();
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        }
    }
    
	// TODO: Implement Command Line Options Parser
	protected void parseOptions(final String[] args) {
		int i = 0;
		
		if (args == null)
			return;
		
		while (i < args.length) {
			log.info("Parsing arg " + i + " " + args[i]);
			if (args[i].equalsIgnoreCase("--noloadbam")) {
				this.loadBAMs = false;
			} else if (args[i].equalsIgnoreCase("-h")) {
				printUsage();
				System.exit(0);
			} else if (i == args.length-1) {
				loadVCFFile(new File(args[i]));
				i++;
			} else if (args[i].equalsIgnoreCase("--prioritizeColumns")) {
				this.prioritizeColumns = args[i+1];
				i++;
			} else if (args[i].equalsIgnoreCase("--hideColumns")) {
				this.hideColumns = args[i+1];
				i++;
			} else if (args[i].equalsIgnoreCase("--loadRows")) {
				this.loadRows = Integer.parseInt(args[i+1]);
				i++;
			} else {
				log.error("Unknown command line option: " + args[i]);
				printUsage();
				System.exit(-1);
			}
			i++;
		}
	}
	
	protected static void printUsage() {
		log.info("java -jar VCFBrowser.jar [options] [vcf file]");
		log.info("     --hideColumns 	    [comma separated list of columns] List of columns to hide in browser");
		log.info("     --loadRows			Specify # of rows to load");
		log.info("     --noloadbam    		Do not prompt for BAM file locations");
		log.info("     --prioritizeColumns  [comma separated list of columns] List of columns to display leftmost");
		log.info("     -h             	    This help message");
	}

	private void onClearRecentHistory() {
		m_appSettings.setMRUList(null);
        saveSettings();
        setMRU();
	}
	
	// Update the menu to reflect the MRU list
    private void setMRU() {
        String[] mruList = m_appSettings.getMRUList();

        // remove all the recent menu items
        mnRecent.removeAll();
        
        // rebuild the standard items
        mnRecent.add(new JSeparator());
        JMenuItem mntmClearHistory = new JMenuItem("Clear History");
        mntmClearHistory.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onClearRecentHistory();
			}
		});
        mnRecent.add(mntmClearHistory);
        
        if (mruList == null)
        	return;
        
        // Set the mru list in the menu
    	for (int i = 0; i < mruList.length; i++) {    		
    		final String text = mruList[i];
    		JMenuItem mruMenuItem = new JMenuItem(mruList[i]);
    		mruMenuItem.addActionListener(new ActionListener() {
    			public void actionPerformed(ActionEvent e) {
    				mruMenuItemActionPerformed(text);
    			}
    		});
    		mnRecent.insert(mruMenuItem, i);
    	} 
    }

    // Update the MRU list by shifting everything down one
    private void updateMRU(String file) {
        // update the list to include the new file
        String[] mruList = m_appSettings.getMRUList();
        
        // if there is no MRU, add this as the first item
        if (mruList == null) {
            mruList = new String[1];
            mruList[0] = file;
        } else {
        	// if there is an MRU, add file as the first item, then remove any other instances of it
        	// then make sure the array is not bigger than mruListSize
        	
        	LinkedList<String> a = new LinkedList<String>(Arrays.asList(mruList));
        	
        	a.add(0, file);
        	int i = 1;
        	while (i < a.size()) {        		
        		if (file.equals(a.get(i)))
        			a.remove(i);
        		else 
        			i++;
        	}
        	
        	while (a.size() > mruListSize) {
        		a.remove(a.size() - 1);
        	}
        	
            mruList = a.toArray(new String[0]);
        }

        m_appSettings.setMRUList(mruList);
        saveSettings();
        setMRU();
    }
    
	protected void mruMenuItemActionPerformed(String mruPathname) {
        File mruFile = new File(mruPathname);
        
        if (!mruFile.exists()) {
            log.error("Unable to open mruFile: " + mruPathname);
            return;
        }

        loadVCFFile(mruFile);		
	}
	
	protected void onRowDoubleClicked(int row) {
		// If not already connected to IGV,
		// determine if we should try to connect and if we should, then connect, and initialize
		// issue navigation commands
		/*
		if (igv_out == null) {
			boolean worthLoading = false;
			for (String bamFile : bamFiles) {
				if (bamFile != null)
					worthLoading = true;
			}
			if (!worthLoading) {
				log.warn("No bam files specified.  Not worth loading IGV");
				return;
			}

			try {
				if (initIGV() == false) {
					return;
				}
			} catch (Exception e) {
				return;
			}
			
		}
		*/
		if (m_bUseIGV) {
//			private Socket igv_socket = null;
//		    private PrintWriter igv_out = null;
//		    private BufferedReader igv_in = null;
			try {
				Socket igv_socket = new Socket("127.0.0.1", 60151);
				if (igv_socket.isConnected()) {
					PrintWriter igv_out = new PrintWriter(igv_socket.getOutputStream(), true);
					//BufferedReader igv_in = new BufferedReader(new InputStreamReader(igv_socket.getInputStream()));
				
					String chr = (String) table.getValueAt(row, 0);
					Integer pos = (Integer) table.getValueAt(row, 1);
					
					log.info("Navigating IGV to " + chr + ": " + pos);
				
					igv_out.println("goto " + chr + ":" + pos);
/*					
					String response = null;
					try {
						igv_in.
						response = igv_in.readLine();
						if (response.equals("OK")) {
							log.info("IGV done");
						} else {
							log.error("IGV Error: " + response);
						}
					} catch (IOException e) {
						log.error("Unable to navigate IGV");
					}
*/					
					igv_out.close();
					//igv_in.close();
					igv_socket.close();
				}				
			} catch (UnknownHostException e1) {
				log.error("Unable to connect to IGV");
				return;
			} catch (IOException e1) {
				log.error("Unable to connect to IGV");
				return;
			}
		}
	}
	
	private class CustomFileFilter extends FileFilter {
		private String fileType;
		private String fileTypeDesc;
		
		public CustomFileFilter(String extension, String description) {
			// ex ".vcf", "VCF File"
			fileType = extension;
			fileTypeDesc = description;
		}
		@Override
		public boolean accept(File arg0) {
			if (arg0.isDirectory())
				return true;
			
			if (arg0.getName().endsWith(fileType))
				return true;
			
			return false;
		}

		@Override
		public String getDescription() {
			return fileTypeDesc;
		}
	}	

	protected void onFileOpenVCFMenuItemSelected(ActionEvent e) {
		final JFileChooser fc = new JFileChooser();
		fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		fc.setFileFilter(new CustomFileFilter(".vcf", "VCF File"));
		if (workingDirectory != null)
			fc.setCurrentDirectory(workingDirectory);
		
		int returnVal = fc.showOpenDialog(this);
		
		if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            //This is where a real application would open the file.
            log.info("Opening: " + file.getName());
            workingDirectory = fc.getCurrentDirectory();
            loadVCFFile(file);
        }
	}
		
	private void onEditFind(ActionEvent evt) {
		filterText = JOptionPane.showInputDialog(this, "Enter Gene name (in regex) to search for or nothing to clear filtering.", filterText);
		
		if ((filterText == null) || (filterText.length() == 0)) {
			if (rowFilter != null) {
				log.info("Removing filtering");
				rowFilter = null;
				filterText = null;
				sorter.setRowFilter(rowFilter);
			}	
		} else {
			log.info("Filtering table on gene name " + filterText);
			try {
				rowFilter = RowFilter.regexFilter(filterText, varTableModel.getColumnIndex("Gene"));
			} catch (java.util.regex.PatternSyntaxException e) {
				log.error("Unable to parse filter regex/text.");
				return;
			}
			sorter.setRowFilter(rowFilter);
		}
	}
	
	private void onEditGoTo(ActionEvent evt) {
		String str = JOptionPane.showInputDialog(this, "Enter line number");
		
		if ((str != null) && (str.length() != 0)) {
			int lineNumber = Integer.parseInt(str);
			table.setRowSelectionInterval(lineNumber - 1, lineNumber - 1);
			Rectangle rect = table.getCellRect(lineNumber, 0, true);
			table.scrollRectToVisible(rect);
		}
	}
	
	private String onGetToolTipText(int colIndex) {
        return varTableModel.getColumnDescription(colIndex);
	}
	
	protected void onConnectToIGV(ActionEvent evt) {
		/*
		// If IGV is already running, simply connect to it.
		try {
			connectToIGV();
		} catch (Exception e) {
			log.error(e.getMessage());
		}*/
		m_bUseIGV = true;
		log.info("Will attempt to use IGV");
	}
	
	protected void onSpecifyBAMFiles(ActionEvent evt) {
        String[] samples = varTableModel.getSampleNames();
        log.info("Found " + samples.length + " samples:");
        for (String sampleName : samples) {
        	log.info("Sample: " + sampleName);
        }
        
        // Get BAM files corresponding to samples
        int i = 0;
        bamFiles = new String[samples.length];
        for (String sample : samples) {
        	//bamFiles[i] = "http://10.228.81.46/igv/igvdata/RD1000/" + sample + "/SID" + sample + ".dedup.indelrealigner.recal.bam";
        	
        	// If the BAM file is specified in the VCF file, read it from there,
        	// else ask the user for the location of the BAM file(s)
        	
        	String bamFile = varTableModel.getBAMFileURI(sample);
        	if (bamFile != null) {
        		bamFiles[i] = bamFile;
        	} else {
	            JOptionPane.showMessageDialog(this, "Select BAM file(s) corresponding to sample " + sample, "Select BAM File(s)", JOptionPane.PLAIN_MESSAGE);
	    		
	        	final JFileChooser fc = new JFileChooser();
	    		fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
	    		fc.setFileFilter(new CustomFileFilter(".bam", "BAM File"));
	    		if (workingDirectory != null)
	    			fc.setCurrentDirectory(workingDirectory);
	    		
	    		int returnVal = fc.showOpenDialog(this);
	    		if (returnVal == JFileChooser.APPROVE_OPTION) {
	                File file = fc.getSelectedFile();
	                bamFiles[i] = file.getAbsolutePath();
	    		} else {
	    			bamFile = JOptionPane.showInputDialog(this, "Enter BAM file URI for sample " + sample);
	    			bamFiles[i] = bamFile;
	    		}
        	}
            log.info(sample + " -> " + bamFiles[i]);
            i++;
        }
	}
	
	protected void onHelpAbout(ActionEvent evt) {
		JOptionPane.showMessageDialog(this, APP_NAME + " Version " + VERSION);
	}
	
	protected void loadVCFFile(File VCF) {
		// If another VCF file is loaded, unload it prior
		// to loading a new one to free up memory.
		// This is basically a reset.
		if (varTableModel != null) {
			varTableModel = null;
			sorter = null;
			rowFilter = null;
			table.setRowSorter(null);
			table.setModel(new DefaultTableModel());
			filterText = null;
			System.gc();
		}
		
		IOUtil.assertFileIsReadable(VCF);
		varTableModel = new VarTableModel(VCF, loadRows);
        table.setModel(varTableModel);

        sorter = new TableRowSorter<VarTableModel>(varTableModel);
		table.setRowSorter(sorter);
		
		if (prioritizeColumns != null) {
			log.info("Re-ordering columns based on priorty string: " + prioritizeColumns);
			String[] columns = prioritizeColumns.split(",");
			int newIndex = 0;
			
			for (String column : columns) {
				// get view index of column
				TableColumn col = null;
				try {
					col = table.getColumn(column);
				} catch (Exception e) {
					log.error("Error getting column " + column + ".  Skipping...");
					log.error(e.getMessage());
					continue;
				}
				int modelIndex = col.getModelIndex();
				int oldIndex = table.convertColumnIndexToView(modelIndex);
				log.info("Moving column " + column + " " + oldIndex + " -> " + newIndex);
				table.moveColumn(oldIndex, newIndex);
				newIndex++;
			}
		}
		
		if (hideColumns != null) {
			log.info("Hiding columns based on hide string: " + hideColumns);
			String[] columns = hideColumns.split(",");
			for (String column : columns) {
				TableColumn tcolumn = null;
				try {
					tcolumn = table.getColumn(column);
				} catch (java.lang.IllegalArgumentException e) {
					log.error("Unknown column: " + column + ", skipping...");
					continue;
				}
				tcolumn.setMinWidth(0);
				tcolumn.setMaxWidth(0);
				tcolumn.setPreferredWidth(0);
			}
		}
		
        try {
			this.setTitle(VCF.getCanonicalPath());
		} catch (IOException e) {
			log.error("Unable to get path of VCF file");
		}
        
        // Get sample names from VCF file
        if (loadBAMs) {
        	onSpecifyBAMFiles(null);
        }

        updateMRU(VCF.getAbsolutePath());
        updateStatus();
	}

	private void updateStatus() {
		int selectedRows = table.getSelectedRowCount();		
        statusLabel.setText("Rows: " + varTableModel.getRowCount() + "     " + "Selected Rows: " + selectedRows);
	}
	/*
	private boolean connectToIGV() throws Exception {

		if (isIGVRunning()) {
			log.info("Initializing IGV communication...");
			igv_socket = new Socket("127.0.0.1", 60151);
			igv_out = new PrintWriter(igv_socket.getOutputStream(), true);
			igv_in = new BufferedReader(new InputStreamReader(igv_socket.getInputStream()));
			log.info("done");	
		}

		
		return false;
	}

	private boolean initIGV() throws Exception {
		if (!isIGVRunning()) {
			if (!startIGV())
				return false;
		}

		if (igv_socket != null) 
			return true;
				
		if (connectToIGV()) {		
			log.info("Setting IGV genome to hg19...");
	        igv_out.println("genome hg19");
	        String response = igv_in.readLine();
	        if (response.equals("OK")) {
	        	log.info("done");
	        } else {
	        	log.error("Error: " + response);
	        	return false;
	        }
	        
	        for (String bamfile : bamFiles) {
	    		log.info("Loading IGV bam file " + bamfile);
	            igv_out.println("load " + bamfile);
	            response = igv_in.readLine();
	            if (response.equals("OK")) {
	            	log.info("done");
	            } else {
	            	log.error("Error: " + response);
	            	return false;
	            }
	        }
	        return true;
		}
		return false;
	}

	private boolean isIGVRunning() {
		log.info("Checking if IGV is already running...");
		try {
			igv_socket = new Socket("127.0.0.1", 60151);
			igv_socket.close();
			igv_socket = null;
		} catch (IOException e) {
			log.info("IGV is not running.");
			return false;
		}
		log.info("IGV is running.");
		return true;
	}
	
	private boolean startIGV() {
		log.info("Starting IGV...");

		if (isIGVRunning()) {
			log.info("IGV is already running");
			return true;
		}
		
		try {
			// TODO: Read IGV location from properties file
			Runtime.getRuntime().exec("/Users/golharr/bin/igv");

			// Wait for IGV to load before trying to communicate with it.
			int connectTry = 0;
			while (connectTry <= 5) {
				Thread.sleep(5000);
				if (isIGVRunning()) {
					log.info("IGV started.");
					return true;
				}
				log.info("Failed to connect to IGV.  Waiting 5 more seconds then trying again.");
				connectTry++;
			}
		} catch (IOException e) {
			log.error("Unable to start IGV: " + e.getMessage());
			return false;
		} catch (InterruptedException e) {
			
		} 
		log.error("Unable to start or connect to IGV");
		return false;
	}
*/
}