#!/apps/sys/python/2.7.8/bin/python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   06/24/15
This script copies readme.txt,*.fastq.gz,*.fq.gz from a source locations to the local NGS library 
while maintaining directory structure.  This script also sets permissions such that data is correctly
copied.  It reads the readme.txt file in the directory to obtain information about the data.
This script sets the permission such that root:bioinfo is the owner of all files and the permissions
are go+rx.
Usage: copyFilesToLocalLibrary.py [-h] <path to directory to copy>
"""
import argparse
import fileinput
import logging
import os
import re
from shutil import copyfile

projectName = ''
diseaseArea = ''
RACode = ''
PI = ''
projectDescription = ''

def readReadmeTextFile(readMeTextFile):
    if not os.path.isfile(readMeTextFile):
        print "Cannot find %s.  Aborting." % readMeTextFile
        exit(-1)
    
    # Parse the readme.txt file for the following fields:
    # ProjectName:
    # DiseaseArea:
    # RACode:
    # PI:
    # ProjectDescription:
    global projectName
    global diseaseArea
    global RACode
    global PI
    global projectDescription
    for line in fileinput.input(readMeTextFile):
        m = re.match('(.+):(.+)', line)
        field = m.group(1)
        value = m.group(2).strip()
        if field == 'ProjectName':
            projectName = value
        elif field == 'DiseaseArea':
            diseaseArea = value
        elif field == 'RACode':
            RACode = value
        elif field == 'PI':
            PI = value
        elif field == 'ProjectDescription':
            projectDescription = value

    # Make sure all fields are properly entered
    if len(projectName) == 0:
        print "No project name specified in readme.txt"
        exit(-1)
    elif len(diseaseArea) == 0:
        print "No disease area specified in readme.txt"
        exit(-1)
    elif len(RACode) == 0:
        print "No RACode specified in readme.txt"
        exit(-1)
    elif len(PI) == 0:
        print "No PI specific in readme.txt"
        exit(-1)
    elif len(projectDescription) == 0:
        print "No Project Description in readme.txt"
        exit(-1)

def copyData():
    # Make destination directory if it doesn't already exist
    destDir = '%s/%s/%s/%s' % (args.destDir, diseaseArea, RACode, projectName)
    print "DestDir: %s" % destDir
    if not os.path.isdir(destDir):
        os.makedirs(destDir, 0755)

    # Copy readme.txt
    readMeTextFile = args.srcDir + '/readme.txt'
    destReadMeTxtfile = destDir + '/readme.txt'
    if os.path.exists(destReadMeTxtfile):
        print "%s already exists.  Skipping." % destReadMeTxtfile
    else:
        copyfile(readMeTextFile, destReadMeTxtfile)
        os.chmod(destReadMeTxtfile, 0444)

    copytree(args.srcDir, destDir)

def copytree(src, dst, symlinks=False, ignore=None):
    names = os.listdir(src)
    if ignore is not None:
        ignored_names = ignore(src, names)
    else:
        ignored_names = set()

    os.makedirs(dst)
    errors = []
    for name in names:
        if name in ignored_names:
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, symlinks, ignore)
            else:
                copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error) as why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except Error as err:
            errors.extend(err.args[0])
    try:
        copystat(src, dst)
    except WindowsError:
        # can't copy file access times on Windows
        pass
    except OSError as why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise Error(errors)

def main():
    print "srcDir: %s" % args.srcDir
    print "destDir: %s" % args.destDir
    print ""

    # Remove the trailing slash if one exists
    if args.srcDir.endswith('/'):
        args.srcDir = args.srcDir[:-1]
    if args.destDir.endswith('/'):
        args.destDir = args.destDir[:-1]

    # Make sure srcDir exists
    if not os.path.isdir(args.srcDir):
        print "%s is not a directory.  Please specify a directory" % args.srcDir
        exit(-1)
    
    # Read project specifics from readme.txt
    readReadmeTextFile(args.srcDir + '/readme.txt')

    # Copy Data
    copyData()

if __name__ == '__main__':
    # Set up logging
    logging.basicConfig(level=logging.DEBUG)
   
    # Parse Command-Line Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('srcDir', help="Source Directory to Copy")
    parser.add_argument('destDir', help="Destination Directory")
    args = parser.parse_args()

    # Do work
    main()
