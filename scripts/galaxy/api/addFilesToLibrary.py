#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   06/24/15
This script adds readme.txt, *.fastq.gz to a Galaxy Library
Usage: addFilesToLibrary [-h] [--api-key <API_KEY>] [--api-url <API_URL>] <path of directory to scan> <library_name>
"""
import ConfigParser
import os
import argparse
import sys
from common import display
from common import submit
import re
from bioblend import galaxy
import time

_debug = 1
    
def uploadFile(fileToUpload, galaxyInstance, galaxyLibrary, destFolder = '/'):        
    # Note: Right now, Galaxy strips .gz files of .gz.  So when searching of files, make sure to compare to data_set file_name

    libraryContents = galaxyInstance.libraries.show_library(galaxyLibrary['id'], contents = True)
    
    # Get the folder
    galaxyFolder_id = None
    for libraryEntry in libraryContents:
        if libraryEntry['name'] == destFolder and libraryEntry['type'] == 'folder':
            galaxyFolder_id = libraryEntry['id']
            break

    # Make sure the file doesn't exist in the destFolder
    for libraryEntry in libraryContents:
        if libraryEntry['type'] == 'file':            
            dataset = galaxyInstance.libraries.show_dataset(galaxyLibrary['id'], libraryEntry['id']) 
            print libraryEntry
            print dataset           
            if fileToUpload == dataset['file_name']:
                print "File already exists in library: %s" % libraryEntry['name']
                return
        
    # Upload file
    print "Uploading file %s -> %s:%s" % (fileToUpload, galaxyLibrary['name'], destFolder)
    result = galaxyInstance.libraries.upload_from_galaxy_filesystem(galaxyLibrary['id'], fileToUpload, galaxyFolder_id, file_type='fastq', link_data_only='link_to_files')
    print result
          
def main():
    if _debug == 1:
        print 'Galaxy API URL: %s' % args.api_url
        print 'Galaxy API Key: %s' % args.api_key
        print 'Path to Upload: %s' % args.pathToUpload
        print 'Library: %s' % args.library

    # 1.  Make sure Galaxy library exist
    # 2.  Scan path for readme.txt and *.fastq.gz and upload files to library

    # 1.
    gi = galaxy.GalaxyInstance(url=args.api_url, key=args.api_key)
    galaxyLibraries = gi.libraries.get_libraries(name=args.library, deleted=False)
    if len(galaxyLibraries) == 0:
        print "Unable to locate library %s" % args.library
        exit(-1)
    else:
        galaxyLibrary = gi.libraries.get_libraries(name=args.library, deleted=False)[0]
      
    # 2.  Scan the path for readme.txt, *.fastq.gz, *.fq.gz and upload to library
    if os.path.isfile(args.pathToUpload):
        uploadFile(args.pathToUpload, gi, galaxyLibrary)
            
    elif os.path.isdir(args.pathToUpload):
        if args.pathToUpload.endswith('/'):
            # Upload files in directory to dest
            # 3.  Scan the directory for *.fastq.gz and add each file 
            for root, dirs, files in os.walk(args.pathToUpload):
                #for dir in dirs:
                for file in files:
                    if (file.endswith('.gz')):
                        fileToUpload = os.path.join(root, file)
                        uploadFile(fileToUpload, gi, galaxyLibrary)
            
        else:
            # Upload directory and contents 
            print "make directory and upload to directory"
            
                            
if __name__ == '__main__':
    # Get defaults from ~/.galaxy.ini
    config = ConfigParser.RawConfigParser()
    if os.path.exists(os.path.expanduser("~/.galaxy.ini")):
        config.read(os.path.expanduser("~/.galaxy.ini"))
        _api_key = config.get('default', 'api_key')
        _api_url = config.get('default', 'api_url')
    else:
        _api_key = None
        _api_url = None

    # Parse Command-Line Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--api-url', help="Galaxy URL", default=_api_url)
    parser.add_argument('--api-key', help="User's Galaxy Key", default=_api_key)
    parser.add_argument('pathToUpload', help="File or Directory to upload")
    parser.add_argument('library', help="Name of Library to add data to")
    args = parser.parse_args()

    # Do work
    main()
