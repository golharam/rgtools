#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   06/24/15
This script create adds *.fastq.gz to a library while maintaining directory structure
Usage: addFilesToLibrary [-h] [--api-key <API_KEY>] [--api-url <API_URL>] <library_name> <path of directory to scan>

Algorithm:
  1.  Make sure library exists
  2.  Scan for and add files
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

def uploadFileToFolder(galaxyInstance, galaxyLibraryId, galaxyFolderId, fileToUpload, fileDestination):
    print "Uploading file %s -> %s" % (fileToUpload, fileDestination)
    
#    data = {}
#    data['folder_id'] = galaxyFolderId
#    data['source'] = 'admin_path'
#    data['link_data'] = True
#    data['preserve_dirs'] = True
#    data['file_type'] = 'fastq'
#    data['upload_dataset'] = fileToUpload
    
#    data['upload_option'] = 'upload_paths'
#    data['filesystem_paths'] = fileToUpload
#    data['link_data_only'] = 'link_to_files'
    libset = galaxyInstance.libraries.upload_from_galaxy_filesystem(galaxyLibraryId, fileToUpload, None, 'fastq', '', 'link_to_files')
#    libset = submit(args.api_key, args.api_url + "/api/libraries/datasets", data, return_formatted = True)
    
    for lib in libset:
        file_metadata = display(args.api_key, args.api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)
        while file_metadata['state'] == 'running' or file_metadata['state'] == 'queued':
            print 'State is %s.  Sleep for 5 seconds.' % file_metadata['state']
            time.sleep(5)
            file_metadata = display(args.api_key, args.api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)
        print 'State is %s' % file_metadata['state']
    
def doesFileExistsInGalaxy(fileDest, galaxyLibraryContents):
    # Because Galaxy strips the .gz from the file name, it will never match
    # We need to strip .gz to match
    fileDest = os.path.splitext(fileDest)[0]

    for entry in galaxyLibraryContents:      
        if entry['name'] == fileDest:
            return True
    return False
    
def main():
    if _debug == 1:
        print 'Galaxy API URL: %s' % args.api_url
        print 'Galaxy API Key: %s' % args.api_key
        print 'Library: %s' % args.library
        
    gi = galaxy.GalaxyInstance(url=args.api_url, key=args.api_key)
            
    # 1.  Get the Library
    galaxyLibrary = gi.libraries.get_libraries(name=args.library, deleted=False)[0]
    galaxyLibraryContents = display(args.api_key, args.api_url + "/api/libraries/%s/contents" % galaxyLibrary['id'], return_formatted = False)

    if os.path.isfile(args.pathToUpload):
        print 'File: %s' % args.pathToUpload
    
        # 2.  Make sure the file doesn't exist in the library
        # Get the base file name
        fileDest = '/' + os.path.basename(args.pathToUpload)
    
        # 3.  Upload the file if it doesn't already exist in the library
        if doesFileExistsInGalaxy(fileDest, galaxyLibraryContents) == True:
            print "%s already exists in Galaxy library.  Skipping." % fileDest
        else:
            uploadFileToFolder(gi, galaxyLibrary['id'], galaxyLibrary['root_folder_id'], args.pathToUpload, fileDest)
            
    elif os.path.isdir(args.pathToUpload):
        print 'Dir: %s' % args.pathToUpload
              
        # 3.  Scan the directory for *.fastq.gz and add each file 
        for root, dirs, files in os.walk(args.pathToUpload):
            #for dir in dirs:
            for file in files:
                if (file.endswith('.fq.gz')):
                    fileToUpload = os.path.join(root, file)
                    fileDest = '/' + os.path.basename(fileToUpload)
    
                    # If the fileDest already exists in the library, skip it
                    if doesFileExistsInGalaxy(fileDest, galaxyLibraryContents) == True:
                        print "%s already exists in Galaxy library.  Skipping." % fileDest
                    else:
                        # Upload the file, then refresh the library contents
                        uploadFileToFolder(gi, galaxyLibrary['id'], galaxyLibrary['root_folder_id'], fileToUpload, fileDest)
                        galaxyLibraryContents = display(args.api_key, args.api_url + "/api/libraries/%s/contents" % galaxyLibrary['id'], return_formatted = False)

                    
            

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

    parser = argparse.ArgumentParser()
    parser.add_argument('--api-url', help="Galaxy URL", default=_api_url)
    parser.add_argument('--api-key', help="User's Galaxy Key", default=_api_key)
    parser.add_argument('library', help="Name of data library to add files to")
    parser.add_argument('pathToUpload', help="File or Directory to upload")
    args = parser.parse_args()

    main()

