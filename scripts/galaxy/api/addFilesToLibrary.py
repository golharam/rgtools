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

_debug = 1

'''
Upload/link a file into Galaxy.
Steps:
1.  Make sure file doesn't already exist.  If it does, skip it.
2.  Upload/Link to file
'''
def getGalaxyLibrary(libraryName):
    libs = display(api_key, api_url + '/api/libraries', return_formatted=False)
    for library in libs:
        if library['name'] == libraryName and library['deleted'] == False:
            return library

    print 'Unable to locate library %s in Galaxy' % libraryName
    sys.exit(1)

'''
1.  If the folder path exists, return it, else
2.  Get its parent folder then create the folder and return the folder
'''
def getGalaxyFolderFromFolderPath(galaxyFolderID, folderPath):
    # Base Case: Folder already exists, just return it
    library_contents = display(api_key, api_url + "/api/libraries/%s/contents" % galaxyLibrary['id'], return_formatted = False)
    for galaxyFolder in library_contents:
        if galaxyFolder['name'] == folderPath:
            return galaxyFolder

    # Iterative Case: 
    folders = folderPath.split('/')
    folder = folders.pop()
    if len(folders) == 0:
        galaxyParentFolder = getGalaxyFolderFromFolderPath(galaxyLibrary, '/')
    else:
        galaxyParentFolder = getGalaxyFolderFromFolderPath(galaxyLibrary, '/'.join(folders))
        
    # Create the folder
    print "Creating folder: %s" % folder
    data = {}
    data['encoded_parent_folder_id'] = galaxyParentFolder['id']
    data['name'] = folder
    result = submit(api_key, api_url + "/api/folders", data, return_formatted = False)
    if not result['id'] == 0:
        return result
    else:
        print "Unable to create folder: %s" % folder
        sys.exit(1)


def uploadFile(galaxyLibrary, filePath, fileDest):
    print "Uploading %s -> %s:%s" % (filePath, galaxyLibrary['name'], fileDest)

    # We have the library
    # We need to folder path from fileDest
    paths = fileDest.split('/')
    file = paths.pop()

    if len(paths) == 0:
        galaxyFolder = getGalaxyFolderFromFolderPath(galaxyLibrary, "/")
    else:
        folderPath = '/'.join(paths)
        galaxyFolder = getGalaxyFolderFromFolderPath(galaxyLibrary, folderPath)
    
    # For each directory in paths, find the corresponding folder in Galaxy.
    # If the folder doesn't exist, create it and continue
    #for path in paths:
    #    folder = getGalaxyFolder(path) 

    print 'Adding %s to %s:%s' % (file, galaxyLibrary['name'], galaxyFolder['name'])
#    data = {}
#    data['folder_id'] = folder['id']
#    data['create_type'] = 'file'
#    data['file_type'] = 'fastq'
#    data['dbkey'] = ''
#    data['upload_option'] = 'upload_paths'
#    data['filesystem_paths'] = os.path.abspath(filePath)
#    data['link_data_only'] = 'link_to_files'

#    libset = submit(_api_key, _api_url + "/api/libraries/%s/contents" % library['id'], data, return_formatted = False)
#    for lib in libset:
#        file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)
#        while file_metadata['state'] == 'running' or file_metadata['state'] == 'queued':
#            print 'State is %s.  Sleep for 5 seconds.' % file_metadata['state']
#            time.sleep(5)
#            file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)

#        print 'State is %s' % file_metadata['state']
    
def main_old():
    global dirPath
    

    # Make sure the directory exists of the files we are going to add.
    if not os.path.isdir(dirPath):
        print "%s is not a directory." % dirPath
        sys.exit(-1)
    else:
        # If there is not a trailing path, add it
        if not dirPath.endswith('/'):
            dirPath += '/'
    
    # Get the library and its contents
    print "Retrieving file list from Galaxy Library %s" % library
    galaxyLibrary = getGalaxyLibrary(library)
    galaxyLibraryContents = display(api_key, api_url + "/api/libraries/%s/contents" % galaxyLibrary['id'], return_formatted = False)

    # Scan diretory for fastq.gz files and upload the library
    print "Scanning for *.fastq.gz in %s" % dirPath
    filesToAdd = {}
    for root, dirs, files in os.walk(dirPath):
        #for dir in dirs:
        #    print os.path.join(root, dir)
        for file in files:
            if (file.endswith('.fastq.gz')):
                fileName = os.path.join(root, file)
                fileDest = fileName.split(dirPath)[1]

                # If the fileDest already exists in the library, skip it
                if fileExistsInGalaxy(fileDest, galaxyLibraryContents) == True:
                    print "%s already exists in Galaxy library.  Skipping."
                else:
                    filesToAdd[fileName] = fileDest

    print "Adding %d files." % len(filesToAdd)

    for filePath in filesToAdd:
        fileDest = filesToAdd[filePath]
        uploadFile(galaxyLibrary, filePath, fileDest)
          
'''
old
'''   
def uploadFileToFolder(galaxyFolderId, fileToUpload, fileDestination):
    print "Uploading file %s -> %s" % (fileToUpload, fileDestination)
    data = {}
    data['folder_id'] = galaxyFolderId
    data['create_type'] = 'file'
    data['file_type'] = 'fastq'
    data['upload_option'] = 'upload_paths'
    data['filesystem_paths'] = fileToUpload
    data['link_data_only'] = 'link_to_files'
    
#    libset = submit(args.api_key, args.api_url + "/api/uuu/%s/contents" % library['id'], data, return_formatted = False)
#    
#    for lib in libset:
#        file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)
#        while file_metadata['state'] == 'running' or file_metadata['state'] == 'queued':
#            print 'State is %s.  Sleep for 5 seconds.' % file_metadata['state']
#            time.sleep(5)
#            file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)
#
#        print 'State is %s' % file_metadata['state']
    
def doesFileExistsInGalaxy(fileDest, galaxyLibraryContents):
    for entry in galaxyLibraryContents:
        if entry['name'] == fileDest:
            return True
    return False
    
def main():
    if _debug == 1:
        print 'Galaxy API URL: %s' % api_url
        print 'Galaxy API Key: %s' % api_key
        print 'Library: %s' % args.library
        print 'File: %s' % args.file
        print ''

    # 1.  Get the Library
    galaxyLibrary = getGalaxyLibrary(args.library)
    
    # 2.  Make sure the file doesn't exist in the library
    galaxyLibraryContents = display(args.api_key, args.api_url + "/api/libraries/%s/contents" % galaxyLibrary['id'], return_formatted = False)
    # Get the base file name
    fileDest = os.path.basename(args.file)
    
    # 3.  Upload the file if it doesn't already exist in the library
    if doesFileExistsInGalaxy(fileDest, galaxyLibraryContents) == True:
        print "%s already exists in Galaxy library.  Skipping."
    else:
        uploadFileToFolder(galaxyLibrary['root_folder_id'], args.file, fileDest)



if __name__ == '__main__':
    # Get defaults from ~/.galaxy.ini
    config = ConfigParser.RawConfigParser()
    if os.path.exists(os.path.expanduser("~/.galaxy.ini")):
        config.read(os.path.expanduser("~/.galaxy.ini"))
        api_key = config.get('default', 'api_key')
        api_url = config.get('default', 'api_url')
    else:
        api_key = None
        api_url = None

    parser = argparse.ArgumentParser()
    parser.add_argument('--api-url', help="Galaxy URL", default=api_url)
    parser.add_argument('--api-key', help="User's Galaxy Key", default=api_key)
    parser.add_argument('library', help="Name of data library to add files to")
    parser.add_argument('file', help="File to upload")
    args = parser.parse_args()

    main()

