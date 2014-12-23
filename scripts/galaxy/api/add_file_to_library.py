#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   12/23/14
This script adds a file to an existing library.
Usage: add_file_to_library.py <API_KEY> <API_URL> <FILE> <LIBRARY>

Algorithm:
1.  Make sure file exists and is readable
2.  Make sure file doesn't already exist in library
3.  Add to file to library as link

"""
import argparse
from string import split
from common import display, submit
import os
import sys
import time

_api_url = ''
_api_key = ''
_file = ''
_libraryPath = ''
_debug = 0
   
def isFileInGalaxyFolder(folder, file):
    folder_contents = display(_api_key, _api_url + "/api/folders/%s/contents" % folder['id'], return_formatted = False)
    for entry in folder_contents['folder_contents']:
        if entry['type'] == 'file':
            file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % entry['id'], return_formatted = False)
            file_fullpath = file_metadata['file_name']     
            if (file == file_fullpath):
                return 1
    return 0
   
""" 
The libraryName can either be just a library, or a library and folder path.
If the libraryName contains a folder path, we need to find the library, then
find the folder within the library to get its path.

Algorithm:
1.  Split up the libraryFolderPath to get the library and folder path as two 
    separate values
2.  Make sure the library exists and get its contents
3.  Find the folder in question and return it.  If it is the root folder, it
    will be "/", else the folder path will be what we are looking for.
""" 
def getGalaxyFolderFromLibrary(library, folderPath):
    library_contents = display(_api_key, _api_url + "/api/libraries/%s/contents" % library['id'], return_formatted = False)
    for entry in library_contents:
        if entry['name'] == folderPath:
            return entry

    print 'Unable to locate folder %s in library %s' % (folderPath, library['name'])
    sys.exit(1)
    
def getGalaxyLibrary(libraryName):  
    libs = display(_api_key, _api_url + '/api/libraries', return_formatted=False)
    for library in libs:
        if library['name'] == libraryName:
            return library
        
    print 'Unable to locate library %s in Galaxy' % libraryName
    sys.exit(1)
    
def main():
    print 'Galaxy API URL: %s' % _api_url
    print 'Galaxy API Key: %s' % _api_key
    print 'File: %s' % _file
    print 'Library: %s' % _libraryPath
    print ''
    if not os.path.isfile(_file):
        print 'Unable to location file'
        sys.exit(1)

    fields = split(_libraryPath, '/')
    libraryName = fields[0]
    if len(fields) == 1:
        folderPath = '/'
    else:
        sep = '/'
        folderPath = '/' + sep.join(fields[1:])
    
    library = getGalaxyLibrary(libraryName)
    folder = getGalaxyFolderFromLibrary(library, folderPath)
         
    if isFileInGalaxyFolder(folder, _file):
        print 'File already exists in Galaxy library'
        sys.exit(1)
    
    print 'Adding %s to %s' % (_file, _libraryPath)
    data = {}
    data['folder_id'] = folder['id']
    data['create_type'] = 'file'    
    data['file_type'] = 'auto'
    data['dbkey'] = ''
    data['upload_option'] = 'upload_paths'
    data['filesystem_paths'] = _file
    data['link_data_only'] = 'link_to_files'
    
    libset = submit(_api_key, _api_url + "/api/libraries/%s/contents" % library['id'], data, return_formatted = False)
    for lib in libset:
        file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)
        while file_metadata['state'] == 'running' or file_metadata['state'] == 'queued':
            print 'State is %s.  Sleep for 5 seconds.' % file_metadata['state'] 
            time.sleep(5)
            file_metadata = display(_api_key, _api_url + '/api/libraries/datasets/%s' % lib['id'], return_formatted = False)

        print 'State is %s' % file_metadata['state']

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("api_key", help="API KEY")
    parser.add_argument('api_url', help='API URL')
    parser.add_argument('file', help="Full path of file to add")
    parser.add_argument('library', help="Galaxy library to add file to")
    args = parser.parse_args()
    _api_key = args.api_key
    _api_url = args.api_url
    _file = args.file
    _libraryPath = args.library
    main()

