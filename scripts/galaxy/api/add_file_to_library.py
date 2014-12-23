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

api_url = ''
api_key = ''
file = ''
library = ''
_debug = 0
   
def getFileMetaData(file_id):
    metadata = display(api_key, api_url + '/api/libraries/datasets/%s' % file_id, return_formatted = False)
    return metadata
    
def isFileInGalaxyLibrary(library_id, file):
    library_contents = display(api_key, api_url + "/api/libraries/%s/contents" % library_id, return_formatted = False)
    for entry in library_contents:
        if entry['type'] == 'file':
            libraryfile_metadata = getFileMetaData(entry['id'])
            libraryfile_fullpath = libraryfile_metadata['file_name']     
            if (file == libraryfile_fullpath):
                return 1
    return 0
    
""" 
The libraryName can either be just a library, or a library and folder path.
If the libraryName contains a folder path, we need to find the library, then
find the folder within the library to get its path.
"""    
def getGalaxyLibrary(libraryName):
    libs = display(api_key, api_url + '/api/libraries', return_formatted=False)
    for library in libs:
        if library['name'] == libraryName:
            return library
    print 'Unable to locate library %s in Galaxy' % libraryName
    return 0
    
def main():
    print 'Galaxy API URL: %s' % api_url
    print 'Galaxy API Key: %s' % api_key
    print 'File: %s' % file
    print 'Library: %s' % library
    print ''
    if not os.path.isfile(file):
        print 'Unable to location file'
        sys.exit(1)
    
    library_data = getGalaxyLibrary(library) 
    if library_data == 0:
        print '%s is not a library in Galaxy' % library
        sys.exit(1)
        
    if isFileInGalaxyLibrary(library_data['id'], file):
        print 'File already exists in Galaxy library'
        sys.exit(1)
    
    print 'Adding %s to %s' % (file, library)
    data = {}
    data['folder_id'] = library_data['root_folder_id']
    data['create_type'] = 'file'    
    data['file_type'] = 'auto'
    data['dbkey'] = ''
    data['upload_option'] = 'upload_paths'
    data['filesystem_paths'] = file
    data['link_data_only'] = 'link_to_files'
    
    libset = submit(api_key, api_url + "/api/libraries/%s/contents" % library_data['id'], data, return_formatted = True)
    print libset

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("api_key", help="API KEY")
    parser.add_argument('api_url', help='API URL')
    parser.add_argument('file', help="Full path of file to add")
    parser.add_argument('library', help="Galaxy library to add file to")
    args = parser.parse_args()
    api_key = args.api_key
    api_url = args.api_url
    file = args.file
    library = args.library
    main()

