#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   12/23/14
This script lists the libraries and contents of libraries available in Galaxy.
Usage: list_libraries.py <API_KEY> <API_URL> <DEBUG>

Algorithm:
  1.  Get a list of libraries
  2.  For each entry
  3.    If entry is a folder
  4.      Call listFolderAndContents(entry)
  5.    If entry is a file
  6.     List file name
"""
import argparse
from string import split
from common import display

api_url = ''
api_key = ''
_debug = 0
   
def main():
    print 'Galaxy API URL: %s' % api_url
    print 'Galaxy API Key: %s' % api_key
    print ''
    print 'Libraries'
    print '---------'
    i = 1;
    libs = display(api_key, api_url + '/api/libraries', return_formatted=False)
    for library in libs:
        # Print the name/description of the library
        if len(library['description']) != 0:
            print '%d. %s - %s' % (i, library['name'], library['description'])
        else:
            print '%d. %s' % (i, library['name'])
        if _debug:
            print '(%s)' % library

        # Print the Library Contents
        # Galaxy returns a list of files, folders, and files/folders within folders.  There is no tree of elements, just a list.
        library_contents = display(api_key, api_url + "/api/libraries/%s/contents" % library['id'], return_formatted = False) 
        for entry in library_contents:
            if entry['name'] == '/':
                continue
            if entry['type'] == 'folder':
                print ' %s' % entry['name']
                if _debug:
                    print ' (%s)' % entry
            if entry['type'] == 'file':
                fields = split(entry['name'], '/')
                print '    %s' % fields[len(fields) - 1]
                if _debug:
                    print '    (%s)' % entry
        print ''
        i+=1
        #for folder in folders:
        #    if folder['type'] == 'folder' and folder['name'] != "/":
	    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("api_key", help="API KEY")
    parser.add_argument('api_url', help='API URL')
    parser.add_argument('debug', help="Print Debug Statement (boolean)")
    args = parser.parse_args()
    _debug = args.debug
    api_key = args.api_key
    api_url = args.api_url
    main()

