#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   06/22/15
This script creates a library in Galaxy.
Usage: create_library.py [-h] [--api-key <API_KEY>] [--api-url <API_URL>] <library_name>

Algorithm:
  1.  Get a list of libraries
  2.  Check if the library already exists
  3.  If it doesn't, create the library
"""
import argparse
from string import split
from common import display, submit
import sys
import ConfigParser
import os

api_url = ''
api_key = ''
library_to_create = ''
_debug = 0
   
def main():
    if _debug == 1:
      print 'Galaxy API URL: %s' % api_url
      print 'Galaxy API Key: %s' % api_key
      print 'Library to create: %s' % library_to_create
      print ''

    if api_url == None or api_key == None:
        print "Galaxy API Key and/or URL was not specified"
        sys.exit(1)

    libs = display(api_key, api_url + '/api/libraries', return_formatted=False)
    for library in libs:
        if library['name'] == library_to_create and library['deleted'] == False:
            print 'Error: Library %s already exists.' % library['name']
            sys.exit(1)

    data = {}
    data['name'] = library_to_create
    
    result = submit(api_key, api_url + "/api/libraries", data, return_formatted = False)
    if not result['id'] == 0:
        print 'Success: Library created.'
    else:
        print 'Error: Failed to create library (%s).' % result['id']


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

    # Get library to create and override defaults with options specified on command lne
    parser = argparse.ArgumentParser()
    parser.add_argument('--api-url', help="Galaxy URL", default=api_url)
    parser.add_argument('--api-key', help="User's Galaxy Key", default=api_key)
    parser.add_argument('library', help="Name of data library to create")
    args = parser.parse_args()

    api_key = args.api_key
    api_url = args.api_url
    library_to_create = args.library
    main()
