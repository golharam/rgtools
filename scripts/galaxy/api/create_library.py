#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date:   12/23/14
This script creates a library in Galaxy.
Usage: create_library.py <API_KEY> <API_URL> <library_name>

Algorithm:
  1.  Get a list of libraries
  2.  Check if the library already exists
  3.  If it doesn't, create the library
"""
import argparse
from string import split
from common import display, submit
import sys

api_url = ''
api_key = ''
library_to_create = ''
_debug = 0
   
def main():
    print 'Galaxy API URL: %s' % api_url
    print 'Galaxy API Key: %s' % api_key
    print 'Library to create: %s' % library_to_create
    print ''

    libs = display(api_key, api_url + '/api/libraries', return_formatted=False)
    for library in libs:
        if library['name'] == library_to_create:
            print 'Library already exists.'
            sys.exit(1)

    data = {}
    data['name'] = library_to_create
    
    result = submit(api_key, api_url + "/api/libraries", data, return_formatted = False)
    if not result['id'] == 0:
        print 'Library created.'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("api_key", help="API KEY")
    parser.add_argument('api_url', help='API URL')
    parser.add_argument('library', help="Library to create")
    args = parser.parse_args()
    api_key = args.api_key
    api_url = args.api_url
    library_to_create = args.library
    main()

