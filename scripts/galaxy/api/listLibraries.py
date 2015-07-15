#!/usr/bin/env python
"""
Author: Ryan Golhar <ryan.golhar@bms.com>
Date: 07/15/15
This script list the available libraries in Galaxy
Usage: listLibraries.py [-h] [--api-key <API_KEY>] [--api-url <API_URL>] [--show-contents]
"""
import argparse
#from string import split
#import sys
from bioblend import galaxy
import ConfigParser
import os

api_url = ''
api_key = ''
show_contents = ''
_debug = 0

def main():
    if _debug == 1:
      print 'Galaxy API URL: %s' % api_url
      print 'Galaxy API Key: %s' % api_key
      print 'Show Library Contents: %s' % show_contents
      print ''

    if api_url == None or api_key == None:
        print "Galaxy API Key and/or URL was not specified"
        sys.exit(1)

    gi = galaxy.GalaxyInstance(url=args.api_url, key=args.api_key)
    galaxyLibraries = gi.libraries.get_libraries(deleted=False)
    for library in galaxyLibraries:
        if library['deleted'] == False:
            print "%s" % library['name']
            if show_contents == True:
                libraryContents = gi.libraries.show_library(library['id'], contents = True)
                for libraryEntry in libraryContents:
                    print "  %s" % libraryEntry['name']
    
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
    parser.add_argument('-l', '--show-contents', help="Show contents of libraries", action='store_true')
    args = parser.parse_args()

    api_key = args.api_key
    api_url = args.api_url
    show_contents = args.show_contents
    main()

