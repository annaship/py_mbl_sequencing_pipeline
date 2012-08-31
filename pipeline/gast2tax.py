#!/usr/bin/env python
#!/usr/local/www/vamps/software/python/bin/python

##!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

import os
from stat import * # ST_SIZE etc
import sys
import shutil
import types

from pipeline.utils import *

from pipeline.gast import Gast

import argparse





import pipeline.constants as C


def run_gast2tax(args):
    import pipeline.gast

    if os.path.exists(args.names_file) and os.path.getsize(args.names_file) > 0:
        mygast = Gast()
        
        ( refdb, taxdb ) = mygast.get_reference_databases(args.dna_region)
        
        #print tax_file
        max_distance = C.max_distance['default']
        if args.dna_region in C.max_distance:
            max_distance = C.max_distance[args.dna_region] 
        
        
        ref_taxa = mygast.load_reftaxa(taxdb)

        mygast.assign_taxonomy(args.gast_dir, args.dna_region, args.names_file, ref_taxa)
        
        
        
        
if __name__ == '__main__':
    usage = """
        usage: ./pipeline-ui.py [options]
        
            options:
                -c/--configuration      configuration file with path  [required]
                
                -f/--config_format      configuration file format: csv or ini [optional (default:csv)]
                
                -p/--platform           Platform: illumina, 454 or ion_torrent [required]
                
                -i/--input_directory    Directory where sequence files can be found [optional (default: ./)]
                
               
             
        """
    if  len(sys.argv) == 1:
        print usage
        sys.exit()
    #THE_DEFAULT_BASE_OUTPUT = '.'

    # required items: configuration file, run and platform only
    # NO DEFAULTS HERE: DO Not give items defaults here as the script needs to look in the ini file as well
    # except steps (status) and loglevel (error) and 
    # see metadata.py and constants.py:  get_command_line_items()
    # BUT general section of ini file must have important things not supplied on command line
    # which means that csv file will require more commandline parameters.
    # NOTE: do not store any of the command line item as store_true or store_false or they
    # may not be able to be overridden buy the config file (ini).
    parser = argparse.ArgumentParser(description='MBL Sequence Pipeline: gast2tax')
    parser.add_argument('-dna', '--dna', required=True,                         dest = "dna_region",
                                                 help = 'Configuration parameters (.ini file) of the run. See README File')
    parser.add_argument("-max", "--max",     required=True,  action="store",              dest = "max_distance", 
                                                    help="unique run number ") 
                                                    
                                                    
    parser.add_argument("-o", "--output_dir",     required=True,  action="store",           dest = "gast_dir",            default = 'status',
                                                help="""
                                                """)
    parser.add_argument("-n", "--names_file",     required=True,  action="store",         dest = "names_file", 
                                                    help="Names file path ")                                              
   
     
    
   
    
    args = parser.parse_args() 
    
    
    run_gast2tax(args)

