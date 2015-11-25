#!/usr/bin/env python


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



def run_gast2tax(args):

    #sys.path.append("/xraid2-2/vampsweb/"+args.site+"/")
    if args.site == 'vamps' or args.site == 'vampsdev' or args.site == 'new_vamps':
        sys.path.append('/groups/vampsweb')
    else:
        sys.path.append('/bioware/linux/seqinfo/bin/python_pipeline')
    #from pipeline.utils import *
    from py_mbl_sequencing_pipeline.pipeline.gast import Gast
    import py_mbl_sequencing_pipeline.pipeline.constants as C

    if os.path.exists(args.names_file) and os.path.getsize(args.names_file) > 0:
        
        class expando(object):pass
        runobj = expando()
        runobj.output_dir = ''
        runobj.use_cluster = True
        runobj.user = ''
        runobj.run = ''
        # for vamps this has to be True and empty list
        print 'VUU ',args.vamps_user_upload
        runobj.vamps_user_upload = args.vamps_user_upload
        runobj.platform = args.platform
        runobj.site = args.site
        runobj.project_dir = ''
        key = args.key
        
        runobj.datasets = []
        
        mygast = Gast(run_object = runobj)
        
        ( refdb, taxdb ) = mygast.get_reference_databases(args.dna_region)
        
        #print tax_file
        max_gast_distance = C.max_gast_distance['default']
        if args.dna_region in C.max_gast_distance:
            max_gast_distance = C.max_gast_distance[args.dna_region] 
        
        
        ref_taxa = mygast.load_reftaxa(taxdb)

        mygast.assign_taxonomy(key, args.gast_dir, args.dna_region, args.names_file, ref_taxa)
                                #key, gast_dir, dna_region, names_file, ref_taxa
        
        
        
if __name__ == '__main__':
    import argparse
    usage = """
        usage: ./gast2tax.py [options]
        
            options:
                -c/--configuration      configuration file with path  [required]
                
                -f/--config_format      configuration file format: csv or ini [optional (default:csv)]
                
                -p/--platform           Platform: illumina, 454 or ion_torrent [required]
                
                -i/--input_directory    Directory where sequence files can be found [optional (default: ./)]
                
               
             
        """
    if  len(sys.argv) == 1:
        print usage
        sys.exit()

    
    parser = argparse.ArgumentParser(description='MBL Sequence Pipeline: gast2tax')
    parser.add_argument('-dna', '--dna', required=True,                         dest = "dna_region",
                                                 help = 'Configuration parameters (.ini file) of the run. See README File')
    parser.add_argument("-max", "--max",     required=True,  action="store",              dest = "max_gast_distance", 
                                                    help="unique run number ") 
                                                    
                                                    
    parser.add_argument("-o", "--output_dir",     required=True,  action="store",           dest = "gast_dir",            default = 'status',
                                                help="""
                                                """)
    parser.add_argument("-n", "--names_file",     required=True,  action="store",         dest = "names_file", 
                                                    help="Names file path ")                                              
    parser.add_argument("-key", "--key",     required=True,  action="store",         dest = "key", 
                                                    help="key: dataset or lane_key") 
    parser.add_argument("-site", "--site",     required=True,  action="store",         dest = "site", 
                                                    help="vamps or vampsdev")  
    parser.add_argument("-platform", "--platform",     required=True,  action="store",         dest = "platform", 
                                                    help="454,illumina...") 
    parser.add_argument("-vamps_user_upload","--vamps_user_upload", required=False,  action="store_true",   dest = "vamps_user_upload", default=False,
                                                        help = 'data comes from vamps upload') 
    
    args = parser.parse_args() 
    
    
    run_gast2tax(args)

