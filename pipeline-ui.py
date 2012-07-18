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
from time import sleep
from pipeline.utils import *
from pipeline.sample import Sample
from pipeline.runconfig import RunConfig
from pipeline.run import Run
from pipeline.chimera import Chimera
from pipeline.gast import Gast
from pipeline.validate import MetadataUtils
from pipeline.pipelinelogging import logger
from pipeline.trim_run import TrimRun
from pipeline.get_ini import readCSV


import logging
import argparse
from pipelineprocessor import process    

# pycogent
import cogent

import pipeline.constants as C

if __name__ == '__main__':
    THE_DEFAULT_BASE_OUTPUT = '.'

    usage = "usage: %prog [options] arg1 arg2"
    parser = argparse.ArgumentParser(description='MBL Sequence Pipeline')
    parser.add_argument('-c', '--configuration', required=True, dest = "configPath",
                                                 help = 'Configuration parameters of the run. See README File')
    parser.add_argument('-f', '--config_format',  required=False,   action="store",   default='csv', dest = "config_file_type",        
                                                 help = 'ini or csv') 
    parser.add_argument("-p", "--platform",     required=True,  action="store",   dest = "platform", 
                                                    help="Platform ")  
    parser.add_argument("-d", "--input_directory",     required=False,  action="store",   dest = "input_dir", default='./',
                                                    help="Directory where sequence files can be found. ")
    parser.add_argument("-r", "--run",     required=True,  action="store",   dest = "run", 
                                                    help="unique run number ")
    parser.add_argument("-ft", "--seq_file_type",     required=False,  action="store",   dest = "input_file_format", default='fasta',
                                                    help="Sequence file type: fasta, fastq or sff ")
    parser.add_argument("-fs", "--seq_file_suffix",     required=False,  action="store",   dest = "input_file_suffix", default='fa.unique',
                                                    help="Sequence file suffix [optional] ") 
    # see note for base_output_dir in runconfig.py  about line: 130                                               
    parser.add_argument("-b", "--baseoutputdir",     required=False,  action="store",  dest = "baseoutputdir", default='./',
                                                help="default: ./")
    parser.add_argument("-s", "--steps",     required=False,  action="store",   dest = "steps", default = 'status',
                                                help="Comma seperated list of steps.  Choices are: test,trim,chimera,status,upload_env454,gast,upload_vamps")
    parser.add_argument('-l', '--loglevel',  required=False,   action="store",   default='ERROR', dest = "loglevel",        
                                                 help = 'Sets logging level...INFO, DEBUG, [ERROR]') 

    
                                                 
    args = parser.parse_args() 
    print "\nLog Level set to:",args.loglevel.upper()
    # deal with logging level
    loggerlevel = logging.ERROR
    if args.loglevel.upper() == 'DEBUG':
        loggerlevel = logging.DEBUG
    elif  args.loglevel.upper() == 'INFO':     
        loggerlevel = logging.INFO
    logger.setLevel(loggerlevel)
    
    """
    TODO: read the config file here, depending on its type
    """
    if args.platform == 'illumina' and args.config_file_type == 'csv':
        v = MetadataUtils()
        # read the csv config file
        my_csv = readCSV(file_path = args.configPath)
        v.validate_illumina_csv(args, my_csv)
        
    elif args.platform == 'illumina' and args.config_file_type == 'ini':
        pass
    elif args.platform == '454' and args.config_file_type == 'csv':
        v = MetadataUtils()
        # read the csv config file
        my_csv = readCSV(file_path = args.configPath)
        v.validate_454_csv(args, my_csv)
        
    elif args.platform == '454' and args.config_file_type == 'ini':
        v = MetadataUtils()
        my_csv = readCSV(file_path = args.configPath)
        v.validate_454_ini(args, my_csv)
        
    elif args.platform == 'ion_torrent' and args.config_file_type == 'csv':
        pass
    elif args.platform == 'ion_torrent' and args.config_file_type == 'ini':
        pass
    else:
        sys.exit("Unknown platform and configFile type for validation")
    
    run = Run(args.configPath, args, os.path.dirname(os.path.realpath(__file__)))    
    
    cfg = None
    if my_csv:
        cfg = my_csv
    # now do all the work
    process(run, args.steps, cfg)

