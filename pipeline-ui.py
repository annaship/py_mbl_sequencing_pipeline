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
from pipeline.pipelinelogging import logger
from pipeline.trim_run import TrimRun

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
    parser.add_argument('-f', '--config_format',  required=False,   action="store",   default='ini', dest = "config_file_type",        
                                                 help = 'ini or csv')                                              
    parser.add_argument("-b", "--baseoutputdir",     required=False,  action="store",  default=THE_DEFAULT_BASE_OUTPUT, dest = "baseoutputdirarg", 
                                                help="Comma seperated list of steps.  Choices are: trim,chimera,gast,vampsupload,all")
    
    parser.add_argument("-s", "--steps",     required=True,  action="store",   dest = "steps", 
                                                help="Comma seperated list of steps.  Choices are: trim,chimera,status,upload_env454,gast,upload_vamps")
    parser.add_argument('-l', '--loglevel',  required=False,   action="store",   default='ERROR', dest = "loglevel",        
                                                 help = 'Sets logging level...INFO, DEBUG, [ERROR]') 
    
                                                 
    args = parser.parse_args() 
    print "Log Level set to:",args.loglevel.upper()
    # deal with logging level
    loggerlevel = logging.ERROR
    if args.loglevel.upper() == 'DEBUG':
        loggerlevel = logging.DEBUG
    elif  args.loglevel.upper() == 'INFO':     
        loggerlevel = logging.INFO
    logger.setLevel(loggerlevel)
    
   
    run = Run(args.configPath, args.baseoutputdirarg, args.config_file_type, os.path.dirname(os.path.realpath(__file__)))  

    
    # read the config file
    
    # now do all the work
    process(run, args.steps)

