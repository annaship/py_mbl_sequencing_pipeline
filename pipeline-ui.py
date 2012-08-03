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
from pipeline.metadata import MetadataUtils
from pipeline.pipelinelogging import logger
from pipeline.trim_run import TrimRun
from pipeline.get_ini import readCSV

import logging
import argparse
from pipelineprocessor import process    

# pycogent
import cogent

import pipeline.constants as C

def get_values(is_vamps, args, config_dict={}):
    collector={}
#     requested_steps = args.steps.split(",") 
#     for step in collector_steps:
#         stepDict = 'pipeline_'+step+'_items'
#         step_collection = globals()[stepDict]
#         print stepDict, C.step_collection
        
    if is_vamps:
        stanza='vamps_pipeline'
        collector['vamps_user_upload'] = True
        for item in C.list_of_vamps_items:
            if item in args and getattr( args, item ) != None:
                collector[item]  = getattr( args, item )
            elif config_dict and item in config_dict[stanza] and config_dict[stanza][item] != '':
                collector[item]  = config_dict[stanza][item]
            elif config_dict and item in config_dict['general'] and config_dict['general'][item] != '':
                collector[item]  = config_dict['general'][item]
            else:
                collector[item] = C.list_of_vamps_items[item]
    else:
        stanza='mbl_pipeline'
        collector['vamps_user_upload'] = False
        for item in C.list_of_mbl_items:
            if item in args and getattr( args, item ) != None:
                collector[item]  = getattr( args, item )
            elif config_dict and item in config_dict[stanza] and config_dict[stanza][item] != '':
                collector[item]  = config_dict[stanza][item]
            elif config_dict and item in config_dict['general'] and config_dict['general'][item] != '':
                collector[item]  = config_dict['general'][item]
            else:
                collector[item] = C.list_of_mbl_items[item]
           
    return collector
    
def validate_args(args):
    """
    # THOUGHTS
    # vamps users
    # single project and dataset
    # Supply an ini file OR commandline (for web interface), but no csv file
    #
    # MBL pipeline
    # REQUIRE a csv file and a ini file
    """
    collector={}
    if args.csvPath:
        print "Must be MBL origin"
        
        if not args.configPath:
            sys.exit("MBL Pipeline: you must supply an ini file with a csv file")
        else:
            config_dict = configDictionaryFromFile_ini(args.configPath)           
            collector= get_values(False, args, config_dict)
    else:
        print "VAMPS User origin"
        list_of_items = C.list_of_vamps_items
        
        if args.configPath:
            config_dict = configDictionaryFromFile_ini(args.configPath)
            # eg dna_region
            # precidence: cl,ini file
            
            collector = get_values(True, args, config_dict)
            
        else:
            # Should never get here because configPath is required
            sys.exit("No config file")
        collector['project'] = collector['project'][:1].capitalize() + collector['project'][1:]
        # these are all the bool items from the ini file
        # they need to be converted fron str to bool here
        if 'from_fasta' in collector and collector['from_fasta'] == 'True':
            collector['from_fasta'] = True
        else:
            collector['from_fasta'] = False
        if 'use_cluster' in collector and collector['use_cluster'] == 'True':
            collector['use_cluster'] = True
        else:
            collector['use_cluster'] = False
        if 'load_vamps_database' in collector and collector['load_vamps_database'] == 'True':
            collector['load_vamps_database'] = True
        else:
            collector['load_vamps_database'] = False 
            
        collector['runcode'] = args.run
        collector['run'] = args.run
        collector['steps'] = args.steps
#     if 'run' in collector and collector['run'] != '':
#         collector['runcode'] = collector['run']
#     elif 'runcode' in collector and collector['runcode'] != '':
#         collector['run'] = collector['runcode']
#     else:
#         sys.exit("No run or runcode found")
    
    collector['loglevel'] = collector['loglevel'].upper()
    collector['datetime'] = str(datetime.date.today())
    return collector
        
        
if __name__ == '__main__':
    usage = """
        usage: ./pipeline-ui.py [options]
        
            options:
                -c/--configuration      configuration file with path  [required]
                
                -f/--config_format      configuration file format: csv or ini [optional (default:csv)]
                
                -p/--platform           Platform: illumina, 454 or ion_torrent [required]
                
                -i/--input_directory    Directory where sequence files can be found [optional (default: ./)]
                
                -r/--run                Run - number or date  [required]
                
                -ft/--seq_file_type     File type for sequences: fasta, fastq or sff 
                                            [optional (default: fasta)]
                                            
                -fs/--seq_file_suffix   File suffix - useful when there are additional files
                                            in the input directory that you don't want to include. [optional (default: fa.unique)]
                                            
                -b/--baseouputdir       Base output directory where the run directory will be found.
                                            The run directory will be created if it is not found.  [optional (default: ./)]
                                            
                -s/--steps              Steps to be performed by this pipeline (comma separated list)
                                            Choices:    validate        - validates your metadata file
                                                        status          - prints out status messages if any
                                                        trim            - trims your sequences
                                                        chimera         - performs chimera check on trimmed sequences
                                                        upload_env454   - Load data into the env454 database
                                                        gast            - assign taxonomy to the trimmed sequences using GAST 
                                                        upload_vamps    - load sequences and taxonomy to VAMPS
                                                        clean           - removes run from database and filesystem
                                                        
                -l/--loglevel           Change the level of logging: info, debug, error   [optional (default: error)]
             
        """
    if  len(sys.argv) == 1:
        print usage
        sys.exit()
    #THE_DEFAULT_BASE_OUTPUT = '.'

    # required items: configuration file, run and config_file_type only
    # DO Not give items defaults here as the script needs to look in the ini file as well
    # except steps (status) and loglevel (error) and baseoutputdir (./)(for ouptput_dir creation)
    # see validate.py:  get_command_line_items()
    # BUT general section of ini file must have important things not supplied on command line
    # which means that csv file will require more commandline parameters.
    parser = argparse.ArgumentParser(description='MBL Sequence Pipeline')
    parser.add_argument('-c', '--configuration', required=True,                         dest = "configPath",
                                                 help = 'Configuration parameters (.ini file) of the run. See README File')
    parser.add_argument("-r", "--run",     required=True,  action="store",              dest = "run", 
                                                    help="unique run number ") 
                                                    
                                                    
    parser.add_argument("-s", "--steps",     required=False,  action="store",           dest = "steps",            default = 'status',
                                                help="""
                                                Comma seperated list of steps.  
                                                Choices are: validate,trim,chimera,status,upload_env454,gast,otu,upload_vamps,clean
                                                """)
                                                
    #################################################################################################################### 
    parser.add_argument('-l', '--loglevel',  required=False,   action="store",          dest = "loglevel",          default='ERROR',       
                                                 help = 'Sets logging level... DEBUG, [INFO], WARNING, ERROR, CRITICAL')
     # see note for base_output_dir in runconfig.py  about line: 130                                               
    parser.add_argument("-b", "--baseoutputdir",     required=False,  action="store",   dest = "baseoutputdir", default='.',
                                                help="default: ./") 
    parser.add_argument("-i", "--input_directory",     required=False,  action="store", dest = "input_dir",   default='',
                                                    help="Directory where sequence files can be found. ")                           
    #################################################################################################################### 
    # Illumina and 454 Specific
    parser.add_argument('-csv', '--csv',            required=False,                         dest = "csvPath",
                                                        help = 'CSV file path. See README File')
    parser.add_argument("-p", "--platform",     required=False,  action="store",         dest = "platform", 
                                                    help="Platform ")                                                  
    parser.add_argument('-f', '--config_format',  required=False,   action="store",     dest = "config_file_type",  
                                                 help = 'ini or csv') 
    
    
    
    parser.add_argument("-ft", "--seq_file_type",     required=False,  action="store",  dest = "input_file_format", default='',
                                                    help="Sequence file type: fasta, fastq or sff ")
    parser.add_argument("-fs", "--seq_file_suffix",     required=False,  action="store",dest = "input_file_suffix", default='',
                                                    help="Sequence file suffix [optional] ") 
    
    
     
    parser.add_argument('-cp', '--compressed',  required=False,   action="store",       dest = "compressed",        default='True',       
                                                 help = 'Make it "False" if illumina fastq files are not compressed with gzip') 
    parser.add_argument('-o', '--output_directory',  required=False,   action="store",  dest = "output_dir",        default='.',       
                                                 help = 'Output directory') 
    parser.add_argument('-h', '--database_host',  required=False,   action="store",  dest = "database_host",        default='vampsdev',       
                                                 help = 'Database host') 
    parser.add_argument('-db', '--database_name',  required=False,   action="store", dest = "database_name",        default='test',       
                                                 help = 'Database name') 
    #
    # VAMPS Specific: all can be in the ini file
    #
    parser.add_argument("-site",  "--site",         required=False,  action="store",   dest = "site", 
                                                        help="""database hostname: vamps or vampsdev
                                                        [default: vampsdev]""")     
    parser.add_argument("-u", "--user",             required=False,  action="store",   dest = "user", 
                                                        help="user name")         
    parser.add_argument("-p", "--project",          required=False,  action='store', dest = "project", 
                                                        help="") 
    parser.add_argument('-d',"--dataset",           required=False,  action="store",   dest = "dataset", 
                                                        help = '')
    parser.add_argument("-load", "--load_database", required=False,  action="store",   dest = "load_db", 
                                                        help = 'VAMPS: load files into vamps db')                                              
    parser.add_argument("-env", "--envsource",      required=False,  action="store",   dest = "env_source_id", 
                                                        help = '')
    parser.add_argument("-uc", "--use_cluster",      required=False,  action="store",   dest = "use_cluster", 
                                                        help = 'if false: the cluster will not be used (for testing)') 
    #DEBUG	Detailed information, typically of interest only when diagnosing problems.
    #INFO	Confirmation that things are working as expected.
    #WARNING	An indication that something unexpected happened, or indicative of some problem in the near future (e.g. 'disk space low'). 
    #           The software is still working as expected.
    #ERROR	Due to a more serious problem, the software has not been able to perform some function.
    #CRITICAL	A serious error, indicating that the program itself may be unable to continue running.
    
    args = parser.parse_args() 

    # this will read the args and ini file and return a dictionary
    data_object = validate_args(args)
    if 'commandline' in data_object and data_object['commandline'] == True:
        for item in data_object:
            print item+' = ',data_object[item]
        answer = raw_input("\n\tDoes this look okay? ('c' to continue; 'q' to quit) ")
        if answer == 'q':
            sys.exit()
        elif answer == 'c':
            pass
        else:
            sys.exit()
            
    # set logging
    
    print "\nLog Level set to:",data_object['loglevel']    
    logger.setLevel(data_object['loglevel'] )
    
    logger.info("Starting pipeline")
    ##############
    #
    #  Test cl parameters
    #
    ############## 
    # CL RULES:
    # for ini file:  (no plurals)
    # 1) CL: input_dir ONLY shall be supplied on CL - no input filenames
    #   
    # 2) All input files should be in the same directory AND of the same format
    #       
    # 3) Supply a input_file_suffix on the CL if there are varying file types in the
    #       input_dir and you only are using some (default will read all files)
    # 4) 
    #
   
    ##############
    #
    # CREATE or FIND OUTPUT DIRECTORY
    # need to look for or create output_dir here
    # base output directory and run are required so need to create output_dir here
    # to write ini file and status file
    ##############
    if not os.path.exists(os.path.join(args.baseoutputdir,args.run)):
        logger.debug("Creating output directory: "+os.path.join(args.baseoutputdir,args.run))
        os.makedirs(os.path.join(args.baseoutputdir,args.run))    
    
   
    
    ##############
    #
    #  VALIDATE THE INI FILE
    #
    ##############  
    # pass True to validate, convert (csv => ini if needed) and write out ini file
    #print 'steps',args.steps
    #if args.steps == 'validate':
    #    v = MetadataUtils(args, validate=True)
    #else:
    v = MetadataUtils(args)
    v.convert_and_save_ini()
    v.validate()
    general_data = v.get_general_data()
    
    answer = v.get_confirmation(args.steps, general_data)
    
    if answer == 'q':
        sys.exit()
    elif answer == 'v':
        # view CONFIG file contents
        fh = open(os.path.join(args.baseoutputdir, args.run, args.run+'.ini'))
        lines = fh.readlines()
        print "\n=== START ===\n"
        for line in lines:
            line = line.strip()
            print line
        print "==== END ====\n"
        sys.exit()
    elif answer != 'c':
        sys.exit()
    ##############
    #
    # CREATE THE RUN OBJECT (see runconfig.py for details)
    #
    ##############     
    run = Run(args.configPath, args, os.path.dirname(os.path.realpath(__file__)))    
    

#    for key in run.samples:
#        print key,run.samples[key].dataset
#    sys.exit()
    ##############
    #
    # now do all the work
    #
    ##############         
    process(run, args.steps)

