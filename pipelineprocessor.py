#!/usr/bin/env python

##!/usr/local/www/vamps/software/python/bin/python

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
from pipeline.vamps import Vamps
from pipeline.pipelinelogging import logger
from pipeline.trim_run import TrimRun

import logging
import json    
from pipeline.fasta_mbl_pipeline import MBLPipelineFastaUtils

TRIM_STEP = "trim"
CHIMERA_STEP = "chimera"
GAST_STEP = "gast"
UPLOAD_VAMPS_STEP  = "upload_vamps"
UPLOAD_ENV454_STEP = "upload_env454"
STATUS_STEP = 'status'

existing_steps = [ TRIM_STEP, CHIMERA_STEP, GAST_STEP, UPLOAD_VAMPS_STEP, UPLOAD_ENV454_STEP, STATUS_STEP ]

# the main loop for performing each of the user's supplied steps
def process(run, steps):
    # create output directory:
    requested_steps = steps.split(",")            
    
    if not os.path.exists(run.output_dir):
        logger.debug("Creating output directory: "+run.output_dir)
        os.makedirs(run.output_dir)      
                    
    # loop through official list...this way we execute the
    # users requested steps in the correct order                
    for step in requested_steps:
        if step not in existing_steps:
            print "Invalid processing step: " + step
            sys.exit()
        else:
            # call the method in here
            step_method = globals()[step]
            step_method(run)



    
    
# perform trim step
# TrimRun.trimrun() does all the work of looping over each input file and sequence in each file
# all the stats are kept in the trimrun object
#
# when complete...write out the datafiles for the most part on a lane/runkey basis
#
def trim(run):
    # (re) create the trim status file
    run.trim_status_file_h = open(run.trim_status_file_name, "w")
    
    # do the trim work
    mytrim = TrimRun(run) 
    
    # pass True to write out the straight fasta file of all trimmed non-deleted seqs
    # Remember: this is before chimera checking
    if run.platform == 'illumina':
        trim_codes = mytrim.trimrun_illumina(True)
    elif run.platform == '454':
        trim_codes = mytrim.trimrun_454(True)
    elif run.platform == 'ion-torrent':
    	trim_codes = mytrim.trimrun_ion_torrent(True)
    else:
        trim_codes = ['ERROR','No Platform Found']
        
    trim_results_dict = {}
    if trim_codes[0] == 'SUCCESS':
        # setup to write the status
        new_lane_keys = trim_codes[2]
        trim_results_dict['status'] = "success"
        trim_results_dict['new_lane_keys'] = new_lane_keys
        logger.debug("Trimming finished successfully")
        # write the data files
        mytrim.write_data_files(new_lane_keys)
        run.trim_status_file_h.write(json.dumps(trim_results_dict))
        run.trim_status_file_h.close()
    else:
        logger.debug("Trimming finished ERROR")
        trim_results_dict['status'] = "error"
        trim_results_dict['code1'] = trim_codes[1]
        trim_results_dict['code2'] = trim_codes[2]
        run.trim_status_file_h.write(json.dumps(trim_results_dict))
        run.trim_status_file_h.close()
        sys.exit()

# chimera assumes that a trim has been run and that there are files
# sitting around that describe the results of each lane:runkey sequences
# it also expectes there to be a trim_status.txt file around
# which should have a json format with status and the run keys listed        
def chimera(run):
    chimera_cluster_ids = [] 
    logger.debug("Starting Chimera Checker")
    # lets read the trim status file out here and keep those details out of the Chimera code
    new_lane_keys = convert_unicode_dictionary_to_str(json.loads(open(run.trim_status_file_name,"r").read()))["new_lane_keys"]
    mychimera = Chimera(run)
    c_den    = mychimera.chimera_denovo(new_lane_keys)
    if c_den[0] == 'SUCCESS':
        chimera_cluster_ids += c_den[2]
        chimera_code='PASS'
    elif c_den[0] == 'NOREGION':
        chimera_code='NOREGION'
    elif c_den[0] == 'FAIL':
        chimera_code = 'FAIL'
    else:
        chimera_code='FAIL'
    
    c_ref    = mychimera.chimera_reference(new_lane_keys)
    
    if c_ref[0] == 'SUCCESS':
        chimera_cluster_ids += c_ref[2]
        chimera_code='PASS'
    elif c_ref[0] == 'NOREGION':
        chimera_code = 'NOREGION'
    elif c_ref[0] == 'FAIL':
        chimera_code='FAIL'
    else:
        chimera_code='FAIL'
    
    #print chimera_cluster_ids
    run.chimera_status_file_h = open(run.chimera_status_file_name,"w")
    if chimera_code == 'PASS':  
        
        chimera_cluster_code = wait_for_cluster_to_finish(chimera_cluster_ids) 
        if chimera_cluster_code[0] == 'SUCCESS':
            logger.info("Chimera checking finished successfully")
            run.chimera_status_file_h.write("CHIMERA SUCCESS\n")
            
            
        else:
            logger.info("Chimera checking Failed")
            run.chimera_status_file_h.write("CHIMERA ERROR: "+str(chimera_cluster_code[1])+" "+str(chimera_cluster_code[2])+"\n")
            sys.exit()
            
    elif chimera_code == 'NOREGION':
        logger.info("No regions found that need chimera checking")
        run.chimera_status_file_h.write("CHIMERA CHECK NOT NEEDED\n")
        
    elif chimera_code == 'FAIL':
        logger.info("Chimera checking Failed")
        run.chimera_status_file_h.write("CHIMERA ERROR: \n")
        sys.exit()
    else:
        logger.info("Chimera checking Failed")
        run.chimera_status_file_h.write("CHIMERA ERROR: \n")
        sys.exit()
    sleep(2)   
    if  chimera_code == 'PASS' and  chimera_cluster_code[0] == 'SUCCESS':
        mychimera.write_chimeras_to_deleted_file(new_lane_keys)
        # should also recreate fasta
        # then read chimera files and place (or replace) any chimeric read_id
        # into the deleted file.
        
        mymblutils = MBLPipelineFastaUtils(new_lane_keys, mychimera.outdir)
        
        # write new cleaned files that remove chimera if apropriate
        # these are in fasta_mbl_pipeline.py
        # the cleaned file are renamed to the original name:
        # lane_key.unique.fa
        # lane_key.trimmed.fa
        # lane_key.names        -- 
        # lane_key.abund.fa     -- this file is for the uclust chimera script
        # lane_key.deleted.txt  -- no change in this file
        # THE ORDER IS IMPORTANT HERE:
        mymblutils.write_clean_fasta_file()
        mymblutils.write_clean_names_file()
        mymblutils.write_clean_uniques_file()
        mymblutils.write_clean_abundance_file()
        # write keys file for each lane_key - same fields as db table? for easy writing
        # write primers file for each lane_key
 
        
        # Write new clean files to the database
        # rawseq table not used
        # trimseq
        # runkeys
        # primers
        # run primers
        mymblutils.write_clean_files_to_database()

def gast(run):  
    
    mygast = Gast(run)
    
    # for vamps 'new_lane_keys' will be prefix 
    # of the uniques and names file
    # that was just created in vamps_gast.py
    if(run.vamps_user_upload):
        lane_keys = [run.user+run.runcode]        
    else:
        lane_keys = convert_unicode_dictionary_to_str(json.loads(open(run.trim_status_file_name,"r").read()))["new_lane_keys"]
        
    mygast.clustergast(lane_keys)
    sleep(5)
    mygast.gast_cleanup(lane_keys)
    sleep(5)
    mygast.gast2tax(lane_keys)

def upload_env454(run):
    print "TODO upload_env454(run)"
    # where are files?
    input_dir = run.input_dir
    print run.input_dir
    # config file has been validated and we know the data is there
    
    # is upload appropriate? that is, are the sequences trimmed?
    # how can I tell?
    # maybe there should be a status file
    # in the output_dir that receives updates during the entire run
    # or should that be in the db or both?
    # for a 454 run the data is in 20100917
    # the seqs are in 1_GATGA.trimmed.fa
    # and the trim data is in the run variable
    
    # presumably when illumina gets going the input_dir
    # will have a 'fa.unique' suffix (Meren's code will do this)
    
    
def upload_vamps(run):
    
    myvamps = Vamps(run)
    
    if(run.vamps_user_upload):
        lane_keys = [run.user+run.runcode]        
    else:
        lane_keys = convert_unicode_dictionary_to_str(json.loads(open(run.trim_status_file_name,"r").read()))["new_lane_keys"]
        
    myvamps.taxonomy(lane_keys)
    #myvamps.sequences(lane_keys)        
    #myvamps.exports(lane_keys)
    #myvamps.projects(lane_keys)
    #myvamps.info(lane_keys)
def status(run):
    print "TODO status(run)"
