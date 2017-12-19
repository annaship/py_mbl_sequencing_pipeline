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
#sys.path.append("/bioware/pythonmodules/illumina-utils/")
sys.path.append("/xraid/bioware/linux/seqinfo/bin")
sys.path.append("/Users/ashipunova/bin/illumina-utils")
sys.path.append("/bioware/pythonmodules/illumina-utils/")


import shutil
import types
from time import sleep, time, gmtime, strftime
import timeit
from pipeline.utils import *
from pipeline.sample import Sample
from pipeline.runconfig import RunConfig
from pipeline.run import Run
from pipeline.chimera import Chimera
from pipeline.gast import Gast
from pipeline.vamps import Vamps
from pipeline.pipelinelogging import logger
from pipeline.trim_run import TrimRun
from pipeline.get_ini import readCSV
from pipeline.metadata import MetadataUtils
from pipeline.illumina_files import IlluminaFiles
from inspect import currentframe, getframeinfo

import logging
import json    
import IlluminaUtils.lib.fastalib as fastalib
#import fastalib as fa
from pipeline.fasta_mbl_pipeline import MBLPipelineFastaUtils
from pipeline.db_upload import MyConnection, dbUpload 
# from pipeline.utils import Dirs


# the main loop for performing each of the user's supplied steps
def process(runobj, steps):
    
    
    
    requested_steps = steps.split(",")            
    if 'clean' in requested_steps and len(requested_steps) > 1:
        sys.exit("The clean step cannot be combined with other steps - Exiting")
       
    # Open run STATUS File here.
    # open in append mode because we may start the run in the middle
    # say at the gast stage and don't want to over write.
    # if we re-run trimming we'll get two trim status reports
    runobj.run_status_file_h = open(runobj.run_status_file_name, "a")
    
    # loop through official list...this way we execute the
    # users requested steps in the correct order 

    for step in C.existing_steps:
        if step in requested_steps:
            # call the method in here
            print 'RUN',step
            step_method = globals()[step]
            step_method(runobj)

def validate(runobj):
    #open_zipped_directory(runobj.run_date, runobj.output_dir)
    #logger.debug("Validating")
    pass
    #v = MetadataUtils(run, validate=True)
    
    #print 'Validates:  Configfile and Run Object'
    #runobj.run_status_file_h.write(strftime("%Y-%m-%d %H:%M:%S", gmtime.time())+"\tConfigFile Validated\n")

    

##########################################################################################
# perform trim step
# TrimRun.trimrun() does all the work of looping over each input file and sequence in each file
# all the stats are kept in the trimrun object
#
# when complete...write out the datafiles for the most part on a lane/runkey basis
#
def trim(runobj):
    # def is in utils.py
    #open_zipped_directory(runobj.run_date, runobj.output_dir)
    # (re) create the trim status file
    runobj.trim_status_file_h = open(runobj.trim_status_file_name, "w")
    idx_keys = get_keys(runobj)
    
    # do the trim work
    mytrim = TrimRun(runobj, idx_keys) 
    
    # pass True to write out the straight fasta file of all trimmed non-deleted seqs
    # Remember: this is before chimera checking
    # trim_codes should alwas be a tuple with 3 elements!
    if runobj.vamps_user_upload:
        trim_codes = mytrim.trimrun_vamps(True)
    else:
        if runobj.platform == 'illumina':
            trim_codes = mytrim.filter_illumina()
    #        trim_codes = mytrim.trim_illumina(file_list = trim_codes[2])
        elif runobj.platform == '454':
            trim_codes = mytrim.trimrun_454(True)
        elif runobj.platform == 'ion-torrent':
            trim_codes = mytrim.trimrun_ion_torrent(True)        
        else:
            trim_codes = ('ERROR','No Platform Found','')
        
    trim_results_dict = {}
    if trim_codes[0] == 'SUCCESS':
        # setup to write the status
        new_lane_keys = trim_codes[2]
        trimmed_seq_count = trim_codes[1]
        if trimmed_seq_count == 0 or trimmed_seq_count == '0':
            trim_results_dict['status'] = "ERROR"
            logger.debug("Trimming finished: ERROR: no seqs passed trim")
        else:
            trim_results_dict['status'] = "success"
            logger.debug("Trimming finished successfully")
        
        trim_results_dict['new_lane_keys'] = new_lane_keys
        trim_results_dict['trimmed_seq_count'] = trimmed_seq_count
        
        # write the data files
        
        mytrim.write_data_files(new_lane_keys)
        runobj.trim_status_file_h.write(json.dumps(trim_results_dict)+"\n")
        runobj.trim_status_file_h.close()
        runobj.run_status_file_h.write(json.dumps(trim_results_dict)+"\n")
        runobj.run_status_file_h.close()
    else:
        logger.debug("Trimming finished ERROR")
        trim_results_dict['status'] = "ERROR"
        trim_results_dict['code1'] = trim_codes[1]
        trim_results_dict['code2'] = trim_codes[2]
        runobj.trim_status_file_h.write(json.dumps(trim_results_dict)+"\n")
        runobj.trim_status_file_h.close()
        runobj.run_status_file_h.write(json.dumps(trim_results_dict)+"\n")
        runobj.run_status_file_h.close()
        sys.exit("Trim Error")
        
        
    # def is in utils.py: truncates and rewrites
    #zip_up_directory(runobj.run_date, runobj.output_dir, 'w')

# chimera assumes that a trim has been run and that there are files
# sitting around that describe the results of each lane:runkey sequences
# it also expectes there to be a trim_status.txt file around
# which should have a json format with status and the run keys listed        
def chimera(runobj):
    chimera_cluster_ids = [] 
    logger.debug("Starting Chimera Checker")
    # lets read the trim status file out here and keep those details out of the Chimera code
    idx_keys = get_keys(runobj)
    #new_lane_keys = convert_unicode_dictionary_to_str(json.loads(open(runobj.trim_status_file_name,"r").read()))["new_lane_keys"]
    # Open run STATUS File here.
    # open in append mode because we may start the run in the middle
    # say at the gast stage and don't want to over write.
    # if we re-run trimming we'll get two trim status reports
    runobj.run_status_file_h = open(runobj.run_status_file_name, "a")
    
    mychimera = Chimera(runobj)
    logger.debug("\nStarting DeNovo Chimera")
    c_den    = mychimera.chimera_denovo()
    logger.debug("Ending DeNovo Chimera")
    if c_den[0] == 'SUCCESS':
        chimera_cluster_ids += c_den[2]   # add a list to a list
        logger.debug("chimera_cluster_ids: "+' '.join(chimera_cluster_ids))
        chimera_code='PASS'
    elif c_den[0] == 'NOREGION':
        chimera_code='NOREGION'
    elif c_den[0] == 'FAIL':
        chimera_code = 'FAIL'
    else:
        chimera_code='FAIL'

    logger.debug("Chimera DeNovo Code: "+chimera_code)
    logger.debug("\nStarting Reference Chimera")
    c_ref    = mychimera.chimera_reference()
    
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
    runobj.chimera_status_file_h = open(runobj.chimera_status_file_name,"w")
    if chimera_code == 'PASS':  
        
        if runobj.use_cluster:
            chimera_cluster_code = wait_for_cluster_to_finish(chimera_cluster_ids) 
            if chimera_cluster_code[0] == 'SUCCESS':
                logger.info("Chimera checking finished successfully")
                runobj.chimera_status_file_h.write("CHIMERA SUCCESS\n")
                runobj.run_status_file_h.write("CHIMERA SUCCESS\n")
                
            else:
                logger.info("3-Chimera checking Failed")
                runobj.chimera_status_file_h.write("3-CHIMERA ERROR: "+str(chimera_cluster_code[1])+" "+str(chimera_cluster_code[2])+"\n")
                runobj.run_status_file_h.write("3-CHIMERA ERROR: "+str(chimera_cluster_code[1])+" "+str(chimera_cluster_code[2])+"\n")
                sys.exit("3-Chimera checking Failed")
        else:
            chimera_cluster_code = ['SUCCESS','Not using cluster']
            logger.info("Chimera checking finished without using cluster")
            runobj.chimera_status_file_h.write("CHIMERA SUCCESS--no cluster\n")
            runobj.run_status_file_h.write("CHIMERA SUCCESS--no cluster\n")
    elif chimera_code == 'NOREGION':
        logger.info("No regions found that need chimera checking")
        runobj.chimera_status_file_h.write("CHIMERA CHECK NOT NEEDED\n")
        runobj.run_status_file_h.write("CHIMERA CHECK NOT NEEDED\n")
        
    elif chimera_code == 'FAIL':
        logger.info("1-Chimera checking Failed")
        runobj.chimera_status_file_h.write("1-CHIMERA ERROR: \n")
        runobj.run_status_file_h.write("1-CHIMERA ERROR: \n")
        sys.exit("1-Chimera Failed")
    else:
        logger.info("2-Chimera checking Failed")
        runobj.chimera_status_file_h.write("2-CHIMERA ERROR: \n")
        runobj.run_status_file_h.write("2-CHIMERA ERROR: \n")
        sys.exit("2-Chimera checking Failed")
    
    sleep(2) 
    
    if  chimera_code == 'PASS' and  chimera_cluster_code[0] == 'SUCCESS':
        logger.info("Writing Chimeras to deleted files")
        mychimera.write_chimeras_to_deleted_file()

        # should also recreate fasta
        # then read chimera files and place (or replace) any chimeric read_id
        # into the deleted file.
        
        mymblutils = MBLPipelineFastaUtils(idx_keys, runobj)
        
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

        #mymblutils.write_clean_files_to_database()
        
    # def is in utils.py: appends
    #zip_up_directory(runobj.run_date, runobj.output_dir, 'a')

def illumina_chimera(runobj):
    utils = PipelneUtils()

    start = time.time()
    mychimera = Chimera(runobj)
#     elapsed = (time.time() - start)
#     print elapsed
    print "Preparing input files (replacing \"frequency:\" with \";size=\" and capitalize reads)"

#     start = time.time()
#     mychimera.illumina_freq_to_size_in_chg()
#     elapsed = (time.time() - start)
#     print "1a) illumina_freq_to_size_in_chg time: %s" % elapsed
    start = time.time()
    mychimera.call_illumina_sed("from_frequency_to_size")
    elapsed = (time.time() - start)
    print "call_illumina_sed from_frequency_to_size time: %s" % elapsed
#     
    print "START chimera checking"
#     c_den = 
    mychimera.chimera_checking()
# #     print "c_den - check denovo res: %s" % c_den
#     print c_den
#     c_den = 
#     mychimera.chimera_checking("ref")
#     print c_den
#     todo: use run_until_done_on_cluster from utils
    """run after cluster is done with it work:"""
    start = time.time()  
    time_before = utils.get_time_now()
    print "time_before = %s" % time_before
    print "Waiting for the cluster..."
    while True:
        if utils.is_local():
            sleep(1)  
            break
      
        else:
            sleep(120)        
            cluster_done = mychimera.check_if_cluster_is_done(time_before)
            print "cluster_done = %s" % cluster_done
            if (cluster_done):
                break
    
    elapsed = (time.time() - start)
    print "Cluster is done with both chimera checkings in: %s" % elapsed     
    
    mychimera.check_if_chimera_dir_empty()
    
    mychimera.illumina_rm_size_files()

#     start = time.time()
#     mychimera.illumina_size_to_freq_in_chimer()
#     elapsed = (time.time() - start)
#     print "2a) illumina_size_to_freq_in_chimer time: %s" % elapsed
    start = time.time()
    mychimera.call_illumina_sed("from_size_to_frequency")
    elapsed = (time.time() - start)
    print "call_illumina_sed from_size_to_frequency time: %s" % elapsed
    
#     start = time.time()
#     print "Check chimeric statistics. If ref > 15% and ratio ref to de-novo > 2 use only de-novo"
#     mychimera.check_chimeric_stats()
#     elapsed = (time.time() - start)
#     print "check_chimeric_stats time: %s" % elapsed
    
    start = time.time()
    print "Creating nonchimeric files in %s" % mychimera.indir
    mychimera.move_out_chimeric()
    elapsed = (time.time() - start)
    print "move_out_chimeric time: %s" % elapsed
    
    
def illumina_chimera_after_cluster(runobj):  
    mychimera = Chimera(runobj)

    mychimera.illumina_rm_size_files()
    start = time.time()
    mychimera.illumina_size_to_freq_in_chimer()
    elapsed = (time.time() - start)
    print "illumina_size_to_freq_in_chimer time: %s" % elapsed
    
#     start = time.time()
#     print "Check chimeric statistics. If ref > 15% and ratio ref to de-novo > 2 use only de-novo"
#     mychimera.check_chimeric_stats()
#     elapsed = (time.time() - start)
#     print "check_chimeric_stats time: %s" % elapsed
    
    start = time.time()
    print "Creating nonchimeric files in %s" % mychimera.indir
    mychimera.move_out_chimeric()
    elapsed = (time.time() - start)
    print "move_out_chimeric time: %s" % elapsed
    print "illumina_chimera_after_cluster time = %s" % str(elapsed)

def illumina_chimera_only(runobj):  
    start = time.time()
    illumina_chimera(runobj)
    elapsed = (time.time() - start)
    print "illumina_chimera_only time = %s" % str(elapsed)

def illumina_files_demultiplex_only(runobj):  
    start = time.time()
    illumina_files = IlluminaFiles(runobj)
    illumina_files.open_dataset_files()
    illumina_files.split_files(compressed = runobj.compressed)
    elapsed = (time.time() - start)
    print "illumina_files demultiplex only time = %s" % str(elapsed)

def illumina_files(runobj):  
    utils = PipelneUtils()
    start = time.time()
#     illumina_files_demultiplex_only(runobj)
    illumina_files = IlluminaFiles(runobj)    
    if runobj.do_perfect: 
#         illumina_files.perfect_reads()
        script_file_name = illumina_files.merge_perfect()
        utils.run_until_done_on_cluster(script_file_name)
        
        script_file_name = illumina_files.trim_primers_perfect()
        utils.run_until_done_on_cluster(script_file_name)

    else:
#         illumina_files.partial_overlap_reads()
#         pass
# TODO: test utils.run_until_done_on_cluster(illumina_files.partial_overlap_reads_cluster())
        #TODO: add cutting to 251
        script_file_name = illumina_files.partial_overlap_reads_cluster()         
        utils.run_until_done_on_cluster(script_file_name)
        
        script_file_name = illumina_files.filter_mismatches_cluster()
        utils.run_until_done_on_cluster(script_file_name)
        
#         illumina_files.filter_mismatches()
#     illumina_files.uniq_fa()
    script_file_name = illumina_files.uniq_fa_cluster()
    utils.run_until_done_on_cluster(script_file_name)
#     illumina_chimera(runobj)
    elapsed = (time.time() - start)
    print "illumina_files time = %s" % str(elapsed)
        
def env454run_info_upload(runobj):
    my_read_csv = dbUpload(runobj)
    wrapped   = wrapper(my_read_csv.put_run_info)
    print "put_run_info time = %s" % timeit.timeit(wrapped, number=1) 
    
def get_sequences(my_env454upload, filenames):
    utils = PipelneUtils()

    sequences = [my_env454upload.make_seq_upper(filename) for filename in filenames]
    return utils.flatten_list_of_lists(sequences)

  
#     
# def env454upload(runobj):  
#     """
#     Run: pipeline dbUpload testing -c test/data/JJH_KCK_EQP_Bv6v4.ini -s env454upload -l debug
#     For now upload only Illumina data to env454 from files, assuming that all run info is already on env454 (run, run_key, dataset, project, run_info_ill tables)
#     Tables:
#     sequence_ill
#     sequence_pdr_info_ill
#     taxonomy
#     sequence_uniq_info_ill
#      
#     """
# 
#     full_upload = True
#     whole_start     = time.time()
# 
#     my_env454upload = dbUpload(runobj)
#     filenames       = my_env454upload.get_fasta_file_names()
#     if not filenames:
#         logger.debug("\nThere is something wrong with fasta files or their names, please check pathes, contents and suffixes in %s." % my_env454upload.fasta_dir)
#   
# #     sequences = get_sequences(my_env454upload, filenames)
#     for filename in filenames:
#         sequences = my_env454upload.make_seq_upper(filename)
#         env454upload_seq(my_env454upload, filename, sequences)
#         wrapped   = wrapper(my_env454upload.get_seq_id_dict, sequences)
#         get_seq_id_dict_time = timeit.timeit(wrapped, number=1)
#         logger.debug("get_seq_id_dict() took %s sec to finish" % get_seq_id_dict_time)
#        
#     total_seq = env454upload_all_but_seq(my_env454upload, filenames, full_upload)
#     my_env454upload.check_seq_upload()
#     logger.debug("total_seq = %s" % total_seq)
#     whole_elapsed = (time.time() - whole_start)
#     print "The whole upload took %s s" % whole_elapsed

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

def get_filename_base_no_suff(filename):
    filename_base_no_suff = "-".join(filename.split("/")[-1].split("-")[:-1])
    if (filename.find("MERGED") > 0):
#       For v4v5 illumia
        filename_base_no_suff   = "_".join(filename.split("/")[-1].split("_")[:3]).replace('MERGED', '1')
#     if len(filename_base_no_suff) == 0:
#         filename_base_no_suff   = "_".join(filename.split("/")[-1].split("_")[:3])
    return filename_base_no_suff 
# 
# def env454upload_no_seq(runobj):  
#     """
#     Run: pipeline dbUpload testing -c test/data/JJH_KCK_EQP_Bv6v4.ini -s env454upload -l debug
#     For now upload only Illumina data to env454 from files, assuming that all run info is already on env454 (run, run_key, dataset, project, run_info_ill tables) 
#     TODO: It's a dupliacate of env454upload! Needed refactoring.
#     """
#     
#     full_upload = False
#     
#     whole_start     = time.time()
# 
#     my_env454upload = dbUpload(runobj)
#     filenames       = my_env454upload.get_fasta_file_names()
#     if not filenames:
#         logger.debug("\nThere is something wrong with fasta files or their names, please check pathes, contents and suffixes in %s." % my_env454upload.fasta_dir)
#   
#     for filename in filenames:
#         sequences = my_env454upload.make_seq_upper(filename)
#         wrapped   = wrapper(my_env454upload.get_seq_id_dict, sequences)
#         get_seq_id_dict_time = timeit.timeit(wrapped, number=1)
#         logger.debug("get_seq_id_dict() took %s sec to finish" % get_seq_id_dict_time)
#         
#     total_seq = env454upload_all_but_seq(my_env454upload, filenames, full_upload)
#     my_env454upload.check_seq_upload()
#     logger.debug("total_seq = %s" % total_seq)
#     whole_elapsed = (time.time() - whole_start)
#     print "The whole_upload took %s s" % whole_elapsed

def env454upload(runobj):  
    full_upload = True
    env454upload_main(runobj, full_upload)

def env454upload_no_seq(runobj):  
    full_upload = False
    env454upload_main(runobj, full_upload)

def env454upload_main(runobj, full_upload):
    """
    Run: pipeline dbUpload testing -c test/data/JJH_KCK_EQP_Bv6v4.ini -s env454upload -l debug
    For now upload only Illumina data to env454 from files, assuming that all run info is already on env454 (run, run_key, dataset, project, run_info_ill tables)
    Tables:
    sequence_ill
    sequence_pdr_info_ill
    taxonomy
    sequence_uniq_info_ill

    """

    whole_start     = time.time()

    my_env454upload = dbUpload(runobj)
    filenames       = my_env454upload.get_fasta_file_names()
    if not filenames:
        logger.debug("\nThere is something wrong with fasta files or their names, please check pathes, contents and suffixes in %s." % my_env454upload.fasta_dir)

#     sequences = get_sequences(my_env454upload, filenames)
    for filename in filenames:
        sequences = my_env454upload.make_seq_upper(filename)
        if full_upload:
            env454upload_seq(my_env454upload, filename, sequences)
        wrapped   = wrapper(my_env454upload.get_seq_id_dict, sequences)
        get_seq_id_dict_time = timeit.timeit(wrapped, number=1)
        logger.debug("get_seq_id_dict() took %s sec to finish" % get_seq_id_dict_time)

    total_seq = env454upload_all_but_seq(my_env454upload, filenames, full_upload)
    my_env454upload.check_seq_upload()
    logger.debug("total_seq = %s" % total_seq)
    whole_elapsed = (time.time() - whole_start)
    print "The whole upload took %s s" % whole_elapsed
  
    
def env454upload_seq(my_env454upload, filename, sequences):
#     for filename in filenames:
    try:
        logger.debug("\n----------------\nfilename = %s" % filename)
#             sequences = my_env454upload.make_seq_upper(filename)
        if not (len(sequences)):
            logger.debug("There are 0 sequences in filename = %s" % filename)
#             continue           
        wrapped = wrapper(my_env454upload.insert_seq, sequences)
        insert_seq_time = timeit.timeit(wrapped, number=1)
        logger.debug("insert_seq() took %s sec to finish" % insert_seq_time)
    except:                       # catch everything
        print "\r[pipelineprocessor] Unexpected:"         # handle unexpected exceptions
        print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
        raise                       # re-throw caught exception   

# TODO: DRY
def env454upload_all_but_seq(my_env454upload, filenames, full_upload):
    insert_pdr_info_time = 0
    insert_taxonomy_time = 0
    insert_sequence_uniq_info_ill_time = 0
    total_seq = 0

    start_c = time.time()
    try:
        for filename in filenames:
            gast_dict             = my_env454upload.get_gast_result(os.path.basename(filename))
            
            filename_base_no_suff = get_filename_base_no_suff(filename)
            
            run_info_ill_id       = my_env454upload.get_run_info_ill_id(filename_base_no_suff)
            
            fasta                 = fastalib.SequenceSource(filename, lazy_init = False) 
            seq_in_file           = fasta.total_seq
            logger.debug("seq_in_file = %s" % seq_in_file)
            my_env454upload.put_seq_statistics_in_file(filename, seq_in_file)
            total_seq += seq_in_file

            start_fasta_next = time.time()
            
            start_prepare_pdf_info_query_time = 0
            start_prepare_pdf_info_query_time = time.time()
            all_insert_pdr_info_sql_to_run = my_env454upload.prepare_pdr_info_upload_query(fasta, run_info_ill_id, gast_dict)
            prepare_pdf_info_query_time = (time.time() - start_prepare_pdf_info_query_time)

            start_prepare_taxonomy_upload_query_time = 0
            start_prepare_taxonomy_upload_query_time = time.time()
            all_insert_taxonomy_sql_to_run = my_env454upload.prepare_taxonomy_upload_query(gast_dict)
            prepare_taxonomy_upload_query_time = (time.time() - start_prepare_taxonomy_upload_query_time)


            if (full_upload):
                insert_pdr_info_time = upload_w_time(my_env454upload, all_insert_pdr_info_sql_to_run)
                
            insert_taxonomy_time = upload_w_time(my_env454upload, all_insert_taxonomy_sql_to_run)
            
            start = time.time()
            my_env454upload.taxonomy.get_taxonomy_id_dict()
            elapsed = (time.time() - start)
            logger.debug("get_taxonomy_ids took %s sec to finish" % elapsed)

            prepare_insert_sequence_uniq_info_ill_sql_time = 0
            start_prepare_insert_sequence_uniq_info_ill_sql_time = time.time()
            all_insert_sequence_uniq_info_ill_sql_to_run = my_env454upload.prepare_insert_sequence_uniq_info_ill_sql(fasta, gast_dict)
            prepare_insert_sequence_uniq_info_ill_sql_time = (time.time() - start_prepare_insert_sequence_uniq_info_ill_sql_time)
            
            insert_sequence_uniq_info_ill_time = upload_w_time(my_env454upload, all_insert_sequence_uniq_info_ill_sql_to_run)

            logger.debug("start_fasta_loop took %s sec to finish" % (time.time() - start_fasta_next))
            logger.debug("prepare_pdf_info_query_time took %s sec to finish" % (prepare_pdf_info_query_time))
            logger.debug("start_prepare_taxonomy_upload_query_time took %s sec to finish" % (prepare_taxonomy_upload_query_time))
            logger.debug("insert_pdr_info() took %s sec to finish" % insert_pdr_info_time)
            logger.debug("insert_taxonomy_time.time() took %s sec to finish" % insert_taxonomy_time)
            
            logger.debug("prepare_insert_sequence_uniq_info_ill_sql_time took %s sec to finish" % (prepare_insert_sequence_uniq_info_ill_sql_time))
            
            logger.debug("insert_sequence_uniq_info_ill() took %s sec to finish" % insert_sequence_uniq_info_ill_time)
        logger.debug("env454upload_all_but_seq() took %s sec to finish" % (time.time() - start_c))
        return total_seq
        
    except:                       # catch everything
        print "\r[pipelineprocessor] Unexpected:"         # handle unexpected exceptions
        print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
        raise                       # re-throw caught exception             
      
def upload_w_time(my_env454upload, sql):
    start = time.time()
    my_env454upload.my_conn.cursor.execute(sql)
    my_env454upload.my_conn.cursor.execute("COMMIT")
    run_time = (time.time() - start)
    return run_time
    
def gast(runobj):  
    
    logger.info("STARTING GAST()")
#     logger.info("vsearch version: " % utils.get_vsearch_version)
    # for vamps 'new_lane_keys' will be prefix 
    # of the uniques and names file
    # that was just created in vamps_gast.py
    # or we can get the 'lane_keys' directly from the config_file
    # for illumina:
    # a unique idx_key is a concatenation of barcode_index and run_key
    # Should return a list not a string
    idx_keys = get_keys(runobj)
    
    # get GAST object
    mygast = Gast(runobj, idx_keys)
    
    
    # Check for unique files and create them if not there
    result_code = mygast.check_for_unique_files(idx_keys)
    runobj.run_status_file_h.write(json.dumps(result_code)+"\n")
    if result_code['status'] == 'ERROR':
        logger.error("uniques not found failed")
        sys.exit("uniques not found failed")
        if runobj.vamps_user_upload:
            write_status_to_vamps_db( runobj.site, runobj.run, "GAST ERROR", "uniques file not found - failed" )
    elif runobj.vamps_user_upload:
        write_status_to_vamps_db( runobj.site, runobj.run, result_code['status'], result_code['message'] )
        
    sleep(5)
    
    # CLUSTERGAST
    result_code = mygast.clustergast()
    runobj.run_status_file_h.write(json.dumps(result_code)+"\n")
    if result_code['status'] == 'ERROR':
        logger.error("clutergast failed")
        sys.exit("clustergast failed")
        if runobj.vamps_user_upload:
            write_status_to_vamps_db( runobj.site, runobj.run, "GAST ERROR", "clustergast failed" )
    elif runobj.vamps_user_upload:
        write_status_to_vamps_db( runobj.site, runobj.run, result_code['status'], result_code['message'] )
        
    sleep(5)
    
    # GAST_CLEANUP
    result_code = mygast.gast_cleanup()
    runobj.run_status_file_h.write(json.dumps(result_code)+"\n")
    if result_code['status'] == 'ERROR':
        logger.error("gast_cleanup failed")        
        sys.exit("gast_cleanup failed")
        if runobj.vamps_user_upload:
            write_status_to_vamps_db( runobj.site, runobj.run, "GAST ERROR", "gast_cleanup failed" )
    elif runobj.vamps_user_upload:
        write_status_to_vamps_db( runobj.site, runobj.run, result_code['status'], result_code['message'] )
        
    sleep(5)
    
    # GAST2TAX
    result_code = mygast.gast2tax()
    runobj.run_status_file_h.write(json.dumps(result_code)+"\n")
    if result_code['status'] == 'ERROR':
        logger.error("gast2tax failed") 
        sys.exit("gast2tax failed")
        if runobj.vamps_user_upload:
            write_status_to_vamps_db( runobj.site, runobj.run, "GAST ERROR", "gast2tax failed" )
    elif runobj.vamps_user_upload:
        write_status_to_vamps_db( runobj.site, runobj.run, result_code['status'], result_code['message'] )
        # write has_tax=1 to INFO-TAX.config
            
def cluster(runobj):
    """
    TO be developed eventually:
        Select otu creation method
        using original trimmed sequences
    """
    pass
    
    
def new_vamps(runobj):
    """
    
    """
    logger.info("STARTING NEW_VAMPS()")
    idx_keys = get_keys(runobj)
    myvamps = Vamps(runobj, idx_keys)
    myvamps.create_vamps_files()
    
    
def vampsupload(runobj):
    """
    Upload data files to VAMPS database
    """
    # for vamps 'new_lane_keys' will be prefix 
    # of the uniques and names file
    # that was just created in vamps_gast.py
    # or we can get the 'lane_keys' directly from the config_file
    # for illumina:
    # a unique idx_key is a concatenation of barcode_index and run_key
    idx_keys = get_keys(runobj)
    
#     if(runobj.vamps_user_upload):
#         idx_keys = [runobj.user+runobj.runcode]        
#     else:
#         idx_keys = convert_unicode_dictionary_to_str(json.loads(open(runobj.trim_status_file_name,"r").read()))["new_lane_keys"]
     
     # NOT NEEDED HERE: Find duplicate project names
     # if vamps user uploads this has already been done and this project is
     # already in vamps_upload_info table
     # if data from a csv file (illumina and 454) this also is not needed
     # as data is checked in metadata.py
    
     
    myvamps = Vamps(runobj, idx_keys)
    # Create files
    myvamps.create_vamps_files()
    # put files in db
    result_code = myvamps.load_vamps_db()
    
    if result_code[:5] == 'ERROR':
        logger.error("load_vamps_db failed") 
        sys.exit("load_vamps_db failed")
        if runobj.vamps_user_upload:
            write_status_to_vamps_db( runobj.site, runobj.run, "GAST_ERROR", result_code )
    elif runobj.vamps_user_upload:
        print "Finished loading VAMPS data",result_code
        write_status_to_vamps_db( runobj.site, runobj.run, 'GAST_SUCCESS', 'Loading VAMPS Finished' )
    # check here for completion of 
    # 1-file creation
    # 2-data appears in vamps
    
        
    
def status(runobj):
    
    f = open(runobj.run_status_file_name)
    lines = f.readlines()
    f.close()
    
    print "="*40
    print "STATUS LOG: "
    for line in lines:
        line =line.strip()
        print line
    print "="*40+"\n"
    
def clean(runobj):
    """
    Removes a run from the database and output directory
    """
    
    answer = raw_input("\npress 'y' to delete the run '"+runobj.run_date+"': ")
    if answer == 'y' or answer == 'Y':
        
        for (archiveDirPath, dirNames, fileNames) in os.walk(runobj.output_dir):
            print "Removing run:",runobj.run_date
            for file in fileNames:
                filePath = os.path.join(runobj.output_dir,file)
                print filePath
                os.remove(os.path.join(runobj.output_dir,file))
                # should we also remove STATUS.txt and *.ini and start again?
                # the directory will remain with an empty STATUS.txt file
                #os.removedirs(runobj.output_dir)

def get_keys(runobj):
    try:
        idx_keys = convert_unicode_dictionary_to_str(json.loads(open(runobj.trim_status_file_name,"r").read()))["new_lane_keys"]
        # {"status": "success", "new_lane_keys": ["1_GATGA"]}
    except:
        # here we have no idx_keys - must create them from run
        # if illumina they are index_runkey_lane concatenation
        # if 454 the are lane_key
        if runobj.vamps_user_upload:
            #print 'KEYS: '+' '.join(runobj.run_keys)
            idx_keys=runobj.samples.keys()
        else:
            if runobj.platform == 'illumina':  
                idx_keys = runobj.idx_keys
                ct = 0
                for h in runobj.samples:
                    logger.debug(h)
#                    logger.debug(h,runobj.samples[h]) #TypeError: not all arguments converted during string formatting
                    ct +=1
            elif runobj.platform == '454':
                idx_keys = runobj.idx_keys
            elif runobj.platform == 'ion_torrent':
                idx_keys = runobj.idx_keys
            else:
                logger.debug("GAST: No keys found - Exiting")
                runobj.run_status_file_h.write("GAST: No keys found - Exiting\n")
                sys.exit()
    if type(idx_keys) is types.StringType:
        return idx_keys.split(',')
    elif type(idx_keys) is types.ListType:
        return idx_keys
    else:
        return None
    return idx_keys

def vamps2upload(runobj):  
    full_upload = True
#     env454upload_main(runobj, full_upload)
    
    whole_start     = time.time()

    my_upload = dbUpload(runobj, db_server="vamps2")
    filenames       = my_upload.get_fasta_file_names()
    if not filenames:
        logger.debug("\nThere is something wrong with fasta files or their names, please check pathes, contents and suffixes in %s." % my_upload.fasta_dir)

#     sequences = get_sequences(my_upload, filenames)
    for filename in filenames:
        sequences = my_upload.make_seq_upper(filename)
        if full_upload:
            env454upload_seq(my_upload, filename, sequences)
        wrapped   = wrapper(my_upload.get_seq_id_dict, sequences)
        get_seq_id_dict_time = timeit.timeit(wrapped, number=1)
        logger.debug("get_seq_id_dict() took %s sec to finish" % get_seq_id_dict_time)

    total_seq = env454upload_all_but_seq(my_upload, filenames, full_upload)
    my_upload.check_seq_upload()
    logger.debug("total_seq = %s" % total_seq)
    whole_elapsed = (time.time() - whole_start)
    print "The whole upload took %s s" % whole_elapsed
  