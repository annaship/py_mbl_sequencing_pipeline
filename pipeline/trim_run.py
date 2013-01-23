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
import sys, os, stat
import shutil
#import hashlib
sys.path.append("/bioware/linux/seqinfo/bin/")
sys.path.append("/bioware/pythonmodules/illumina-utils/")
sys.path.append("/Users/ashipunova/bin/illumina-utils")
sys.path.append('/bioware/linux/seqinfo/bin/python_pipeline/py_mbl_sequencing_pipeline')

from suites.primer import PrimerSuite 
from pipeline.primer_utils import *
from pipeline.utils import *
from fastalib import *
from fastqlib import *
from pipeline.Fasta import sfasta
from pipeline.anchortrimming_mbl import *
from pipeline.utils import Dirs, PipelneUtils
from pipeline.pipelinelogging import logger
from Bio import SeqIO
import subprocess


DELETED_RUNKEY      = 'Runkey: no key found'
DELETED_LANEKEY     = 'Lanekey: unknown lane/key'
DELETED_PROXIMAL    = 'Proximal: no proximal primer found'
DELETED_DISTAL      = 'Distal: no distal primer found'
DELETED_N           = 'N'
DELETED_QUALITY     = 'Quality'
DELETED_NO_INSERT   = 'No Insert'
DELETED_MINIMUM_LENGTH      = 'Minimum Length'
DELETED_MAXIMUM_LENGTH      = 'Maximum Length'

class TrimRun( object ):
    
    """ Define here"""
    def __init__(self, run_object = None, idx_keys=None):
    
        self.runobj     = run_object
        self.outdir     = self.runobj.output_dir
        
        self.indir     = self.runobj.input_dir
        
        if self.runobj.vamps_user_upload:
            site = self.runobj.site
            dir_prefix=self.runobj.user+'_'+self.runobj.run+'_trim'
        else:
            site = ''
            dir_prefix = self.runobj.run
            
        dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, site = site) 
        
        self.analysis_dir = dirs.check_dir(dirs.analysis_dir)        
        self.trimming_dir = dirs.check_dir(dirs.trimming_dir)
        
        
        # do something with 'run'.
        self.run = self.runobj.run
        logger.debug("Rundate:" + str(self.run))
        self.use_cluster = self.runobj.use_cluster
        self.idx_keys = idx_keys
            
        self.seqDirs            = {}
        self.dna_regions        = {}
        self.taxonomic_domain   = {}
        self.expanded_primers   = {}
        self.proximal_primers   = {}
        self.distal_primers     = {}
        self.anchor_name            = {}
        self.adtnl_anchors ={}
        self.anchors ={}
        self.psuite             = {}
        self.id_list_passed     = {}
        self.deleted_ids        = {}
        self.trimmed_ids        = {}
        self.uniques            = {}
        self.names              = {}
        self.uniquefa ={}
        self.abundfa ={}
        self.fa ={}
        self.statsFileName     = 'run_trim_stats'
        
        os.environ['SGE_ROOT']='/usr/local/sge'
        os.environ['SGE_CELL']='grendel'
        path = os.environ['PATH']
        os.environ['PATH'] = '/usr/local/sge/bin/lx24-amd64:'+path
        
        self.runbin={}
#         print 'SAMPLES: '
#         for idx_key in self.runobj.samples:
#                 print idx_key,self.runobj.samples[idx_key]
#         print 'SAMPLES: '    
        if self.runobj.vamps_user_upload:
            for idx_key in self.runobj.samples:
                sample = self.runobj.samples[idx_key]
                # strip off surrounding single or double quotes
                self.seqDirs[idx_key]          = sample.direction
                self.dna_regions[idx_key]      = sample.dna_region
                self.taxonomic_domain[idx_key] = sample.taxonomic_domain
                
                # this should be defaiult 'suite' of anchors
                # but also in ini file: anchor=XXXXX will be added to defaults
                self.anchor_name[idx_key]         = sample.anchor
                self.adtnl_anchors[idx_key]         = sample.stop_sequences  #list
               
                self.anchors[idx_key]           = {}
                
                #self.id_list_all[idx_key]      = []
                self.id_list_passed[idx_key]    = []
                self.deleted_ids[idx_key]       = {}
                self.deleted_ids['nokey']       = {}
                self.trimmed_ids[idx_key]       = {}
                self.uniques[idx_key]           = {}
                self.names[idx_key]             = {}
                self.fa[idx_key]       = FastaOutput(os.path.join(self.trimming_dir, idx_key) + ".trimmed.fa")
                
      
      
                #####################
                #
                #  PrimerSuite Class
                #
                #####################
                print self.taxonomic_domain[idx_key],self.dna_regions[idx_key],idx_key
                
                self.psuite[idx_key] = PrimerSuite(self.runobj, self.taxonomic_domain[idx_key],self.dna_regions[idx_key],idx_key)
                #self.runbin['psuite'][idx_key]= PrimerSuite(self.taxonomic_domain[idx_key],self.dna_regions[idx_key])
                
                if(self.seqDirs[idx_key] == 'F' or self.seqDirs[idx_key] == 'B'):
                    self.proximal_primers[idx_key] = self.psuite[idx_key].primer_expanded_seq_list['F']
                    self.distal_primers[idx_key]   = self.psuite[idx_key].primer_expanded_seq_list['R']
                    if self.anchor_name[idx_key]:
                        self.anchors[idx_key] = get_anchor_list(self.runobj, self.anchor_name[idx_key], self.adtnl_anchors[idx_key])
     
                    
                if(self.seqDirs[idx_key] == 'R' or self.seqDirs[idx_key] == 'B'):
                    self.proximal_primers[idx_key] = [revcomp(primer_seq) for primer_seq in self.psuite[idx_key].primer_expanded_seq_list['R'] ]  
                    self.distal_primers[idx_key]   = [revcomp(primer_seq) for primer_seq in self.psuite[idx_key].primer_expanded_seq_list['F'] ] 
                    if self.anchor_name[idx_key]:
                        self.anchors[idx_key] = [revcomp( anchor ) for anchor in get_anchor_list(self.runobj, self.anchor_name[idx_key], self.adtnl_anchors[idx_key]) ]
                
                    if len(self.proximal_primers[idx_key]) == 0 and len(self.distal_primers[idx_key]) == 0:
                        logger.debug("**** Didn't find any primers that match any of the domain/regions in the lane/key sections")
                
        else:            
            if self.runobj.platform == '454':
                for idx_key in self.idx_keys:
                    
                    sample = self.runobj.samples[idx_key]
                    # strip off surrounding single or double quotes
                    self.seqDirs[idx_key]          = sample.direction
                    self.dna_regions[idx_key]      = sample.dna_region
                    self.taxonomic_domain[idx_key] = sample.taxonomic_domain
                    
                    # this should be defaiult 'suite' of anchors
                    # but also in ini file: anchor=XXXXX will be added to defaults
                    self.anchor_name[idx_key]         = sample.anchor
                    self.adtnl_anchors[idx_key]         = sample.stop_sequences  #list
                   
                    self.anchors[idx_key]           = {}
                    
                    #self.id_list_all[idx_key]      = []
                    self.id_list_passed[idx_key]    = []
                    self.deleted_ids[idx_key]       = {}
                    self.deleted_ids['nokey']       = {}
                    self.trimmed_ids[idx_key]       = {}
                    self.uniques[idx_key]           = {}
                    self.names[idx_key]             = {}
                    self.fa[idx_key]       = FastaOutput(os.path.join(self.trimming_dir, idx_key) + ".trimmed.fa")
                    
          
          
                    #####################
                    #
                    #  PrimerSuite Class
                    #
                    #####################
                    self.psuite[idx_key] = PrimerSuite(self.runobj, self.taxonomic_domain[idx_key],self.dna_regions[idx_key],idx_key)
                    #self.runbin['psuite'][idx_key]= PrimerSuite(self.taxonomic_domain[idx_key],self.dna_regions[idx_key])
                    
                    if(self.seqDirs[idx_key] == 'F' or self.seqDirs[idx_key] == 'B'):
                        self.proximal_primers[idx_key] = self.psuite[idx_key].primer_expanded_seq_list['F']
                        self.distal_primers[idx_key]   = self.psuite[idx_key].primer_expanded_seq_list['R']
                        if self.anchor_name[idx_key]:
                            self.anchors[idx_key] = get_anchor_list(self.runobj, self.anchor_name[idx_key], self.adtnl_anchors[idx_key])
         
                        
                    if(self.seqDirs[idx_key] == 'R' or self.seqDirs[idx_key] == 'B'):
                        self.proximal_primers[idx_key] = [revcomp(primer_seq) for primer_seq in self.psuite[idx_key].primer_expanded_seq_list['R'] ]  
                        self.distal_primers[idx_key]   = [revcomp(primer_seq) for primer_seq in self.psuite[idx_key].primer_expanded_seq_list['F'] ] 
                        if self.anchor_name[idx_key]:
                            self.anchors[idx_key] = [revcomp( anchor ) for anchor in get_anchor_list(self.runobj, self.anchor_name[idx_key], self.adtnl_anchors[idx_key]) ]
                    
                    if len(self.proximal_primers[idx_key]) == 0 and len(self.distal_primers[idx_key]) == 0:
                        logger.debug("**** Didn't find any primers that match any of the domain/regions in the lane/key sections")
            elif self.runobj.platform == 'illumina':
                # create our directories for each key
                pass

                    
                    
    def trimrun_ion_torrent(self, write_files = False):
        return ('TODO','Needs Completion','Needs Completion')
        
    def trimrun_vamps(self, write_files = False):
        
        self.stats_fp = open( os.path.join(self.outdir, self.statsFileName),"w" )
        self.number_of_raw_sequences  = 0
        self.number_of_good_sequences = 0
        success_code = ()
        self.deleted_count_for_nokey =0
        self.deleted_count_for_proximal =0
        self.deleted_count_for_distal =0
        self.deleted_count_for_n =0
        self.deleted_count_for_quality =0
        self.deleted_count_for_no_insert =0
        self.deleted_count_for_minimum_length =0
        self.deleted_count_for_unknown_lane_runkey = 0

        # input type is defined in the config file and can be sff, fasta, fasta-mbl or fastq
        # need to figure out how if to get quality data with fasta    
        # save a list of read_ids: all_ids, passed_ids
        
        for file_info in self.runobj.input_file_info.values():  # usually just one file: a list of one?    
            file_format = file_info["format"]  
            
            # the illumina fastq format needs to be parsed as fastq format then specially
            # post processed :(
            if file_format == "fastq-illumina":
                parsing_format = "fastq"
            elif file_format == "fastq-sanger":
                parsing_format = "fastq"
            elif file_format == "sff":
                parsing_format = "sff-trim"
            elif file_format == "fasta-mbl":
                parsing_format = "fasta"
            else:
                parsing_format = file_format
                
            
            
            
            file_path = os.path.join(self.indir,file_info["name"])
            print 'FILE',file_path
            logger.debug(file_info["name"]+' parsing format: '+parsing_format)
            # sff and fasta (non-mbl) get their lane info from this .ini field
            lane = file_info["lane"] 
                
                
            for record in SeqIO.parse(file_path, parsing_format):            
                self.number_of_raw_sequences += 1
                
                # get ids
                if file_format == "fasta-mbl":
                    id_line_parts = record.description.split(' ')
                    id = id_line_parts[0]
                    lane = id_line_parts[2]
                elif file_format == "fastq-sanger":                    
                    id = record.id
                elif file_format == "fastq-illumina":
                    # will need lots of other stuff here for fastq-illumina
                    id = record.id
                else:
                    id = record.id
                # should merge these with above if/else
                
                if file_format == "fasta" or file_format == "fasta-mbl":
                    q_scores = ''
                    # for vamps user upload:
                    if os.path.exists(self.indir + "/qualfile_qual_clean") and os.path.getsize(self.indir + "/qualfile_qual_clean") > 0:
                        # for vamps uploads use '_clean' file 
                        # format is on one line ( created from reg fasta in upload_file.php: 
                        #  >FRZPY5Q02G73IH	37 37 37 37 37 37 37 37 37 37
                        q_scores = subprocess.check_output("grep "+ id +" " + self.indir + "/qualfile_qual_clean", shell=True).strip().split("\t")[1].split()                    
                        q_scores = [int(q) for q in q_scores]
                        
                    else:
                        logger.debug("No usable qual file found")
                elif file_format == "fastq-sanger":
                    q_scores = record.letter_annotations["phred_quality"]  
                elif file_format == "fastq-illumina":
                    q_scores = record.letter_annotations["phred_quality"]  
                else:
                    q_scores = record.letter_annotations["phred_quality"]  
                
                ########################################################
                #
                # DO TRIM: Trim Each Sequence
                #
                ######################################
                seq = record.seq.tostring().upper()
                
                trim_data = self.do_trim(id, lane, seq, q_scores)
                ######################################
                #
                #
                ########################################################

                deleted         = trim_data.get('deleted',None)
                lane_tag        = trim_data.get('idx_key',None)
                seq             = trim_data.get('trimmed_sequence',None)        
                delete_reason   = trim_data.get('delete_reason',None)
                
                # print out a fasta file for each lane_tag                
                if(lane_tag and not deleted):
                    # should have all this data
                    exact_left      = trim_data.get('exact_left',None)   
                    exact_right     = trim_data.get('exact_right',None)        
                    primer_name     = trim_data.get('primer_name',None)
                    # for the names file of unique ids
                    first_id_for_seq = self.uniques[lane_tag].get(seq, None)
                    if first_id_for_seq == None:
                        self.uniques[lane_tag][seq] = first_id_for_seq = id
                    try:
                        self.names[lane_tag][first_id_for_seq].append(id)
                    except KeyError:
                        self.names[lane_tag][id] = [id]
                    
                    if write_files:        
                        self.fa[lane_tag].write_id(id)
                        self.fa[lane_tag].write_seq(seq) 

                    self.number_of_good_sequences += 1
                    logger.debug(record.id + " Passed")
                else:
                    logger.debug(id + " 2-Deleted:" + delete_reason)
                        
        # count up some things
        count_uniques = 0
        good_idx_keys = []
        for idx_key in self.runobj.run_keys:
            count = len(self.uniques[idx_key])
            if count > 0:
                good_idx_keys.append(idx_key)
            count_uniques = count_uniques + count               
        print 'Seqs Passed: ',self.number_of_good_sequences
        return ('SUCCESS',self.number_of_good_sequences,good_idx_keys)
        
    def filter_illumina(self, write_files = False):
        from pipeline.illumina_filtering import IlluminaFiltering
        file_list = []
        # setup illumina filtering
        iFilter = IlluminaFiltering(self.runobj)
        for file in self.runobj.files_list.split(','):
            # these are available in fastq_trimmer_by_quality.py
            #   infile=None,            outfile=None, 
            #   format='sanger',        wsize=1,        wstep=1,            trim_ends='53', 
            #   agg_action='min',       exc_count=0,    score_comp='>=',    qual_score=0,   
            #   filter_first50=False,   filter_Ns=False,filter_Nx=0,        failed_fastq=False,
            #   length=0,               trim=0,         clip=0,             keep_zero_length=False
            #
            # Enter any parameter from the list above or the default (shown) will be used
            # output file is created by the infile name
            file_list.append( iFilter.trim_by_quality(infile = file, length=75, trim=0, clip=0, filter_first50=True, filter_Ns=True, failed_fastq=True) )

        return ('SUCCESS','',file_list)
        
    def trim_illumina(self, write_files = False, file_list = []):
        pass
#        import pipeline.analyze_illumina_v6_overlaps as illumina_trimmer
#        
#        
#        #probj = PerfectReads(self.runobj)
#        #create a configDict
#        # that looks like Merens from analyze-illumina-v6-overlap:
#        #-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<------------
#        #    [general]
#        #    project_name = (project name)
#        #    researcher_email = (your e-mail address)
#        #    input_directory = (directory from which the input files will be read from)
#        #    output_directory = (directory to store output)
#        #    
#        #    
#        #    [files]
#        #    pair_1 = (pair 1 files, comma separated)
#        #    pair_2 = (pair 2 files, comma separated. must be ordered based on pair 1 files)
#        #-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<------------
#        class expando(object):pass
#        config = expando()
#        config.input_directory = os.path.join(self.runobj.output_dir,"illumina_filtered")
#        # project should come from csv file but they are per dataset using idx_keys
#        # self.runobj.samples[key].project
#        config.project_name = "myproject"
#        trim_dir = os.path.join(self.runobj.output_dir,"trim")
#        if not os.path.exists(trim_dir):
#            os.mkdir(trim_dir)
#        config.output_directory = trim_dir
#        config.pair_1 = []
#        config.pair_2 = []
#        projects = []
#        # FILE SET:
#        # v6_Amplicon_IDX4_TGACCA_L003_R1_001.fastq.gz
#        # v6_Amplicon_IDX4_TGACCA_L003_R2_001.fastq.gz
#        # file_prefix = v6_Amplicon_IDX4_TGACCA_L003
#        file_prefixes=[]
#        for key in self.idx_keys:
#            project = self.runobj.samples[key].project
#            file_prefixes.append(self.runobj.samples[key].file_prefix)
#        print file_prefixes
#        for p in sorted(file_prefixes):
#            if p.find('_R1_') > 0:
#                #print 'R1 ',os.path.join(config.input_directory,file)
#                file = p+".filtered.fastq"
#                config.pair_1.append(os.path.join(config.input_directory,file))
#            if p.find('_R2_') > 0:
#                #print 'R2 ',os.path.join(config.input_directory,file)
#                file = p+".filtered.fastq"
#                config.pair_2.append(os.path.join(config.input_directory,file))
#                
#                
#        sys.exit()
#        
#        
#        for file in sorted(file_list):
#            if file.find('_R1_') > 0:
#                #print 'R1 ',os.path.join(config.input_directory,file)
#                config.pair_1.append(os.path.join(config.input_directory,file))
#            if file.find('_R2_') > 0:
#                #print 'R2 ',os.path.join(config.input_directory,file)
#                config.pair_2.append(os.path.join(config.input_directory,file))
#        
#        illumina_trimmer.main(config) 
#        
#        return ('SUCCESS','TODO',self.idx_keys)
        
        
        
    #
    #  initialize counters for deletes and open stats file
    #    for each input file
    #      for each sequence
    #         call do_trim(seq...)  (this routine will keep stats on per runkey/lane basis)
    #         write out sequence data
    def trimrun_454(self, write_files = False):
        
        self.stats_fp = open( os.path.join(self.outdir, self.statsFileName),"w" )
        self.number_of_raw_sequences  = 0
        self.number_of_good_sequences = 0
        success_code = ()
        self.deleted_count_for_nokey =0
        self.deleted_count_for_proximal =0
        self.deleted_count_for_distal =0
        self.deleted_count_for_n =0
        self.deleted_count_for_quality =0
        self.deleted_count_for_no_insert =0
        self.deleted_count_for_minimum_length =0
        self.deleted_count_for_unknown_lane_runkey = 0

        # input type is defined in the config file and can be sff, fasta, fasta-mbl or fastq
        # need to figure out how if to get quality data with fasta    
        # save a list of read_ids: all_ids, passed_ids
        
        for file_info in self.runobj.input_file_info.values():  # usually just one file: a list of one?    
            file_format = file_info["format"]     
            # the illumina fastq format needs to be parsed as fastq format then specially
            # post processed :(
            if file_format == "fastq-illumina":
                parsing_format = "fastq"
            elif file_format == "fastq-sanger":
                parsing_format = "fastq"
            elif file_format == "sff":
                parsing_format = "sff-trim"
            elif file_format == "fasta-mbl":
                parsing_format = "fasta"
            else:
                parsing_format = file_format
            file_path = os.path.join(self.indir,file_info["name"])
            logger.debug(file_info["name"]+' '+parsing_format)
            # sff and fasta (non-mbl) get their lane info from this .ini field
            lane = file_info["lane"]   
            for record in SeqIO.parse(file_path, parsing_format):            
                self.number_of_raw_sequences += 1
                
                # get ids
                if file_format == "fasta-mbl":
                    id_line_parts = record.description.split(' ')
                    id = id_line_parts[0]
                    lane = id_line_parts[2]
                elif file_format == "fastq-sanger":                    
                    id = record.id
                elif file_format == "fastq-illumina":
                	# will need lots of other stuff here for fastq-illumina
                	id = record.id
                else:
                	id = record.id
                # should merge these with above if/else
                
                if file_format == "fasta" or file_format == "fasta-mbl":
                    q_scores = ''
                    # for vamps user upload:
                    if os.path.exists(self.indir + "/qualfile_qual_clean"):
                        # for vamps uploads use '_clean' file 
                        # format is on one line ( created from reg fasta in upload_file.php: 
                        #  >FRZPY5Q02G73IH	37 37 37 37 37 37 37 37 37 37
                        q_scores = subprocess.check_output("grep "+ id +" " + self.indir + "/qualfile_qual_clean", shell=True).strip().split("\t")[1].split()                    
                        q_scores = [int(q) for q in q_scores]
                        
                    else:
                        logger.debug("No qual file found")
                elif file_format == "fastq-sanger":
                    q_scores = record.letter_annotations["phred_quality"]  
                elif file_format == "fastq-illumina":
                    q_scores = record.letter_annotations["phred_quality"]  
                else:
                    q_scores = record.letter_annotations["phred_quality"]  
                
                
                # trim each sequence
                seq = record.seq.tostring().upper()
                trim_data = self.do_trim(id, lane, seq, q_scores)
                

                deleted         = trim_data.get('deleted',None)
                lane_tag        = trim_data.get('idx_key',None)
                seq             = trim_data.get('trimmed_sequence',None)        
                delete_reason   = trim_data.get('delete_reason',None)
                
                # print out a fasta file for each lane_tag                
                if(lane_tag and not deleted):
                    # should have all this data
                    exact_left      = trim_data.get('exact_left',None)   
                    exact_right     = trim_data.get('exact_right',None)        
                    primer_name     = trim_data.get('primer_name',None)
                    # for the names file of unique ids
                    first_id_for_seq = self.uniques[lane_tag].get(seq, None)
                    if first_id_for_seq == None:
                        self.uniques[lane_tag][seq] = first_id_for_seq = id
                    try:
                        self.names[lane_tag][first_id_for_seq].append(id)
                    except KeyError:
                        self.names[lane_tag][id] = [id]
                    
                    if write_files:        
                        self.fa[lane_tag].write_id(id)
                        self.fa[lane_tag].write_seq(seq) 

                    self.number_of_good_sequences += 1
                    logger.debug(record.id + " Passed")
                else:
                    logger.debug(id + " 2-Deleted:" + delete_reason)
                        
        # count up some things
        count_uniques = 0
        good_idx_keys = []
        for idx_key in self.runobj.run_keys:
            count = len(self.uniques[idx_key])
            if count > 0:
                good_idx_keys.append(idx_key)
            count_uniques = count_uniques + count               

        return ('SUCCESS',self.number_of_good_sequences,good_idx_keys)
 
    ###################################################################
    #
    #   do_trim - Trims each sequence
    # 0. check N's
    # 1. find runkey/lane
    # 2. determine direction
    # 3. trim proximal
    # 4. trim anchor optional
    # 5. trim distal
    # 6. check length and other things
    ################################################################### 
    def do_trim(self, read_id, lane, raw_sequence, quality_scores=None):
    
        trim_collector = {}    
        logger.debug("\nTrimming read " + read_id)
        trimmed_sequence = raw_sequence
        
        delete_reason = ''
        logger.debug(read_id  + '\t' + raw_sequence)

        # check for N's
        countNs = check_for_Ns(raw_sequence)
        
        tag      = ''
        lane_tag = ''
        # try to locate and remove the runkey
        if self.runobj.force_runkey != None:
            tag = self.runobj.force_runkey
            trimmed_sequence = raw_sequence[len(tag):]
        else:
            if self.runobj.vamps_user_upload:
                keys = [i.split('_')[1] for i in self.idx_keys]
                tag, trimmed_sequence = remove_runkey(raw_sequence, keys)
            else:
                tag, trimmed_sequence = remove_runkey(raw_sequence, self.idx_keys)

        # did we find a run key
        if( not tag ):
            delete_reason = DELETED_RUNKEY
            self.deleted_count_for_nokey += 1
            self.deleted_ids['nokey'][read_id] = delete_reason
            logger.debug("deleted: " + tag)
            return {"deleted" : True,"delete_reason" : delete_reason}
        else:            
            lane_tag = lane + '_' + tag
            # is there a run_key with this tag?  Samples are keyed by this
            if self.runobj.samples.get(lane_tag,None)==None:
                # bad lane_runkey 
                delete_reason = DELETED_LANEKEY
                self.deleted_count_for_unknown_lane_runkey += 1
                logger.debug("lane/runkey error got a runkey of: " + tag + " but could not find entry in Samples for lane_tag: " + lane_tag + " deleted: " + tag)
                return {"deleted" : True,"delete_reason" : delete_reason}
            else:
                logger.debug("lane and tag found: " + lane_tag)
                logger.debug("key trimmed seq: " + trimmed_sequence)
            
        seq_direction = self.seqDirs[lane_tag] 
            
        ###################################################################
        #
        #   5-find primers
        #     defs:
        #       proximal === closest to runkey
        #       distal   === furthest from run_key
        #       forward 'F' === sequences was read R->L  ie v3v5
        #       reverse 'R' === sequence was read L->R   ie v4v6
        ###################################################################
        exactRight = ''
        exactLeft=''
        orientation=''
        primer_name=''
        proximals_exist = True if len(self.proximal_primers[lane_tag]) else False
        if proximals_exist and (seq_direction == 'F' or seq_direction == 'B'):
            exactLeft, offset, trimmed_sequence \
                  = trim_proximal_primer( self.proximal_primers[lane_tag], trimmed_sequence)
            primer_name = ''
            if ( exactLeft):  
                orientation = '+'
                primer_name = self.psuite[lane_tag].primer_names[exactLeft] if not offset else "offset " + primer_name
                logger.debug('ExactLeft(Proximal Primer Found):' + str(exactLeft) + " offset: " + str(offset) + " orientation: " +  str(orientation) + " primer_name: " + primer_name)
                logger.debug('proximal trimmed seq: ' + trimmed_sequence)
                
        # try to trim the proximal if we have some to trim
        if proximals_exist and not exactLeft and (seq_direction == 'R' or seq_direction == 'B'):
            exactLeft, offset, trimmed_sequence \
                  = trim_proximal_primer( self.proximal_primers[lane_tag], trimmed_sequence) 
            if ( exactLeft):  
                orientation = '+'
                # find the original name of this primer....need to revert to what
                # it used to be
                primer_name = self.psuite[lane_tag].primer_names_by_reverse_complement[exactLeft] if not offset else "offset " + primer_name
                logger.debug('ExactLeft(Proximal Primer Found):' + str(exactLeft) + " offset: " + str(offset) + " orientation: " +  str(orientation) + " primer_name: " + primer_name)
                logger.debug('proximal trimmed seq: ' + trimmed_sequence)

        # if we could not find a proximal but we wanted to find them then fail the sequence                
        if not exactLeft and proximals_exist:
                logger.debug('proximal primer not found...deleted reason: proximal')                        
                delete_reason = DELETED_PROXIMAL
                self.deleted_count_for_proximal += 1
                                       
        # try to trim the anchor if we have enough sequence left at this point 
        
        if len(trimmed_sequence) < int(self.runobj.minimumLength):
            delete_reason = DELETED_MINIMUM_LENGTH
            self.deleted_count_for_minimum_length += 1
        else:
            # possible anchor trimming
            if len(self.anchors[lane_tag]) > 0:
                anchor_name = self.anchor_name[lane_tag]
                logger.debug( 'Have Anchor name: ' + anchor_name)                    
                exactRight, exactTrimmedOff, trimmed_sequence \
                            = trim_anchor(anchor_name, self.anchors[lane_tag], self.runobj.anchors[anchor_name], trimmed_sequence)
                if exactRight: 
                    anchor_found = True 
                    logger.debug( 'found exactRight-anchor: ' + exactRight)
                    logger.debug('after anchor trimming seq: ' + trimmed_sequence)
                    if len(trimmed_sequence) < int(self.runobj.minimumLength):
                        delete_reason = DELETED_MINIMUM_LENGTH
                        self.deleted_count_for_minimum_length += 1
                    
            # try to find and trim distal primer if we have them and we did not already find an anchor and clip the sequence
            if not exactRight and len(self.distal_primers[lane_tag]) > 0:                
                exactRight, exactTrimmedOff, trimmed_sequence \
                            = trim_distal_primer( self.distal_primers[lane_tag], trimmed_sequence)
                if exactRight: 
                    exact_found = True
                    logger.debug( 'found distal-primer: ' + exactRight)
                    logger.debug('after distal primer trimming seq: ' + trimmed_sequence)
            
            if not exactRight and self.runobj.require_distal == True:
                logger.debug('deleted: no distal found')                        
                delete_reason = DELETED_DISTAL
                self.deleted_count_for_distal += 1
                    
        ###################################################################
        #
        #   7-check length and other sundry things
        #
        ###################################################################
        logger.debug('length raw:' + str(len(raw_sequence)) + ' length trimmed: ' + str(len(trimmed_sequence)))
        if( delete_reason != DELETED_RUNKEY and not trimmed_sequence ):
            delete_reason = DELETED_NO_INSERT
            self.deleted_count_for_no_insert += 1
        elif( delete_reason != DELETED_RUNKEY and len(trimmed_sequence) < int(self.runobj.minimumLength)):
            delete_reason = DELETED_MINIMUM_LENGTH
            self.deleted_count_for_minimum_length += 1
        
        if ( (not delete_reason) and (countNs > C.maxN) ):     
            delete_reason = DELETED_N
            self.deleted_count_for_n += 1
            logger.debug('deleted N')
        
        # save passed ids
        if(not delete_reason):
            self.id_list_passed[lane_tag].append(read_id)
        else:
            self.deleted_ids[lane_tag][read_id] = delete_reason

        
        ###################################################################
        #
        #   8-check quality
        #
        ###################################################################
        if(quality_scores and not delete_reason):
            average_score = check_for_quality(raw_sequence, trimmed_sequence, quality_scores)
            if average_score < self.runobj.minAvgQual:
                delete_reason = DELETED_QUALITY
                self.deleted_count_for_quality += 1
                
                logger.debug('deleted quality')
                
        # reverse complementing for reverse reads
        if seq_direction == 'R':
            logger.debug("doing reverse of: " + trimmed_sequence)
            trimmed_sequence = revcomp(trimmed_sequence)
            logger.debug("reverse is      : " + trimmed_sequence)
        
        trim_collector['raw_sequence']      = raw_sequence
        trim_collector['trimmed_sequence']  = trimmed_sequence
        trim_collector['exact_left']        = exactLeft
        trim_collector['exact_right']       = exactRight
        trim_collector['deleted']           = False if (not delete_reason)  else True
        trim_collector['delete_reason']     = delete_reason
        trim_collector['primer_name']       = primer_name
        trim_collector['idx_key']          = lane_tag
        trim_collector['orientation']       = orientation
       
        return trim_collector

                
    def write_data_files(self, idx_keys): 
        
        ###################################################################
        #
        #   10-print out files of unique trimmed sequence
        #           and deleted ids
        #    for each amplification:
        #       fasta file          is this needed?
        #       names file          (just like mothur)
        #       uniques fasta file  (just like mothur)
        #       abundance fasta file -formated for usearch (-uc): sorted most abund first
        #       deleted read_ids file
        #   Also there should be one file each for rawaseq, trimseq, primers and runkeys
        #       Not divided by amplification for rapidly puting in db.
        #   These are printed in out directory: './out'+rundate
        #
        ###################################################################  
        #rawseqFileName = self.outdir + '/rawseq_file.txt'
        #trimseqFileName = self.outdir + '/trimseq_file.txt'
        #f_rawseq  = open(rawseqFileName, "w")
        #f_trimseq = open(trimseqFileName,"w")
        if self.runobj.platform == 'illumina':
            return
        

        for idx_key in self.runobj.run_keys:
            self.fa[idx_key].close()  
            base_file_name = os.path.join(self.trimming_dir,idx_key)
            uniquesFileName = base_file_name + ".unique.fa"
            abundFileName   = base_file_name + ".abund.fa"
            namesFileName   = base_file_name + ".names"
            delFileName     = base_file_name + ".deleted.txt"
            # clean out old files if they exists
            remove_file(uniquesFileName)
            remove_file(abundFileName)
            remove_file(namesFileName)
            remove_file(delFileName)
            f_names = open(namesFileName,"w") 
            #  if we order the uniques by length of self.uniques[idx_key][seq] then we have abundance file
            
            # Write abund.fa file
            # mysort returns a list of tuples: (read_id, count, seq) sorted highest to lowest freq
            try:
                sorted_uniques = mysort( self.uniques[idx_key], self.names[idx_key] )
                for item in sorted_uniques:
                    read_id = item[0]
                    count = item[1]
                    seq = item[2]
                    
                    sfastaRead = read_id + ";size="+str(count)                
                    abundfa = sfasta(sfastaRead, seq)
                    abundfa.write(abundFileName,'a')

            except:
                print "**********fail abund **************"
                success_code = ('FAIL','abund',idx_key)
               
            # Write uniques.fa file
            #print 'UNIQUES',self.uniques
            try:
                for seq in self.uniques[idx_key]:
                    read_id = self.uniques[idx_key][seq]
                    print uniquesFileName,read_id,seq
                    uniquefa = sfasta(read_id, seq)
                    uniquefa.write(uniquesFileName,'a')
                logger.debug("\nwrote uniques file " + uniquesFileName)
 
            except:
                success_code = ('FAIL','unique',idx_key)    
                
            # Write names file
            try:
                for id in self.names[idx_key]:
                    others = ','.join(self.names[idx_key][id])                
                    f_names.write(id+"\t"+others+"\n")
                f_names.close()
                logger.debug("wrote names file " + namesFileName)
            except:
                success_code = ('FAIL','names',idx_key) 
                
            # Write deleted.txt file   
            if idx_key in self.deleted_ids and self.deleted_ids[idx_key]:
                f_del   = open(delFileName,  "w") 
                reason_counts = {}
                for id in self.deleted_ids[idx_key]:
                    reason = self.deleted_ids[idx_key][id]                
                    f_del.write(id+"\t"+reason+"\n")
                    current_count = reason_counts.get(reason, 0)
                    reason_counts[reason] = current_count + 1
                # now write out some stats
                f_del.write("\nTotal Passed Reads in this lane/key: " + str(len(self.names[idx_key])) + "\n")
                if(len(self.names[idx_key]) > 0):
                    for key,value in reason_counts.items():
                        f_del.write(" " + key + ": " + str(value) + " " + str(float(value*100.0)/float(len(self.names[idx_key]))) + "% of total \n")                
                else:
                    pass
                f_del.close()
                logger.debug("wrote deleted file: "  + delFileName)
            
            
            
        # print out readids that failed the key test: one file only
        if 'nokey' in self.deleted_ids and self.deleted_ids['nokey']:
            nokeyFileName = os.path.join(self.trimming_dir,'nokey.deleted.txt')
            f_del = open(nokeyFileName,"w")
            for id in  self.deleted_ids['nokey']:                
                f_del.write(id+"\tnokey\n")
            f_del.close()
            
        if True: 
            print   
            print 'Output Directory:', './'+self.outdir
            print self.number_of_raw_sequences,  "raw sequences read"  
            pct = '%.1f' % ((float(self.number_of_good_sequences)/self.number_of_raw_sequences) *100)
            print self.number_of_good_sequences, "sequences passed" ,pct+'%'
            print "Unique Counts:"
            count_uniques = 0
            good_idx_keys = []
            for idx_key in self.runobj.run_keys:
                count = len(self.uniques[idx_key])
                if count > 0:
                    good_idx_keys.append(idx_key)
                count_uniques = count_uniques + count
                print "   ",idx_key,self.dna_regions[idx_key],count
            print "   Total Uniques:",count_uniques
 
 
        #####
        #
        #  Write to stats file for this run
        #
        self.stats_fp.write("Run: "+self.run+"\n")
        self.stats_fp.write("Unique Counts:\n")
        #stats_fp.write("Run: "+self.run)
        count_uniques = 0
        for idx_key in self.runobj.run_keys:
            count = len(self.uniques[idx_key])
            count_uniques = count_uniques + count
            self.stats_fp.write("   " + str(count)+"\t"+idx_key+"\n")
        
        self.stats_fp.write("Total Uniques: "+str(count_uniques)+"\n")
        self.stats_fp.write("\nDeleted Counts (before chimera check):\n")
        self.stats_fp.write("   deleted_count_for_nokey:\t" + str(self.deleted_count_for_nokey) + "\n")
        self.stats_fp.write("   deleted_count_for_proximal:\t" + str(self.deleted_count_for_proximal)+ "\n")
        self.stats_fp.write("   deleted_count_for_distal:\t" + str(self.deleted_count_for_distal)+ "\n")
        self.stats_fp.write("   deleted_count_for_n:\t" + str(self.deleted_count_for_n)+ "\n")
        self.stats_fp.write("   deleted_count_for_quality:\t" + str(self.deleted_count_for_quality)+ "\n")
        self.stats_fp.write("   deleted_count_for_no_insert:\t" + str(self.deleted_count_for_no_insert)+ "\n")
        self.stats_fp.write("   deleted_count_for_minimum_length:\t" + str(self.deleted_count_for_minimum_length)+ "\n")
        self.stats_fp.write("   deleted_count_for_unknown_lane_runkey:\t" + str(self.deleted_count_for_unknown_lane_runkey)+ "\n")

        self.stats_fp.write("Total Deleted: "+str(self.number_of_raw_sequences-self.number_of_good_sequences)+"\n")    
        self.stats_fp.close()
        
        success_code=''
        if not success_code:
            success_code = ('SUCCESS ' + str(good_idx_keys))
        return success_code
        
def remove_file(file):
    if os.path.exists(file):
        os.remove(file)

# eventually parse the id line for illumina data
# this is scooped from meren's code!!!
def parse_sequence_id(self, seq_id):
    id_fields = seq_id.split(':')
    parsed_fields = {}
    parsed_fields['machine_name']   = id_fields[0]
    parsed_fields['run_id']         = id_fields[1]
    parsed_fields['flowcell_id']    = id_fields[2]
    parsed_fields['lane_number']    = id_fields[3]
    parsed_fields['tile_number']    = id_fields[4]
    parsed_fields['x_coord']        = id_fields[5]
    if len(id_fields[6].split()):
        parsed_fields['y_coord']    = id_fields[6].split()[0]
        parsed_fields['pair_no']    = id_fields[6].split()[1]
    else:
        parsed_fields['y_coord']    = id_fields[6]
        parsed_fields['pair_no']    = None
    parsed_fields['quality_passed'] = id_fields[7] == 'Y'
    parsed_fields['control_bits_on'] = id_fields[8]
    parsed_fields['index_sequence']  = id_fields[9]
    return parsed_fields

    
if __name__=='__main__':
    pass