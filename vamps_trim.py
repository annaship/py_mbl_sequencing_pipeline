#!/usr/bin/env python

#!/bioware/python/bin/python

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
import random
import csv
from time import sleep
sys.path.append( '/bioware/python/lib/python2.7/site-packages/' )
import MySQLdb
import ConMySQL
import datetime
import subprocess
from pipeline.pipelinelogging import logger
from pipeline.run import Run
from pipelineprocessor import process
import logging

class FastaReader:
    def __init__(self,file_name=None):
        self.file_name = file_name
        self.h = open(self.file_name)
        self.seq = ''
        self.id = None
        self.revcomp_seq = None
        self.base_counts = None

    def next(self): 
        def read_id():
            return self.h.readline().strip()[1:]

        def read_seq():
            ret = ''
            while True:
                line = self.h.readline()
                
                while len(line) and not len(line.strip()):
                    # found empty line(s)
                    line = self.h.readline()
                
                if not len(line):
                    # EOF
                    break
                
                if line.startswith('>'):
                    # found new defline: move back to the start
                    self.h.seek(-len(line), os.SEEK_CUR)
                    break
                    
                else:
                    ret += line.strip()
                    
            return ret
        
        self.id = read_id()
        self.seq = read_seq()
        
        
        if self.id:
            return True  
            
class FastaCleanReader:
    """
      This class reads a cleaned fasta file which was created specifically for upload
      in upload_file.php. The format looks like this:
      
        >GFADN4I02FRD2O	TGGGCGTAAAGCGCGCGCAGGTGGTTCCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGA	no info
        >GFADN4I02F6C1Y	TGGGCGTAAAACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGG	no info
        >GFADN4I02G7HL4	TGGGCGTAAAGGGAGCGCAGGCGGGAAACTAAGCGGATCTTAAAAGTGCGGGGCTCAACCCCGTGATGGGGTCCGAACTGGTTTTCTTGAGTGC	no info
        >GFADN4I02JLBKU	TGGGTTTAAAGGAGCGTAGATGGATGTTTA
        
    """
    def __init__(self,file_name=None):
        self.file_name = file_name
        self.h = open(self.file_name)
        self.seq = ''
        self.id = None
        self.revcomp_seq = None
        self.base_counts = None

    def next(self):
       
        line = self.h.readline().strip().split()
        if not len(line):
            # EOF
            return            
            
        def read_id():
            #return self.h.readline().strip()[1:]
            return line[0][1:]

        def read_seq():
            return line[1]
            
           
        
        self.id = read_id()
        self.seq = read_seq()
        
        
        if self.id:
            return True    
            
def get_primers(file_base):
    fh = open(file_base+'/primfile_prim_clean','r')
    primers = {}
    for line in fh:
        items = line.strip().split("\t")
        primers[items[0]] = {'name':items[0],'direction':items[1],'sequence':items[2]}
    fh.close()
    return primers
        
def get_metadata(file_base):  
    fh = open(file_base+'/metafile_meta_clean','r')
    metadata = {}
    for line in fh:
        items = line.strip().split("\t")
        #metadata[items[0]] = {'lanekey':'1:'+items[0],'direction':items[1],'dna_region':items[2],'project':items[3],'dataset':items[4]}
        metadata['1_'+items[0]] = {  'key':items[0],        'project':items[1],             'dataset':items[2], 
                                'dna_region':items[3],          'domain':items[4],		 		'direction':items[5],
                                'project_title':items[6],       'project_description':items[7], 'dataset_description':items[8], 
                                'env_sample_source_id':items[9]
                                }
    fh.close()
    return metadata
    

def trim_file(myobject):
    """
      Doc string
    """
    require_distal  = myobject['require_distal']
    minlength       = myobject['minlength']
    maxlength       = myobject['maxlength']
    user            = myobject['user']
    runcode         = myobject['runcode']
    site            = myobject['site']
    file_type       = myobject['file_type']
    file_base       = myobject['file_base']
    datetime        = myobject['datetime']
    use_cluster         = myobject['use_cluster']
    
    primers_obj = get_primers(file_base)
    metadata_obj = get_metadata(file_base)
    # use the files from file_base directory
    # but we get the primers and keys from the database
    # which were stored there during the loading phase
    
    
    if file_type == 'fasta' or file_type == 'fasta_clean':
        # if upload file was fasta then the script upload_file.php->file_checker
        # created a '_clean' file that is convered back into a regular fasta file here
        # for mothur to unique
        file_to_trim = file_base+'/fasta_file.fa'
        fh = open(file_to_trim,'w')
        try:
            infile = file_base+'/seqfile.fa_clean'
            f = FastaCleanReader(infile)
        except:
            infile = file_base+'/seqfile_seq_clean'
            f = FastaCleanReader(infile)
        # create new fasta here for mothur to unique,names
        while f.next():
            read_id = f.id
            seq = f.seq
            #print read_id,seq
            fh.write('>'+read_id+"\n"+seq+"\n")
        
        fh.close()
        
        # get qual file if present
        # USES: clean qual file in trim_run.py
#         if os.path.exists( file_base + "/qualfile_qual_clean"):
#             qualfile_to_trim = file_base+'/fasta_file.qual'
#             fh = open(qualfile_to_trim,'w')
#             infile = file_base + "/qualfile_qual_clean"
#             f = FastaCleanReader(infile)
#             # create new fasta quality file here for trimming
#             while f.next():
#                 read_id = f.id
#                 seq = f.seq
#                 #print read_id,seq
#                 fh.write('>'+read_id+"\n"+seq+"\n")
#             
#             fh.close()
        
        # create unique and names file (for fasta file only)
        mothur_cmd = "/bioware/mothur/mothur \"#unique.seqs(fasta="+file_to_trim+"); \" ";

        subprocess.call(mothur_cmd, shell=True)
        if not os.path.exists( file_base+"/fasta_file.unique.fa" ):    
            print "Uniques fasta file: fasta_file.unique.fa, is not created. Exiting\n";
            os.exit()
    
    
        if not os.path.exists( file_base+"/fasta_file.names" ):   
            print "Names file: fasta_file.names, is not created. Exiting\n";
            os.exit()
    
        
    elif file_type[:5] == 'fastq':     
        infile = file_base+'/seqfile.fq' 
        file_to_trim = file_base+'/seqfile.fq'
    elif file_type == 'sff':
        infile = file_base+'/seqfile.sff'
        file_to_trim = file_base+'/seqfile.sff'
    else:
        logger.debug("vamps_trim.py : Input filetype ERROR "+file_type)
                
    

    ######### Create a Run here for the uploaded data ###############
    #
    # need to create a 'Run' here and then feed it to trim_run in the py pipeline
    # A run object emulates the ini file
    # and has an output directory, general section, rundate, and a list of samples.
    # A sample has direction,dna_region,taxonomic_domain,anchor,stop_sequences
    #
    ########################
    myRunDict = {}
    for r in metadata_obj:
        myRunDict[r] = {}    
    
    #lanekeys = [metadata_obj[key]['lanekey'] for key in metadata_obj]
    
    
    myRunDict['general'] = {'run_date':datetime,            'input_dir':file_base,
                            'platform':'vamps',             'require_distal':require_distal,
                            'input_file_formats':file_type, 'user':user,
                            'input_file_lane':'1',          'vamps_user_upload':True,
                            'gast_data_source':'database',  'minimumLength':minlength,
                            'maximumLength':maxlength,      'run':runcode,
                            'use_cluster':use_cluster,      'site':site,
                            'load_vamps_database':True,     'output_dir':file_base,
                            'input_files':file_to_trim,     'files_list':[file_to_trim]
                            }
    #'input_file_names':file_to_trim,
    f_primers=''
    r_primers=''
    
    #print myRunDict
    for p in primers_obj:        
        if primers_obj[p]['direction'] == 'F':
            f_primers += primers_obj[p]['sequence']+','
        elif primers_obj[p]['direction'] == 'R':
            r_primers += primers_obj[p]['sequence']+','
    f_primers = f_primers[:-1]  
    r_primers = r_primers[:-1] 
    for r in metadata_obj:
        
        #myRunDict[metadata_obj[r]['lanekey']]['data_owner'] = user
        # r       = 1_AGTC
        
        myRunDict[r]['forward_primers'] = f_primers
        myRunDict[r]['reverse_primers'] = r_primers
        myRunDict[r]['key'] = metadata_obj[r]['key']
        myRunDict[r]['direction'] = metadata_obj[r]['direction']
        myRunDict[r]['project'] = metadata_obj[r]['project']
        myRunDict[r]['dataset'] = metadata_obj[r]['dataset']
        myRunDict[r]['dna_region'] = metadata_obj[r]['dna_region']
        myRunDict[r]['taxonomic_domain'] = metadata_obj[r]['domain']
        myRunDict[r]['project_description'] = metadata_obj[r]['project_description']
        myRunDict[r]['project_title'] = metadata_obj[r]['project_title']
        myRunDict[r]['dataset_description'] = metadata_obj[r]['dataset_description']
        myRunDict[r]['env_sample_source_id'] = metadata_obj[r]['env_sample_source_id']
        
        # turn off looking for mbl primers as well as the uploaded oned
        myRunDict[r]['use_mbl_primers'] = '0'
        
	
    
    #run = Run(myRunDict, file_base, "/xraid2-2/vampsweb/"+site)
    #for i in myRunDict:
    #    print i,myRunDict[i]
    run = Run(myRunDict,  "/xraid2-2/vampsweb/"+site)
    # output_dir is created in run so add it to dict here
    #print 'samples',run.samples
    myRunDict['output_dir'] = run.output_dir
    #print myRunDict 
    #run = Run(args.configPath, args.baseoutputdirarg, os.path.dirname(os.path.realpath(__file__)))  
    #print 'OUT dir ',run.output_dir
    # now do all the work
    # steps: trim,chimera,gast,vampsupload
    steps = 'trim'
    process(run, steps)
    
    return myRunDict
    
    
def push_to_database(myobject, runDict):
    dir = runDict['output_dir']
    conn = myobject['db_conn']
    cursor = conn.cursor()
    # here we put the data if any into trimseq
    # we do not save primer information
    for name in runDict:
        if name.startswith('1') or name.startswith('2'):
            # change ':' to '_'
            try:
            	temp = name.split('_')
            	filename = name
            	key = temp[1]
            except:
            	filename = name
            	key = name
            
            file_path = dir+'/'+filename+'.trimmed.fa'
            
            if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                project = runDict[name]['project']
                dataset = runDict[name]['dataset']
                dna_source = runDict[name]['dna_region']
                domain = runDict[name]['taxonomic_domain']
                run = myobject['runcode']
                date = myobject['datetime']
                user = myobject['user']
                f = FastaReader(file_path)
                while f.next():
                    length = str(len(f.seq))
                    trimseq_insert = "insert ignore into vamps_upload_trimseq (read_id,sequence,project,dataset,run_key, \
                                    deleted,source,length,run,entry_date,domain,user) \
                                   VALUES('"+f.id+"','"+f.seq+"','"+project+"','"+dataset+"','"+key+"','0', \
                                          '"+dna_source+"','"+length+"','"+run+"','"+date+"','"+domain+"','"+user+"')"
                                          
                    cursor.execute(trimseq_insert)
    
    # this is required for innodb tables
    conn.commit()            
    cursor.close()
    
if __name__ == '__main__':
    import argparse
    
    # DEFAULTS
    site = 'vampsdev'
    
    user = ''  
    minlength = 30
    maxlength = ""
    

    data_object = {}
    

    
    myusage = """usage: vamps_trim.py  [options]
         
         Put user created otus into the vamps database. The OTUs must be in the
         form of a matrix file and a taxonomy file.
         
         where
            
            -site
            
            
            
           
            -r,   code
                                [default: vampsdev]
            
            -u, --user       Needed for otu naming and db access.
                             Will be retrieved from .dbconf if not supplied
            -nd
            -min
            -max
            
    
    
    """
    parser = argparse.ArgumentParser(description="" ,usage=myusage)                 
    
    
                                                     
    parser.add_argument("-site",                required=True,  action="store",   dest = "site", 
                                                    help="""database hostname: vamps or vampsdev
                                                        [default: vampsdev]""")  
    parser.add_argument("-r", "--runcode",      required=True,  action="store",   dest = "runcode", 
                                                    help="")  
 
    parser.add_argument("-u", "--user",         required=True,  action="store",   dest = "user", 
                                                    help="user name")  
    
    ## optional       
    parser.add_argument("-nd", "--no_distal",       required=False,  action='store_false', dest = "require_distal", 
                                                    default=True,    help="") 
    parser.add_argument('-min',"--minlength",       required=False, action="store",   dest = "minlength", 
                                                    help = '')                                                 
    parser.add_argument("-max","--maxlength",         required=False,  action="store",   dest = "maxlength", 
                                                    help="")             
    parser.add_argument("-file_type",               required=True,  action="store",   dest = "file_type", default='fasta',
                                                    help="sff, fasta or fastq")                                                     
    parser.add_argument('-file_base',               required=True, action="store",   dest = "file_base", 
                                                    help = '') 
    parser.add_argument("-cl", "--use_cluster",     required=False,  action="store",   dest = "use_cluster", default=True,
                                                        help = '')                                                     
    logger.info("Starting vamps_trim.py")
    args = parser.parse_args()
    
    data_object['datetime'] = str(datetime.date.today())
    
    data_object['runcode'] = args.runcode
    data_object['site'] =   args.site
    data_object['user'] =   args.user
    data_object['require_distal'] = args.require_distal
    data_object['use_cluster']      = args.use_cluster
    if data_object['use_cluster'] == 'True' or data_object['use_cluster'] ==  'true':
        data_object['use_cluster'] = True
    else:
        data_object['use_cluster'] = False
    if args.minlength:
        minlength = args.minlength
    data_object['minlength'] = minlength
    
    if args.maxlength:
        maxlength = args.maxlength
    data_object['maxlength'] = maxlength
    
    
    data_object['file_type'] = args.file_type
    data_object['file_base'] = args.file_base
    
    if data_object['site'] == 'vamps':
        #db_host = 'vampsdb'
        db_host = 'bpcdb2'
        #db_name = 'vamps'
        db_name = 'vamps_user_uploads'
        db_home = '/xraid2-2/vampsweb/vamps/'
    else:
        #db_host = 'vampsdev'
        db_host = 'bpcdb2'
        #db_name = 'vamps'
        db_name = 'vampsdev_user_uploads'
        db_home = '/xraid2-2/vampsweb/vampsdev/'
    
    loggerlevel = logging.DEBUG
    
    
    logger.setLevel(loggerlevel) 
    
    obj=ConMySQL.New(db_host, db_name, db_home)
    
    data_object['db_conn'] = obj.get_conn() 
    
    runDict = trim_file(data_object)
    
    push_to_database(data_object, runDict)
    
    data_object['db_conn'].close()    
    logger.info("Finished:trim "+ args.runcode)