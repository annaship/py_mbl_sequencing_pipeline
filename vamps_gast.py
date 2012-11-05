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
import time
import random
import csv
from time import sleep

sys.path.append( '/bioware/python/lib/python2.7/site-packages/' )
sys.path.append("/xraid/bioware/linux/seqinfo/bin")

import datetime
import subprocess
from pipeline.pipelinelogging import logger
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
            

    

def start_gast(myobject):
    """
      Doc string
    """
    project     = myobject['project']
    dataset     = myobject['dataset']
    dna_region  = myobject['dna_region']
    domain      = myobject['domain']
    runcode     = myobject['runcode']
    site        = myobject['site']
    #user_cursor = myobject['user_cursor']
    datetime    = myobject['datetime']
    user        = myobject['user']
    from_fasta  = myobject['from_fasta']
    load_db     = myobject['load_db']
    env_source_id   = myobject['env_source_id'] 
    steps   = myobject['steps'] 
    fasta_file_from_cl  = myobject['fasta_file']
    use_cluster         = myobject['use_cluster']
    #myobject['baseoutputdir']
    seq_count   = 0
    site_base   = '/xraid2-2/vampsweb/'+site
    file_prefix = user+runcode
    
    output_dir = myobject['output_dir']
    #output_dir  = os.path.join(site_base, 'tmp',user+"_"+runcode+'_gast')
    
    # use the files from file_base directory
    # but we get the primers and keys from the database
    # which were stored there during the loading phase      

    # check for directory:  user_runcode
    #    if present use the data from there
    #    if not: go to the database
    if os.path.exists(output_dir):
        print "files path exists:",output_dir
        #gast_input_source = 'files'
        #file_base = output_dir
        # This may be a mobedac upload and we should try to use the files here
        # rather than look to the database for data

    else:
        output_dir  = os.path.join(site_base, 'tmp',user+"_"+runcode+'_gast')
        print "Files path doesn't exist: attempting to get data from database"
        print "Creating directory",output_dir
        os.mkdir(output_dir)
        
    from pipeline.run import Run
    from pipelineprocessor import process
    myRunDict = {}
    # this is a minimal run dictionary for the general stanza
    myRunDict['general'] = {'run_date':datetime,                'vamps_user_upload':True, 
                            'gast_input_source':'database',     'input_file_names':'vamps_upload', 
                            'input_file_lanes':'1',             'input_file_formats':'fasta',
                            'run':runcode,                      'use_cluster':use_cluster, 
                            'platform':'vamps', 
                            'user':user,                        'site':site,
                            'load_vamps_database':True,
                            'input_files':None,                 'files_list':[],
                            'output_dir':output_dir,            'file_prefix':file_prefix}
    #print myRunDict
    #
    #
    #
    run = Run(myRunDict, "/xraid2-2/vampsweb/"+site)
    #
    #
    #
    # pack the things we'll need for GAST
    run.project = project
    run.dataset = dataset
    run.load_db = load_db
    run.env_source_id=env_source_id
    run.site = site
    run.from_fasta = from_fasta
    run.fasta_file_from_cl=fasta_file_from_cl
    run.runcode = runcode
    run.user = user
    run.samples = {}
    run.dna_region = dna_region
    
    #run.basedir = file_base
#    fastaunique_cmd = '/bioware/bin/fastaunique'
    fastaunique_cmd = 'fastaunique'
    if run.from_fasta:
        print run.from_fasta
        # copy file to
        fasta_file = os.path.join(output_dir,run.user+run.runcode+'.fa')
        shutil.copyfile(run.fasta_file_from_cl, fasta_file)
        grep_cmd = ['grep','-c','>',fasta_file]
        run.dataset_count = subprocess.check_output(grep_cmd).strip()
        
    else:
        # from database
        from pipeline.db_upload import MyConnection
        if site == 'vamps':
            db_host_user    = 'bpcdb2'
            db_name_user    = 'vamps_user_uploads'
        else:
            db_host_user    = 'bpcdb2'
            db_name_user    = 'vampsdev_user_uploads'
        myconn = MyConnection(host=db_host_user,db=db_name_user)
        
        
        # should create the fasta file and names file here and not in gast.py 
        ds_list = []
        if dataset:
            ds_list = [dataset]
            query ="select read_id,sequence,dataset from vamps_upload_trimseq where project='"+project+"' and dataset='"+dataset+"' and user='"+user+"' "    
            print query
            rows = myconn.execute_fetch_select(query)
            fasta_file = os.path.join(output_dir, 'fasta.fa')
            unique_file = os.path.join(output_dir, 'unique.fa')
            names_file = os.path.join(output_dir, 'names')
    
            fh = open(fasta_file, 'w')
            
            if not rows:
                print "No data found using query:", query
                
            for r in rows:
                id  = r[0]
                seq = r[1]                   
                fh.write(">"+id+"\n"+seq+"\n")
                
            fh.close()
            fastaunique_cmd = fastaunique_cmd +" -x -i "+fasta_file+" -o "+unique_file+" -n "+names_file 
            subprocess.call(fastaunique_cmd, shell=True)
        
        else:
            # looks for vamps_projects_datasets_pipe in vamps_user_uploads
           
            q0 = "select distinct dataset from vamps_projects_datasets_pipe where project='"+project+"' and dataset != '' and dataset != 'NoKey'"
            print q0
            dsrows = myconn.execute_fetch_select(q0)
            if not dsrows:
                print "No datasets found using query:", q0
                sys.exit()
            for ds in dsrows:                
                ds  = ds[0]
                ds_list.append(ds)
                
                
                query ="select read_id, sequence, dataset from vamps_upload_trimseq where project='"+project+"' and dataset='"+ds+"' and user='"+user+"' "    
                print query
                rows = myconn.execute_fetch_select(query)
                
                
                ds_dir = os.path.join(output_dir, ds)
                if os.path.exists(ds_dir):
                    # Start with and empty directory
                    shutil.rmtree(ds_dir, True)
                    os.mkdir(ds_dir)
                else:
                    os.mkdir(ds_dir)
                fasta_file = os.path.join(output_dir, ds,  'fasta.fa')
                unique_file = os.path.join(output_dir, ds, 'unique.fa')
                names_file = os.path.join(output_dir,  ds, 'names')
                #dataset_file=os.path.join(output_dir,   'datasets')
                fh = open(fasta_file, 'w')
                
                if not rows:
                    print "No data found using query:", query
                    
                for r in rows:
                    id  = r[0]
                    seq = r[1]  
                    ds  = r[2]
                    fh.write(">"+id+"\n"+seq+"\n")
                    
                fh.close()
                
            
                fastaunique_call = fastaunique_cmd +" "+fasta_file+" -o "+unique_file+" -n "+names_file + " -f"
            
                subprocess.call(fastaunique_call, shell=True)
    run.datasets = ds_list
    
    
    ###############################################################
    # This starts the MBL GAST python pipeline at the GAST STEP
    #
    # now do all the work
    # possible steps: trim,chimera,gast,vampsupload
    
    process(run, steps)
    print "done with gast"
    # is this sequential - YES??
    
   

        
if __name__ == '__main__':
    import argparse
    
    # DEFAULTS
    site = 'vampsdev'
    
    user = ''  
    minlength = 30
    maxlength = ""
    

    data_object = {}
    

    
    myusage = """usage: vamps_gast.py  [options]
         
         
         
         where
            
            -site       vamps or vampsdev        
           
            -r,   code
            
            -u, --user       Needed for otu naming and db access.
                             
            -project
            -dataset
            -reg, --dna_region
            -dom, --domain
            --from_fasta  T/F
            --fasta_file
    
    
    """
    parser = argparse.ArgumentParser(description="" ,usage=myusage)                 
    
    
                                                     
    parser.add_argument("-site",                    required=True,  action="store",   dest = "site", 
                                                        help="""database hostname: vamps or vampsdev
                                                        [default: vampsdev]""")  
    parser.add_argument("-r", "--runcode",          required=True,  action="store",   dest = "runcode", 
                                                        help="")  
    parser.add_argument("-u", "--user",             required=True,  action="store",   dest = "user", 
                                                        help="user name")         
    parser.add_argument("-p", "--project",          required=True,  action='store', dest = "project", 
                                                        help="") 
    parser.add_argument('-d',"--dataset",           required=False,  action="store",   dest = "dataset", default='',
                                                        help = '')                                                 
    parser.add_argument('-reg',"--dna_region",       required=True, action="store",   dest = "dna_region", 
                                                        help = '') 
    parser.add_argument('-dom',"--domain",          required=True,  action="store",   dest = "domain", 
                                                        help = '') 
    parser.add_argument("--from_fasta",             required=False,  action="store_true",   dest = "from_fasta", default=False,
                                                        help = '')                                                 
    parser.add_argument("--fasta_file",             required=False,  action="store",   dest = "fasta_file", default='',
                                                        help = '')  
    parser.add_argument("-b", "--baseoutputdir",    required=False,  action="store",   dest = "baseoutputdir", default='/xraid2-2/vampsweb/vampsdev/tmp',
                                                        help = '')
    parser.add_argument("-out", "--output_dir",    required=False,  action="store",   dest = "output_dir", default='',
                                                        help = '')                                                          
    parser.add_argument("-load", "--load_database", required=False,  action="store_true",   dest = "load_db", default=True,
                                                        help = '')                                              
    parser.add_argument("-env", "--envsource",      required=False,  action="store",   dest = "env_source_id", default='100',
                                                        help = '')
    parser.add_argument("-cl", "--use_cluster",     required=False,  action="store",   dest = "use_cluster", default=True,
                                                        help = '')                                                
    
    parser.add_argument('-l', '--loglevel',         required=False,   action="store",          dest = "loglevel",          default='ERROR',       
                                                        help = 'Sets logging level... DEBUG, INFO, WARNING, ERROR, CRITICAL')                                             
                                                
                                                
    steps ='gast,vampsupload'
    #steps ='gast'
    #steps ='vampsupload'
    args = parser.parse_args()
    
    # set logging
    loggerlevel = args.loglevel.upper()
    print "\nLog Level set to:",loggerlevel    
    logger.setLevel(loggerlevel)
    
    logger.info("Starting vamps_gast.py")
    
    # fill command line object
    data_object['datetime'] = str(datetime.date.today())
    data_object['runcode']    = args.runcode
    data_object['site']       = args.site
    data_object['user']       = args.user
    data_object['project']    = args.project[:1].capitalize() + args.project[1:]
    data_object['dataset']    = args.dataset
    data_object['dna_region'] = args.dna_region
    data_object['domain']     = args.domain
    data_object['from_fasta']     = args.from_fasta
    data_object['fasta_file']     = args.fasta_file
    data_object['baseoutputdir'] = args.baseoutputdir
    data_object['output_dir'] = args.output_dir
    data_object['load_db']      = args.load_db
    data_object['env_source_id']      = args.env_source_id
    data_object['steps']      = steps
    
    data_object['use_cluster']      = args.use_cluster
    if data_object['use_cluster'] == 'True' or data_object['use_cluster'] ==  'true':
        data_object['use_cluster'] = True
    else:
        data_object['use_cluster'] = False
#     if data_object['site'] == 'vamps':
#         db_host_vamps   = 'vampsdb'
#         db_host_user    = 'bpcdb2'
#         db_name_vamps   = 'vamps'
#         db_name_user    = 'vamps_user_uploads'
#         db_home         = '/xraid2-2/vampsweb/vamps/'
#     else:
#         db_host_vamps   = 'vampsdev'
#         db_host_user    = 'bpcdb2'
#         db_name_vamps   = 'vamps'
#         db_name_user    = 'vampsdev_user_uploads'
#         db_home         = '/xraid2-2/vampsweb/vampsdev/'
#     
#     
#     vamps_obj=ConMySQL.New(db_host_vamps, db_name_vamps, db_home)
#     data_object['vamps_cursor'] = vamps_obj.get_cursor()
#     user_obj=ConMySQL.New(db_host_user, db_name_user, db_home)
#     data_object['user_cursor'] = user_obj.get_cursor()
    

    
    
    start_gast(data_object)
        
    logger.info("Finished:gast "+args.runcode)