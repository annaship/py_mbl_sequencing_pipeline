#!/usr/bin/env python

##!/usr/local/www/vamps/software/python/bin/python
##!/usr/bin/env python
##!/usr/local/epd_python/bin/python
##!/bioware/python/bin/python
##!/usr/bin/env python
##!/usr/local/www/vamps/software/python/bin/python
###!/usr/bin/env python

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
from pipeline.pipelinelogging import logger

# from pipeline.utils import *
# from pipeline.sample import Sample
# from pipeline.runconfig import RunConfig
# from pipeline.run import Run
# from pipeline.trim_run import TrimRun
# from pipeline.chimera import Chimera
# from pipeline.gast import Gast
# from pipeline.vamps import Vamps
# from pipeline.fasta_mbl_pipeline import MBLPipelineFastaUtils
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
            
#import pipeline.constants as C
def load_raw_seqs(myobject):
    """
        Function to load raw sequences into rawseq table before trimming.
        Called from upload_file.php class
        
    """
    # in file may be plain text fasta, plain text fastq, or sff 
    # if file type is fasta it will have been cleaned in file_checker
    infile = myobject['infile']
    type = myobject['type']
    cursor  = myobject['cursor']
    user    = myobject['user']
    runcode = myobject['runcode']
    date    = myobject['datetime']
    file_type = myobject['file_type']
    
    if file_type == 'fasta':
        parse_type = 'fasta'
        f = FastaCleanReader(infile)
        # create new fasta here for mothur to unique,names
        while f.next():
            read_id = f.id
            seq = f.seq
            #print read_id,seq
            
            
    else:
    
        if file_type == 'sff':
            parse_type = 'sff-trim'
        else:
            parse_type = 'fastq'
            # this will be a clean file in 
        
        from Bio import SeqIO
        
        for record in SeqIO.parse(infile, parse_type):                
            
            id = record.id    
            seq = str(record.seq)
            length = str(len(seq))
            cursor.execute("insert ignore into vamps_upload_rawseq \
                            (read_id,sequence,length,user,run,entry_date) \
                            VALUES('"+id+"', '"+seq+"', '"+length+"', '"+user+"', '"+runcode+"', '"+date+"') ")
                                
            # file can only be fasta file
    
    
def load_trimmed_seqs(myobject):
    """
        Function to load trimmed sequences into trimseq table.
        Input file is a '_clean' file: similar to a fasta but on one line:
        >30     pre-trimmed     0               unknown -	-	-	0	0	AAGGCTTGACATATAGAGGAAACGTCTGGAA
        >21     pre-trimmed     0               unknown -	-	-	0	0	TGGGCTCGAATGGTATTGGACTGCTGGTGAA
        Called from upload_file.php class
    """
    infile  = myobject['infile']
    # should be seqfile.fa_clean
    project = myobject['project']
    dataset = myobject['dataset']
    dna_region  = myobject['dna_region']
    runcode = myobject['runcode']
    user    = myobject['user']
    date    = myobject['datetime']
    cursor  = myobject['cursor']
    domain  = myobject['domain']
    if not domain:
        if dna_region == 'v9':
            domain = 'eukarya'
        elif dna_region[-1:] == 'a':
            domain = 'archaea'
            dna_region = dna_region[:-1]
        else:
            domain = 'bacteria'
    
    
    insert_table = 'vamps_upload_trimseq';
    loadquery =  "load data local infile '"+infile+" '\n"
    loadquery += " replace into table "+insert_table+" \n"
    loadquery += " fields terminated by '\t' "
    loadquery += " lines starting by '>'\n"
    loadquery += " set entry_date = '"+date+"', \n"
    loadquery += " source = '"+dna_region+"', \n"
    loadquery += " length = char_length(sequence), \n"
    loadquery += " project = '"+project+"', \n"
    loadquery += " dataset = '"+dataset+"', \n"
    loadquery += " run = '"+runcode+"', \n"
    loadquery += " user = '"+user+"', \n"
    loadquery += " domain = '"+domain+"', \n"
    loadquery += " countN = '0', \n"
    loadquery += " delete_reason = ''"
    print loadquery
    cursor.execute(loadquery)
    
    
        
if __name__ == '__main__':
    import argparse
    
    # DEFAULTS
    site = 'vampsdev'
    user = ''  
    
    
    
    data_object = {}
    file_type = "fasta"
    seq_file = ""
    sep = "--"
    dna_region = 'v6'
    project = ""
    dataset = ""
    
    myusage = """usage: otu_matrix2db.py -im matrixFile -it taxFile [options]
         
         Put user created otus into the vamps database. The OTUs must be in the
         form of a matrix file and a taxonomy file.
         
         where
            -i, --infile The name of the input file file.  [required]
            
            -t, --type    raw or trim.   [required]
            
            -p, --project    The name of the project.   [optional]
            
            -d, --dataset     The name of the dataset. 
                                [default: unknown]
            -dna_region       
                                   
            --site            vamps or vampsdev.
                                [default: vampsdev]
            -r, -runcode            runcode
                                [default: random number]
            -u, --user       
            
    
    
    """
    parser = argparse.ArgumentParser(description="" ,usage=myusage)
    parser.add_argument('-i', '--infile',       required=True, action="store",   dest = "infile", 
                                                    help = '')                                                 
    parser.add_argument("-t", "--upload_type",         required=True,  action="store",   dest = "type", 
                                                    help="raw or trimmed")                                   
    
    
                                                     
    parser.add_argument("-site",                required=True,  action="store",   dest = "site", 
                                                    help="""""")  
    parser.add_argument("-r", "--runcode",      required=True,  action="store",   dest = "runcode", 
                                                    help="""""")  
 
    parser.add_argument("-u", "--user",         required=True,  action="store",   dest = "user", 
                                                    help="user name")  
    parser.add_argument("-file_type",           required=True,  action="store",   dest = "file_type", 
                                                    help="sff, fasta or fastq") 
    parser.add_argument('-file_base',               required=True, action="store",   dest = "file_base", 
                                                    help = 'where the files are loacated')                                                   
                                                    
## optional       
    parser.add_argument("-dna_region",       required=False,  action="store",   dest = "dna_region", default='unknown',
                                                    help="") 
    parser.add_argument("-domain",       required=False,  action="store",   dest = "domain", default='unknown',
                                                    help="")                                                 
    parser.add_argument('-d', '--dataset',      required=False, action="store",   dest = "dataset",  
                                                    help = '')                                                 
    parser.add_argument("-p", "--project",      required=False,  action="store",   dest = "project", 
                                                    help="")                                                   
                                                              
                                                        
    
    logger.info("Starting vamps_load.py")
    args = parser.parse_args()
    
    data_object['infile'] = args.infile
    data_object['datetime'] = str(datetime.date.today())
    data_object['type'] = args.type
    data_object['runcode'] = args.runcode
    data_object['site'] =  args.site
    data_object['user'] =  args.user
    data_object['file_base'] =  args.file_base
    data_object['file_type'] =  args.file_type
    
    if args.dna_region:
        dna_region = args.dna_region
    data_object['dna_region'] = dna_region
    
    if args.domain:
        domain = args.domain
    data_object['domain'] = domain
    
    if args.project:
        project = args.project
    data_object['project'] = project
    
    if args.dataset:
        dataset = args.dataset
    data_object['dataset'] = dataset
    
    
    
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
    
    
    obj=ConMySQL.New(db_host, db_name, db_home)
    data_object['cursor'] = obj.get_cursor()
    
    if data_object['type'] == 'trimmed':        
        load_trimmed_seqs(data_object)
        
    elif data_object['type'] == 'raw':        
        load_raw_seqs(data_object)
        
    data_object['cursor'].close()    
    logger.info("Finished:load "+args.runcode ) 
    