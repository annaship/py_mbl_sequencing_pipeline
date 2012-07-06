#import sys
#import fastalib as u
#import MySQLdb
#from os import listdir, walk
#from os.path import isfile, join
import csv

from pipeline.pipelinelogging import logger
#import logging
#import constants as C
            
class readCSV:
    """read a CSV submition data"""
    Name = "readCSV"
    """
    TODO: run_key_id into run_info_ill
    write into ini:
    for item in thelist:
      thefile.write("%s\n" % item)
    """
    def __init__(self, run = None, file_path = None):
        if run:
            self.run         = run
            self.outdir      = run.output_dir
            try:
                self.basedir = run.basedir
            except:
                self.basedir = self.outdir
            self.rundate     = self.run.run_date
            self.use_cluster = 1
            self.csv_dir   = self.run.input_dir + "/csv/" 
        self.filenames   = []
        #Tables:
        self.dataset_table_name                = "dataset"
        self.dna_region_table_name             = "dna_region"
        self.primer_suite_table_name           = "primer_suite"
        self.project_table_name                = "project"
        self.rank_table_name                   = "rank"
        self.run_table_name                    = "run"
        self.run_info_ill_table_name           = "run_info_ill"
        self.run_key_table_name                = "run_key"
        self.sequence_ill_table_name           = "sequence_ill"
        self.sequence_pdr_info_ill_table_name  = "sequence_pdr_info_ill"
        self.sequence_uniq_info_ill_table_name = "sequence_uniq_info_ill"
        self.taxonomy_table_name               = "taxonomy"
        #Fields:
        self.sequence_field_name               = "sequence_comp" 


    def read_csv(self):
        content = {}
        headers = None

        inputFile = '/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/csv/bpc_metadata_JCR_SPO_Bv6_1.csv'
#            spamReader = csv.reader(open('/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/csv/bpc_metadata_JCR_SPO_Bv6_1.csv', 'rb'), delimiter=' ', quotechar='|')
        reader = csv.reader(open(inputFile, 'rb'), delimiter=',')
        rownum = 0
        for row in reader:
            if reader.line_num == 1:
                """
                If we are on the first line, create the headers list from the first row
                """
                headers = row
            else:
                """
                Create the sub-dictionary by using the zip() function.
                """
                content[rownum] = dict(zip(headers, row))
#                content[row[0]]['Stabling'] = [s.strip() for s in content[row[0]]['Stabling'].split(',')]
            rownum += 1

#        
        """
        print dir(myReader)
        ['__class__', '__delattr__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__iter__', '__new__', '__reduce__',
         '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'dialect', 'line_num', 'next']
        """ 
        return content
    
    def create_conf(self):
        pass
