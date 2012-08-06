#import sys
#import fastalib as u
#import MySQLdb
#from os import listdir, walk
#from os.path import isfile, join
import csv
import re

from pipeline.pipelinelogging import logger
#import logging
#import constants as C
            
class readCSV:
    """read CSV submition data"""
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
        self.file_path = file_path
        self.filenames = []
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
        self.end_commas                        = 0


    def read_csv(self):
        content = {}
        headers = None
        reader  = csv.reader(open(self.file_path, 'rb'), delimiter=',')
        rownum  = 0
        empty_commas_len = 0
        for row in reader:
            if reader.line_num == 1:
                """
                If we are on the first line, create the headers list from the first row
                """
                if not row[-1]:
                    row = self.empty_ends_columns(row)                    
                headers = row
            else:
                """
                Create the sub-dictionary by using the zip() function.
                """
                "TODO: remove only the same amount of empties as header had (end_commas)"
                if self.end_commas:
                    row = self.empty_ends_columns(row)
                content[rownum] = dict(zip(headers, row))
            rownum += 1

        """
        print dir(myReader)
        ['__class__', '__delattr__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__iter__', '__new__', '__reduce__',
         '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'dialect', 'line_num', 'next']
        """ 
        return content
    
    def empty_ends_columns(self, row):
        while row[-1] is '':
            row.pop()
            self.end_commas += 1
        if self.end_commas:
            return row
    
    def create_conf(self):
        pass
