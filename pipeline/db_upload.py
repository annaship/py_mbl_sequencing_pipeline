import sys
import fastalib as u
import MySQLdb
from os import listdir, walk
from os.path import isfile, join
import csv

from pipeline.pipelinelogging import logger
#import logging
import constants as C

class MyConnection:
    """
    Connection to env454
    By default takes parameters from "db_conn.conf", host = "newbpcdb2"
    if different: change db_conn.conf and use my_conn = MyConnection(config_file_name, server_name)
    db_conn.conf has lines as: bpcdb2:bpcdb2:3306:my_password                                                                                                            
    """
        
    """
    TODO: change hardcoded values to args: file_name="db_conn.conf", server_name="newbpcdb2_ill", user = "ashipunova", 
    """
    def __init__(self, file_name="db_conn.conf", server_name="newbpcdb2_ill"):
        self.file_name   = file_name
        self.server_name = server_name
        self.conn        = None
        self.cursor      = None
        self.rows        = 0
        
        try:
            content = [line.strip() for line in open(self.file_name).readlines()]
    #        print dir(MySQLdb)
    #        ['BINARY', 'Binary', 'Connect', 'Connection', 'DATE', 'DATETIME', 'DBAPISet', 'DataError', 'DatabaseError', 'Date', 'DateFromTicks', 'Error', 'FIELD_TYPE', 'IntegrityError', 'InterfaceError', 'InternalError', 'MySQLError', 'NULL', 'NUMBER', 'NotSupportedError', 'OperationalError', 'ProgrammingError', 'ROWID', 'STRING', 'TIME', 'TIMESTAMP', 'Time', 'TimeFromTicks', 'Timestamp', 'TimestampFromTicks', 'Warning', '__all__', '__author__', '__builtins__', '__doc__', '__file__', '__name__', '__package__', '__path__', '__revision__', '__version__', '_mysql', 'apilevel', 'connect', 'connection', 'constants', 'debug', 'escape', 'escape_dict', 'escape_sequence', 'escape_string', 'get_client_info', 'paramstyle', 'release', 'result', 'server_end', 'server_init', 'string_literal', 'test_DBAPISet_set_equality', 'test_DBAPISet_set_equality_membership', 'test_DBAPISet_set_inequality', 'test_DBAPISet_set_inequality_membership', 'thread_safe', 'threadsafety', 'times', 'version_info']
    
            for line in content:
                fields = line.split(':')
                if fields[0] == self.server_name:
                    print "server_name = " + str(self.server_name)
                    print "=" * 40
                    # conn = MySQLdb.connect (host = str(fields[1]), port = int(fields[2]), user = "ashipunova", passwd = str(fields[3]), db = "env454")
                    # self.conn = MySQLdb.connect (host = str(fields[1]), port = int(fields[2]), user = "ashipunova", passwd = str(fields[3]), db = "env454")
                    self.conn = MySQLdb.connect (host = str(fields[1]), port = int(fields[2]), user = "ashipunova", passwd = str(fields[3]), db = str(fields[4]))
                    self.cursor = self.conn.cursor()       
        except MySQLdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            raise
        except:                       # catch everything
            print "Unexpected:"         # handle unexpected exceptions
            print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
            raise                       # re-throw caught exception   

    def execute_fetch_select(self, sql):
        if self.cursor:
            self.cursor.execute(sql)
            res = self.cursor.fetchall ()
            return res

    def execute_insert(self, sql):
        if self.cursor:
            self.cursor.execute(sql)
            self.conn.commit()
#            if (self.conn.affected_rows()):
            if (self.conn.insert_id):
                self.rows +=1
#        logger.debug("rows = "  + str(self.rows))
 

class dbUpload:
    """db upload methods"""
    Name = "dbUpload"
    """
    TODO: change hardcoded values to args: server_name="newbpcdb2_ill",
        self.sequence_table_name = "sequence_ill", 
        self.sequence_field_name = "sequence_comp"  
    TODO: run_key_id into run_info_ill
    Order:
        # insert_seq(fasta.seq)
        # insert_pdr_info()
        # gast
        # insert_taxonomy()
        # insert_sequence_uniq_info_ill()

    """
    def __init__(self, run = None):

        self.run 	= run
        self.outdir = run.output_dir
        try:
            self.basedir = run.basedir
        except:
            self.basedir = self.outdir
        self.rundate     = self.run.run_date
        self.use_cluster = 1
        self.fasta_dir   = self.run.input_dir + "/fasta/" 
        self.gast_dir    = self.run.input_dir + "/gast/"
        self.filenames   = []
        self.my_conn     = MyConnection(server_name = 'newbpcdb2_ill')  
        self.sequence_table_name = "sequence_ill" 
        self.sequence_field_name = "sequence_comp" 

#        self.refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
   
    def get_fasta_file_names(self, fasta_dir):
        for (dirpath, dirname, files) in walk(fasta_dir):
            return files
        
    def get_run_info_ill_id(self, filename_base):
        my_sql = """SELECT run_info_ill_id FROM run_info_ill WHERE file_prefix = '%s'""" % (filename_base)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])
        
    def insert_seq(self, seq):
        my_sql = """INSERT IGNORE INTO %s (%s) VALUES (COMPRESS('%s'))""" % (self.sequence_table_name, self.sequence_field_name, seq)
        self.my_conn.execute_insert(my_sql)
    
    def insert_pdr_info(self, fasta, run_info_ill_id):
        # ------- insert sequence info per run/project/dataset --------
        seq_count       = int(fasta.id.split('|')[1].split(':')[1])
        my_sql          = """SELECT sequence_ill_id FROM sequence_ill WHERE COMPRESS('%s') = sequence_comp""" % (fasta.seq)
        res             = self.my_conn.execute_fetch_select(my_sql)
        sequence_ill_id = int(res[0][0])
        
        my_sql = """INSERT IGNORE INTO sequence_pdr_info_ill (run_info_ill_id, sequence_ill_id, seq_count) 
                    VALUES (%s, %s, %s)""" % (run_info_ill_id, sequence_ill_id, seq_count)
        self.my_conn.execute_insert(my_sql)
 
    def get_gasta_result(self, filename):
        gast_file_name = self.gast_dir + filename + '.gast'
#        print "gast_file_name = %s" % gast_file_name
        with open(gast_file_name) as fd:
            gast_dict = dict([(l.split("\t")[0], l.split("\t")[1:]) for l in fd])    
        return gast_dict

    def insert_taxonomy(self, fasta, gast_dict):
        (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast_dict[fasta.id]
        my_sql = """INSERT IGNORE INTO taxonomy (taxonomy) VALUES ('%s')""" % (taxonomy.rstrip())
        self.my_conn.execute_insert(my_sql)
        
    def insert_sequence_uniq_info_ill(self, fasta, gast_dict):
        (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast_dict[fasta.id]
        my_sql = """INSERT IGNORE INTO sequence_uniq_info_ill (sequence_ill_id, taxonomy_id, gast_distance, refssu_count, rank_id, refhvr_ids) VALUES
               (
                (SELECT sequence_ill_id FROM sequence_ill WHERE COMPRESS('%s') = sequence_comp),
                (SELECT taxonomy_id     FROM taxonomy WHERE taxonomy = '%s'),
                '%s',
                '%s',
                (SELECT rank_id         FROM rank WHERE rank = '%s'),
                '%s'                
               )
               """ % (fasta.seq, taxonomy, distance, refssu_count, rank, refhvr_ids.rstrip())
        self.my_conn.execute_insert(my_sql)
    
    def put_run_info(self):
        pass
#        my_csv = readCSV(self.run)
#        content = my_csv.read_csv()
#        print "\nList of lines"
#        print content.keys()
#        # To see the list of statistics available for each line
#        for k, v in content.items():
#            print k, v['dataset'], v 
#            182 H38 {'platform': 'Illumina', 'run_key': 'NNNNGCTAC', 'lane': '4', 'run': '20120613', 'IDX': 'GGCTAC', 'dna_region': 'v6', 'vamps_user': 'jreveillaud', 'adaptor': '', 'barcode': '', 'seq_operator': 'JV', 'overlap': 'complete', 'dataset': 'H38', 'project': 'JCR_SPO_Bv6', 'read_length': '101', 'file_prefix': 'H38', 'primer_suite': 'Bacterial v6 Suite', 'tubelabel': 'H38', 'amp_operator': 'JR', 'insert_size': '230'}


            