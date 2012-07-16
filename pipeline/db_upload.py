import sys
import fastalib as u
import MySQLdb
from os import listdir, walk
from os.path import isfile, join
import csv
from pipeline.validate import Validate

from pipeline.pipelinelogging import logger
#import logging
import constants as C

class MyConnection:
    """
    Connection to env454
    Takes parameters from ~/.my.cnf, default host = "newbpcdb2", db="illumina_reads"
    if different: use my_conn = MyConnection(host, db)
    """
    def __init__(self, host="newbpcdb2", db="illumina_reads"):
        self.conn        = None
        self.cursor      = None
        self.rows        = 0
        
        try:
            print "=" * 40
            print "host = " + str(host)
            print "=" * 40

            self.conn=MySQLdb.connect(host=host, db=db, read_default_file="~/.my.cnf")
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
    def __init__(self, run = None, cfg = None):

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
        self.my_conn     = MyConnection(host = 'newbpcdb2', db="illumina_reads")  
        self.sequence_table_name = "sequence_ill" 
        self.sequence_field_name = "sequence_comp" 
        self.my_csv      = cfg 

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
        content = self.my_csv.read_csv()
#        182 H38 {'platform': 'Illumina', 'run_key': 'NNNNGCTAC', 'lane': '4', 'run': '20120613', 'IDX': 'GGCTAC', 'dna_region': 'v6', 'vamps_user': 'jreveillaud', 'adaptor': '', 'barcode': '', 'seq_operator': 'JV', 'overlap': 'complete', 'dataset': 'H38', 'project': 'JCR_SPO_Bv6', 'read_length': '101', 'file_prefix': 'H38', 'primer_suite': 'Bacterial v6 Suite', 'tubelabel': 'H38', 'amp_operator': 'JR', 'insert_size': '230'}
        
#        --------- bulk_inserts --------- 
        
        bulk_data = ['run_key', 'run', 'dna_region']
        for datum in bulk_data:
            values = list(set([content[entry][datum] for entry in content]))
            self.insert_bulk_data(datum, values)
            
#        --------- indiv_inserts --------- 

        for k, v in content.items():
#            print k, v['dataset'], v
            self.get_contact_v_info()
            self.insert_contact()
            contact_id = self.get_contact_id(v['data_owner'])
            self.insert_project(v, contact_id)
            self.insert_dataset(v) 

            self.insert_run_info(v)
#            self.insert_primer()
#        return content[1].keys()

    def insert_bulk_data(self, key, values):
        query_tmpl = "INSERT IGNORE INTO %s (%s) VALUES (%s)"
        val_tmpl   = "'%s'"
        my_sql     = query_tmpl % (key, key, '), ('.join([val_tmpl % key for key in values]))
        self.my_conn.execute_insert(my_sql)
    
    def get_contact_v_info(self):
        """
        TODO: get info from Hilary? from vamps?
        """
        pass
    def insert_contact(self):
        pass
    def get_contact_id(self, data_owner):
        my_sql = """SELECT contact_id FROM contact WHERE vamps_name = '%s'""" % (data_owner)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])        

    def insert_project(self, content_row, contact_id):
        """
        TODO: get title, project_description, funding, env_sample_source_id
        """
        my_sql = """INSERT IGNORE INTO project (project, rev_project_name, env_sample_source_id, contact_id) VALUES
        ('%s', reverse('%s'), 0, %s)
        """ % (content_row['project'], content_row['project'], contact_id)
        self.my_conn.execute_insert(my_sql)

    def insert_dataset(self, content_row):
        """
        TODO: get dataset_description
        """        
        my_sql = """INSERT IGNORE INTO dataset (dataset, dataset_description) VALUES
        ('%s', '')
        """ % (content_row['dataset'])
        self.my_conn.execute_insert(my_sql)
    
    def get_id(self, table_name, value):
        id_name = table_name + '_id'
        my_sql  = """SELECT %s FROM %s WHERE %s = '%s'""" % (id_name, table_name, table_name, value)
        res     = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])     
    
    def insert_run_info(self, content_row):
        run_key_id      = self.get_id('run_key',      content_row['run_key'])
        run_id          = self.get_id('run',          content_row['run'])
        dataset_id      = self.get_id('dataset',      content_row['dataset'])
        project_id      = self.get_id('project',      content_row['project'])
        dna_region_id   = self.get_id('dna_region',   content_row['dna_region'])
        primer_suite_id = self.get_id('primer_suite', content_row['primer_suite'])
        if (content_row['overlap'] == 'complete'):
            overlap = 0
        
        my_sql = """INSERT IGNORE INTO run_info_ill (run_key_id, run_id, lane, dataset_id, project_id, tubelabel, barcode, 
                                                    adaptor, dna_region_id, amp_operator, seq_operator, barcode_index, overlap, insert_size, 
                                                    file_prefix, read_length, primer_suite_id) 
                                            VALUES (%s, %s, %s, %s, %s, '%s', '%s',  
                                                    '%s', %s, '%s', '%s', '%s', %s, %s, 
                                                    '%s', %s, %s)
        """ % (run_key_id, run_id, content_row['lane'], dataset_id, project_id, content_row['tubelabel'], content_row['barcode'], 
               content_row['adaptor'], dna_region_id, content_row['amp_operator'], content_row['seq_operator'], content_row['barcode_index'], overlap, content_row['insert_size'],
                                                    content_row['file_prefix'], content_row['read_length'], primer_suite_id)
        self.my_conn.execute_insert(my_sql)

    def insert_primer(self):
        pass
        


            
