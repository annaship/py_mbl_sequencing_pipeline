import unittest
import sys
import os
sys.path.append("../")
#from mock import Mock 
import pipeline.db_upload as dbup
sys.path.append("/bioware/pythonmodules/illumina-utils/")
sys.path.append("/Users/ashipunova/bin/illumina-utils")
import fastalib as u

from pipeline.run import Run
import test.run_object_factory as fake_data_object


class DbUloadTestCase(unittest.TestCase): 
    @classmethod  
    def setUpClass(cls):
        cls._connection = dbup.MyConnection(host = "vampsdev", db = "test")
        msql = "SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;" 
        cls._connection.execute_no_fetch(msql) 
        msql = "SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;" 
        cls._connection.execute_no_fetch(msql) 
        msql = "SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL';" 
        cls._connection.execute_no_fetch(msql) 
        
        data_object = fake_data_object.data_object
        pi_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
        cls._runobj = Run(data_object, pi_path)    
        cls._my_db_upload = dbup.dbUpload(cls._runobj)

        cls.filenames   = []
        cls.seq_id_dict = {}
        fasta_file_path = "./test/sample_data/illumina/result/20120614/analysis/perfect_reads/SMPL53_3-PERFECT_reads.fa.unique"
        cls.fasta           = u.SequenceSource(fasta_file_path, lazy_init = False) 


    @classmethod  
    def tearDownClass(cls):
        msql = "SET SQL_MODE=@OLD_SQL_MODE"
        cls._connection.execute_no_fetch(msql) 
        msql = "SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;"
        cls._connection.execute_no_fetch(msql) 
        msql = "SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;"
#        msql = "SET SQL_MODE=@OLD_SQL_MODE; SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS; SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;"
        cls._connection.execute_no_fetch(msql) 
        print "\nDone!"
        
#    def test_1(self):
#        print "URA"
#
    "Run setUp to clean db and fill out run info"
    @unittest.skip("Needs clean db")    
    def test_a_setUpCleanDb(self):
        for table_name in ["test.dataset", "test.run_key", "test.run", 
                      "test.dna_region", "test.project", "test.dataset", "test.run_info_ill", 
                      "test.sequence_ill", "test.sequence_pdr_info_ill", "test.sequence_uniq_info_ill", "test.taxonomy"]:
            truncate_test_db_sql = "TRUNCATE %s;" % table_name
            self._connection.execute_no_fetch(truncate_test_db_sql)
    
    @unittest.skip("Needs clean db")
    def test_b_setUpRunInfo(self):
        my_read_csv = dbup.dbUpload(self._runobj)
        my_read_csv.put_run_info()
        print "done with put_run_info" 
        
    @unittest.skip("Needs clean db")
    def test_c_execute_fetch_select(self): 
        msql = 'INSERT INTO run_info_ill VALUES ("1", "1529", "2164", "8", "6951", "2411", "83", "", "", "19", "JV", "JV", "GCCTAA", "0", "230", "6_FP1BermC_6_14_10_CGCTC", "101", "23")'
        self._connection.execute_no_fetch(msql) 
        
        table_name = "run_info_ill"
        id_name = table_name + "_id"
        sql = "select %s from %s where %s = 1" % (id_name, table_name, id_name)
        res = self._connection.execute_fetch_select(sql)
        self.assertEqual(int(res[0][0]), 1)
#        
    @unittest.skip("Needs clean db")
    def test_d_execute_no_fetch(self):
        taxonomy = "Blah; Blah; Blah"
        sql = """INSERT IGNORE INTO taxonomy (taxonomy) VALUES ('%s')""" % (taxonomy.rstrip())
        res = self._connection.execute_no_fetch(sql)
        "taxonomy not exists"
        self.assertEqual(res, 1)

    @unittest.skip("Run after the previous one")
    def test_e_taxonomy_exists(self):
        taxonomy = "Blah; Blah; Blah"
        sql = """INSERT IGNORE INTO taxonomy (taxonomy) VALUES ('%s')""" % (taxonomy.rstrip())
        res = self._connection.execute_no_fetch(sql)
        "taxonomy exists, nothing inserted"
        self.assertEqual(res, 0)
    
        "FIrst do: illumina_files time = 136.972903013"
    
    def test_f_get_fasta_file_names(self):
        filenames = self._my_db_upload.get_fasta_file_names()
        file_names_list = fake_data_object.file_names_list
        self.assertEqual(filenames, file_names_list)
    
    def test_g_get_run_info_ill_id(self):
        filename_base   = "SMPL31_3"
        run_info_ill_id = self._my_db_upload.get_run_info_ill_id(filename_base)        
        
        sql = "SELECT run_info_ill_id FROM run_info_ill WHERE file_prefix = '%s'" % (filename_base)
        res = self._connection.execute_fetch_select(sql)
        self.assertEqual(run_info_ill_id, int(res[0][0]))
    
    def test_h_insert_seq(self, sequences = ['TACCCTTGACATCATCAGAACTTGTCAGAGATGACTCGGTGCCTTCGGGAACTGATAGAC']):
        self._my_db_upload.insert_seq(sequences)
        
        for seq in sequences:
            sql = "SELECT sequence_ill_id FROM sequence_ill WHERE uncompress(sequence_comp) = '%s'" % (seq)
            res = self._connection.execute_fetch_select(sql)
            self.assertEqual(int(res[0][0]), 1)
            break
    
    def test_i_insert_pdr_info(self):
        sql = "truncate sequence_pdr_info_ill"
        self._connection.execute_no_fetch(sql)
        self._connection.execute_fetch_select(sql)        
        
        self._my_db_upload.seq_id_dict = {'TACCCTTGACATCATCAGAACTTGTCAGAGATGACTCGGTGCCTTCGGGAACTGATAGAC': 1}
        
        sql = "SELECT run_info_ill_id FROM run_info_ill WHERE file_prefix = 'SMPL31_3'"
        res = self._connection.execute_fetch_select(sql)        
        run_info_ill_id = int(res[0][0])
        
        self.fasta.seq  = "TACCCTTGACATCATCAGAACTTGTCAGAGATGACTCGGTGCCTTCGGGAACTGATAGAC"
        self.fasta.id   = "A5BCDEF3:25:Z987YXWUQ:3:1101:4387:2211 1:N:0:ATCACG|frequency:1"

        res_id = self._my_db_upload.insert_pdr_info(self.fasta, run_info_ill_id)
        self.assertEqual(res_id, 1)

#    def test_get_gasta_result(self, filename):

"        inset test if taxonomy exists"
"        inset test if taxonomy not exists"
"""
insert dataset
insert run_key
insert run
insert dna_region
insert project
insert dataset
insert run_info_ill
insert sequence_ill
insert sequence_pdr_info_ill
insert sequence_uniq_info_ill
insert taxonomy

methods:
__init__(self, run = None) 
    execute_fetch_select(self, sql) 
    execute_no_fetch(self, sql) 
    get_fasta_file_names(self) 
    get_run_info_ill_id(self, filename_base) 
    insert_seq(self, sequences) 
get_seq_id_dict(self, sequences) 
get_id(self, table_name, value) 
get_sequence_id(self, seq) 
insert_pdr_info(self, fasta, run_info_ill_id) 
get_gasta_result(self, filename) 
insert_taxonomy(self, fasta, gast_dict) 
insert_sequence_uniq_info_ill(self, fasta, gast_dict) 
put_run_info(self, content = None) 
insert_bulk_data(self, key, values) 
get_contact_v_info(self) 
insert_contact(self) 
get_contact_id(self, data_owner) 
insert_rundate(self) 
insert_project(self, content_row, contact_id) 
insert_dataset(self, content_row) 
insert_run_info(self, content_row) 
insert_primer(self) 
del_sequence_uniq_info(self) 
del_sequences(self) 
del_sequence_pdr_info(self) 
del_run_info(self) 
count_sequence_pdr_info_ill(self) 
count_seq_from_file(self) 
check_seq_upload(self) 
put_seq_statistics_in_file(self, filename, seq_in_file)
"""

if __name__ == '__main__':
    unittest.main()

