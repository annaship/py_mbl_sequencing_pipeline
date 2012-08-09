import unittest
import sys
sys.path.append("../")
#from mock import Mock 
import pipeline.db_upload as dbup

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
#        for table_name in ["test.dataset", "test.run_key", "test.run", 
#                      "test.dna_region", "test.project", "test.dataset", "test.run_info_ill", 
#                      "test.sequence_ill", "test.sequence_pdr_info_ill", "test.sequence_uniq_info_ill", "test.taxonomy"]:
#            truncate_test_db_sql = "TRUNCATE %s;" % table_name
#            cls._connection.execute_no_fetch(truncate_test_db_sql)
        table_name = "test.dataset"
#        truncate_test_db_sql = "TRUNCATE test.dataset;"
        truncate_test_db_sql = "TRUNCATE %s;" % table_name
        cls._connection.execute_no_fetch(truncate_test_db_sql)

#        msql = "SET SQL_MODE=@OLD_SQL_MODE; SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS; SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;"
#        msql = "SELECT * from dataset"
#        cls._connection.execute_no_fetch(msql) 
#        cls._content = {'182': {'platform': 'Illumina', 'run_key': 'NNNNGCTAC', 'lane': '4', 'run': '20120613', 'IDX': 'GGCTAC', 'dna_region': 'v6', 'vamps_user': 'jreveillaud', 'adaptor': '', 'barcode': '', 'seq_operator': 'JV', 'overlap': 'complete', 'dataset': 'H38', 'project': 'JCR_SPO_Bv6', 'read_length': '101', 'file_prefix': 'H38', 'primer_suite': 'Bacterial v6 Suite', 'tubelabel': 'H38', 'amp_operator': 'JR', 'insert_size': '230'}}
        
    def test_1(self):
        print "URA"
#    def test_put_run_info(self):
#        content = {'182': {'platform': 'Illumina', 'run_key': 'NNNNGCTAC', 'lane': '4', 'run': '20120613', 'IDX': 'GGCTAC', 'dna_region': 'v6', 'vamps_user': 'jreveillaud', 'adaptor': '', 'barcode': '', 'seq_operator': 'JV', 'overlap': 'complete', 'dataset': 'H38', 'project': 'JCR_SPO_Bv6', 'read_length': '101', 'file_prefix': 'H38', 'primer_suite': 'Bacterial v6 Suite', 'tubelabel': 'H38', 'amp_operator': 'JR', 'insert_size': '230'}}
#        self._connection.put_run_info(content)
#
#    def test_execute_fetch_select(self): 
#        table_name = "run_info_ill"
#        id_name = table_name + "_id"
#        sql = "select %s from %s where %s = 1" % (id_name, table_name, id_name)
#        res = self._connection.execute_fetch_select(sql)
#        self.assertEqual(int(res[0][0]), 1)
#        
#    def test_execute_no_fetch(self):
#        taxonomy = "Blah; Blah; Blah"
#        sql = """INSERT IGNORE INTO taxonomy (taxonomy) VALUES ('%s')""" % (taxonomy.rstrip())
#        res = self._connection.execute_no_fetch(sql)
#        self.assertEqual(res, 1)
        
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
execute_fetch_select(self, sql) 
execute_no_fetch(self, sql) 
__init__(self, run = None) 
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

