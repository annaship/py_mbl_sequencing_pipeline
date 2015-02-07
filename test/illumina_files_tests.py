import unittest
import sys
import os
import shutil
sys.path.append("../")

import pipeline.illumina_files as ill_f

from pipeline.run import Run
# local: import test.test_factory as fake_data_object
import test_factory as fake_data_object

"""
to run: python pipeline/test/illumina_files_tests.py -v
"""

class IlluminaFilesTestCase(unittest.TestCase): 
    @classmethod  
    def setUpClass(cls):
        data_object   = fake_data_object.data_object
        #local: root_dir      = '/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test'
        root_dir = '/workspace/ashipunova/illumina/pipeline/py_mbl_sequencing_pipeline/test'
	cls.file_path = os.path.join(root_dir, data_object['general']['platform'], data_object['general']['run'], 'lane_1/analysis') 
        if os.path.isdir(cls.file_path):
            shutil.rmtree(cls.file_path)         
        pi_path             = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
        cls._runobj         = Run(data_object, pi_path)    
        cls._illumina_files = ill_f.IlluminaFiles(cls._runobj)
        cls._illumina_files.open_dataset_files()        
        cls.dataset_emails  = fake_data_object.dataset_emails

    @classmethod  
    def tearDownClass(cls):
        print "\nDone!"
        
    "Run setUp to clean db and fill out run info"
        
    def test_01_get_all_files(self):
        res = self._illumina_files.get_all_files()
        self.assertEqual(res[os.path.join(self.file_path, 'unknown.fastq')], (os.path.join(self.file_path, 'unknown'), '.fastq'))

    def test_02_create_inis(self):
        self._illumina_files.dataset_emails = self.dataset_emails
        self._illumina_files.create_inis()
        ini_files = len([f for f in os.listdir(self.file_path) if f.endswith('.ini') and os.path.isfile(os.path.join(self.file_path, f))])
        self.assertEqual(ini_files, 10)
        
        with open(os.path.join(self.file_path, "CGATGT_NNNNTCAGC_3.ini")) as file:
            for line in file.readlines():
                if line.strip() == "project_name = CGATGT_NNNNTCAGC_3":
                    print "\nValid ini\n"

    def test_03_perfect_reads(self):
        self._illumina_files.dataset_emails = self.dataset_emails
        self._illumina_files.perfect_reads()
        f_path = os.path.join(self.file_path, "reads_overlap")
        files_amount  = len([name for name in os.listdir(f_path) if os.path.isfile(os.path.join(f_path, name))])
        self.assertEqual(files_amount, 70)

    def test_04_partial_overlap(self):
        self._illumina_files.dataset_emails = self.dataset_emails
        self._illumina_files.partial_overlap_reads()
        f_path = os.path.join(self.file_path, "reads_overlap")
        files_amount  = len([name for name in os.listdir(f_path) if os.path.isfile(os.path.join(f_path, name))])
        self.assertEqual(files_amount, 110)

    def test_05a_filter_mismatches(self):
        self._illumina_files.filter_mismatches()
        f_path = os.path.join(self.file_path, "reads_overlap")
        files_amount = len([name for name in os.listdir(f_path) if os.path.isfile(os.path.join(f_path, name))])
        merged_files = len([f for f in os.listdir(f_path) if f.endswith('-MERGED_FILTERED') and os.path.isfile(os.path.join(f_path, f))])       
        self.assertEqual(files_amount, 150)
        self.assertEqual(merged_files, 10)
        
    def test_05_uniq_fa(self):
        self._illumina_files.uniq_fa()
        f_path = os.path.join(self.file_path, "reads_overlap")
        files_amount  = len([name for name in os.listdir(f_path) if os.path.isfile(os.path.join(f_path, name))])
        uniq_files = len([f for f in os.listdir(f_path) if f.endswith('-PERFECT_reads.fa.unique') and os.path.isfile(os.path.join(f_path, f))])       
        self.assertEqual(files_amount, 150)
        self.assertEqual(uniq_files, 10)
        
    def test_06_get_fastq_file_names(self):
        self.maxDiff = None
        (in_files_r1, in_files_r2) = self._illumina_files.get_fastq_file_names(self._runobj.input_dir)

#        in_files_r1_compare = './results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX1_ATCACG_L003_R1_001.filtered.fastq' in in_files_r1
        r1file_name = os.path.join(self._runobj.input_dir, 'Sample_v6_Amplicon_IDX1/IDX1_ATCACG_L003_R1_001.fastq')
#        in_files_r1_compare = "./test/sample_data/illumina/Project_J_v6_30/Sample_v6_Amplicon_IDX1/v6_Amplicon_IDX1_ATCACG_L003_R1_001.fastq" in in_files_r1
#        in_files_r2_compare = ('./test/sample_data/illumina/Project_J_v6_30/Sample_v6_Amplicon_IDX1/v6_Amplicon_IDX1_ATCACG_L003_R2_001.fastq.gz' in in_files_r2) or ('./test/sample_data/illumina/Project_J_v6_30/Sample_v6_Amplicon_IDX2/v6_Amplicon_IDX2_CGATGT_L003_R2_001.fastq' in in_files_r2)
        in_files_r1_compare = r1file_name in in_files_r1

        self.assertTrue(in_files_r1_compare)
#        self.assertTrue(in_files_r2_compare)    

"""
    def split_files(self, compressed = False):
        def create_out_dir(self, dirname):
    def open_dataset_files(self):
    def close_dataset_files(self):
        def get_all_files(self):
        def create_inis(self):
        def perfect_reads(self):
        def uniq_fa(self):
    def open_write_close(self, ini_file_name, text):
        def get_fastq_file_names(self, f_input_file_path):
    def read1(self, files_r1, compressed):
    def read2(self, files_r2, compressed):
"""

if __name__ == '__main__':
    unittest.main()

