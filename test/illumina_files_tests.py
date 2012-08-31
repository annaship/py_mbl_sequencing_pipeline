import unittest
import sys
import os
import shutil
sys.path.append("../")

import pipeline.illumina_files as ill_f

from pipeline.run import Run
import test.test_factory as fake_data_object

"""
to run: python pipeline/test/illumina_files_tests.py -v
"""

class IlluminaFilesTestCase(unittest.TestCase): 
    @classmethod  
    def setUpClass(cls):
        if os.path.exists("test/sample_data/illumina/result/20120614"):
            shutil.rmtree("test/sample_data/illumina/result/20120614")
        data_object = fake_data_object.data_object
        pi_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
        cls._runobj = Run(data_object, pi_path)    
        cls._illumina_files = ill_f.IlluminaFiles(cls._runobj)
        cls.file_path = "test/sample_data/illumina/result/20120614/analysis"
        cls.dataset_emails = fake_data_object.dataset_emails

    @classmethod  
    def tearDownClass(cls):
        print "\nDone!"
        
    "Run setUp to clean db and fill out run info"
        
    def test_01_create_out_dir(self):
        file_path  = self._illumina_files.create_out_dir(os.path.join(self._runobj.output_dir, "analysis"))
        self.assertEqual(file_path, self.file_path)
        
#    def test_02_open_dataset_files(self):
#        self._illumina_files.open_dataset_files()
#        print len([name for name in os.listdir(self.file_path) if os.path.isfile(name)])

    def test_03_get_all_files(self):
        res = self._illumina_files.get_all_files()
        self.assertEqual(res['test/sample_data/illumina/result/20120614/analysis/unknown.fastq'], ('test/sample_data/illumina/result/20120614/analysis/unknown', '.fastq'))

    def test_04_create_inis(self):
        self._illumina_files.dataset_emails = self.dataset_emails
        self._illumina_files.create_inis()
        ini_files = len([f for f in os.listdir(self.file_path) if f.endswith('.ini') and os.path.isfile(os.path.join(self.file_path, f))])
        self.assertEqual(ini_files, 10)
        
        with open(os.path.join(self.file_path, "SMPL90_3.ini")) as file:
            for line in file.readlines():
                if line.strip() == "project_name = SMPL90_3":
                    print "\nValid ini\n"

    def test_04_perfect_reads(self):
        self._illumina_files.dataset_emails = self.dataset_emails
        self._illumina_files.perfect_reads()
        f_path = os.path.join(self.file_path, "perfect_reads")
        files_amount  = len([name for name in os.listdir(f_path) if os.path.isfile(os.path.join(f_path, name))])
        self.assertEqual(files_amount, 70)
        
    def test_05_uniq_fa(self):
        self._illumina_files.uniq_fa()
        f_path = os.path.join(self.file_path, "perfect_reads")
        files_amount  = len([name for name in os.listdir(f_path) if os.path.isfile(os.path.join(f_path, name))])
        uniq_files = len([f for f in os.listdir(f_path) if f.endswith('-PERFECT_reads.fa.unique') and os.path.isfile(os.path.join(f_path, f))])       
        self.assertEqual(files_amount, 90)
        self.assertEqual(uniq_files, 10)
        
    def test_06_get_fastq_file_names(self):
        (in_files_r1, in_files_r2) = self._illumina_files.get_fastq_file_names(self._runobj.input_dir)
#        in_files_r1 = ['./results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX1_ATCACG_L003_R1_001.filtered.fastq', './results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX1_ATCACG_L003_R1_001.filtered.fastq.failed', './results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX2_CGATGT_L003_R1_001.filtered.fastq', './results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX2_CGATGT_L003_R1_001.filtered.fastq.failed']

        in_files_r1_compare = './results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX1_ATCACG_L003_R1_001.filtered.fastq' in in_files_r1
#        ['./test/sample_data/illumina/Project_J_v6_30/Sample_v6_Amplicon_IDX1/v6_Amplicon_IDX1_ATCACG_L003_R1_001.fastq', './test/sample_data/illumina/Project_J_v6_30/Sample_v6_Amplicon_IDX2/v6_Amplicon_IDX2_CGATGT_L003_R1_001.fastq']    
        in_files_r2_compare = ('./results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX1_ATCACG_L003_R2_001.filtered.fastq' in in_files_r2) or ('./results/illumina_filtering/123001/illumina_filtered/v6_Amplicon_IDX1_ATCACG_L003_R2_001.fastq' in in_files_r2)

        self.assertTrue(in_files_r1_compare)
        self.assertTrue(in_files_r2_compare)    

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

