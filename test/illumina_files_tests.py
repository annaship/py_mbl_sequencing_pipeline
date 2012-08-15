import unittest
import sys
import os
import shutil
sys.path.append("../")

import pipeline.illumina_files as ill_f

from pipeline.run import Run
import test.test_factory as fake_data_object



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


    

"""
    def split_files(self, compressed = False):
    def create_out_dir(self, dirname):
    def open_dataset_files(self):
    def close_dataset_files(self):
    def get_all_files(self):
    def perfect_reads(self):
    def uniq_fa(self):
    def create_inis(self):
    def open_write_close(self, ini_file_name, text):
    def get_fastq_file_names(self, f_input_file_path):
    def read1(self, files_r1, compressed):
    def read2(self, files_r2, compressed):
"""

if __name__ == '__main__':
    unittest.main()

