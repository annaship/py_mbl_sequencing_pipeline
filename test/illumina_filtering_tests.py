import unittest
import sys
import os
import shutil
sys.path.append("../")

import pipeline.illumina_filtering as ill_f

from pipeline.run import Run
import test.test_factory as fake_data_object

"""
to run: python pipeline/test/illumina_filtering_tests.py -v
"""

class IlluminaFilesTestCase(unittest.TestCase): 
    @classmethod  
    def setUpClass(cls):
        if os.path.exists("test/sample_data/illumina/result/20120614"):
            shutil.rmtree("test/sample_data/illumina/result/20120614")
        data_object = fake_data_object.data_object
        pi_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
        cls._runobj = Run(data_object, pi_path)    
        if not os.path.exists(cls._runobj.output_dir):
            os.mkdir(cls._runobj.output_dir)
        cls._illumina_filtering = ill_f.IlluminaFiltering(cls._runobj)

    @classmethod  
    def tearDownClass(cls):
        print "\nDone!"
        
    "Run setUp to clean db and fill out run info"
        
    def test_01_compare(self):
        result_true = self._illumina_filtering.compare(1, "<", 2)
        self.assertTrue(result_true)
        result_false = self._illumina_filtering.compare(1, ">", 2)
        self.assertEqual(result_false, False)
        

if __name__ == '__main__':
    unittest.main()

