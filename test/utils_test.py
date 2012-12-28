import unittest
import sys
import os
import shutil
sys.path.append("../")

import pipeline.utils as utils

from pipeline.run import Run
import pipeline.constants as C
import test.test_factory as fake_data_object

"""
to run: python pipeline/test/utils_tests.py -v
"""

class UtilsTestCase(unittest.TestCase): 
    @classmethod  
    def setUpClass(cls):          
        is_user_upload = False #we never call pipeline-ui.py to do vamps user upload.
        data_object = fake_data_object.data_object
        cls._dirs = utils.Dirs(is_user_upload, data_object['general']['run'], data_object['general']['platform'])    
        root_dir  = '/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test'
        cls.file_path = os.path.join(root_dir, data_object['general']['platform'], data_object['general']['run']) 
        if os.path.isdir(cls.file_path):
            shutil.rmtree(cls.file_path) 

    @classmethod  
    def tearDownClass(cls):
        print "\nDone!"
        
    """check_and_make_dir: if exists, yes, no, not exists
    check_dir
    get_path (is_user_upload? yes, no)
    check_and_make_output_dir
    create_all_output_dirs
    """
    
    def test_01_check_and_make_output_dir_if_not_exists(self):
        self._dirs.check_and_make_output_dir()
            
        self.assertEqual(self._dirs.output_dir, self.file_path)
        self.assertTrue(os.path.isdir(self._dirs.output_dir))

    def test_02_check_and_make_output_dir_if_exists(self):
        if not os.path.exists(self.file_path):
            os.mkdir(self.file_path)
        self._dirs.check_and_make_output_dir()
            
        self.assertEqual(self._dirs.output_dir, self.file_path)
        self.assertTrue(os.path.isdir(self._dirs.output_dir))
        
    def test_03_create_all_output_dirs(self):
        self._dirs.create_all_output_dirs()

        analysis_dir      = os.path.join(self._dirs.output_dir, C.subdirs['analysis_dir'])
        reads_overlap_dir = os.path.join(analysis_dir, C.subdirs['reads_overlap_dir'])

        self.assertEqual(self._dirs.analysis_dir, analysis_dir)
        self.assertTrue(os.path.isdir(self._dirs.analysis_dir))
        
        self.assertEqual(self._dirs.reads_overlap_dir, reads_overlap_dir)
        self.assertTrue(os.path.isdir(self._dirs.reads_overlap_dir))


if __name__ == '__main__':
    unittest.main()

