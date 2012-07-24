import unittest
from  pipelineprocessor import process
from pipeline.run import Run
#from pipeline.runconfig import configDictionaryFromFile
from pipeline.validate import MetadataUtils
from pipeline.utils import convert_unicode_dictionary_to_str
import os
from pipeline.pipelinelogging import logger
import logging
import json
from Bio import SeqIO


class TestPipeline(unittest.TestCase):
    #BASE_OUTPUT = 'pipeline_output/' 
    def __init__(self):
        self.baseoutputdir = 'pipeline_output/' 
        self.platform = 'illumina'
        self.config_file_type = 'ini'
        
    def setUpForward(self):
        self.configPath = "test/data/trim_test_forward.ini"
        m=MetadataUtils(self)
        config_dict = m.create_dictionary_from_ini()
        print 'configDict',config_dict
        self.run = Run(config_dict, self, self.baseoutputdir)
        process(self.run,"trim")
        self.expected = self.get_expected_results('test/data/test_trim_forward.results')

    def setUpReverse(self):
        self.configPath = "test/data/trim_test_reverse.ini"
        m=MetadataUtils(self)
        config_dict = m.create_dictionary_from_ini()
        self.run = Run(config_dict, self, self.baseoutputdir)
        process(self.run,"trim")
        self.expected = self.get_expected_results('test/data/test_trim_reverse.results')

    def test_all(self):
        logger.setLevel(logging.DEBUG)    
        os.chdir("..")
        self.setUpForward()
        self.run_tests()
        self.setUpReverse()
        self.run_tests()        
    
    def run_tests(self):
        self.run_test_abundance_file()
        self.run_test_no_key_deletions_file()
        self.run_test_deleted_reasons_files()
        self.run_test_names_files()
        self.run_test_trimmed_files()
        self.run_test_unique_files()
        self.run_test_status_file()
        
    def run_test_status_file(self):
        trim_status_dict = convert_unicode_dictionary_to_str(json.loads(open(os.path.join(self.BASE_OUTPUT, os.path.join(self.run.run_date,"trim_status.txt")), 'r').read()))
        expected_trim_status = self.expected['trim_status']
        self.assertDictEqual(trim_status_dict, expected_trim_status)
        
    def run_test_unique_files(self):
        unique_dict = self.get_unique_info()
        self.assertDictEqual(unique_dict, self.expected['unique'])
        
    def run_test_trimmed_files(self):
        trimmed_dict = self.get_trimmed_info()
        self.assertDictEqual(trimmed_dict, self.expected['trimmed'])
        
    def run_test_names_files(self):
        names_dict = self.get_names_info()
        self.assertDictEqual(names_dict, self.expected['names'])
        
    def run_test_deleted_reasons_files(self):
        delete_reasons_dict = self.get_deleted_seqs_info()
        self.assertDictEqual(delete_reasons_dict, self.expected['delete_reasons'])
        
    def run_test_abundance_file(self):
        abundance_data = self.get_abundance_info()
        expected_abundance = self.expected['abundance']
        self.assertDictEqual(abundance_data, expected_abundance)
        
    def run_test_no_key_deletions_file(self):
        no_keys_by_lane_dict = self.get_no_key_deleted_info(self.BASE_OUTPUT + self.run.run_date)
        self.assertSetEqual(set(no_keys_by_lane_dict.keys()), set(self.expected['no_key_deletions']))

    def get_expected_results(self, expected_results_file):
        return convert_unicode_dictionary_to_str(json.loads(open(expected_results_file,'r').read()))

    def get_names_info(self):
        base_path = self.BASE_OUTPUT + self.run.run_date
        lane_keys = self.run.run_keys
        names_dict = {}
        for lane in lane_keys:
            names_dict[lane] = curr_lane_dict = {}
            no_key_file = open(os.path.join(base_path, lane + ".names"), 'r')
            for line in no_key_file:
                parts = line.strip().split('\t')
                curr_lane_dict[parts[0]] = parts[1].split(",")
        return names_dict
    
    def get_no_key_deleted_info(self, base_path):
        no_key_dict = {}
        no_key_file = open(os.path.join(base_path, "nokey.deleted.txt"), 'r')
        for line in no_key_file:
            parts = line.strip().split('\t')
            no_key_dict[parts[0]] = parts[1]
        return no_key_dict

    def get_deleted_seqs_info(self):
        base_path = self.BASE_OUTPUT + self.run.run_date
        lane_keys = self.run.run_keys
        deleted_dict = {}
        for lane in lane_keys:
            deleted_dict[lane] = curr_lane_dict = {}
            no_key_file = open(os.path.join(base_path, lane + ".deleted.txt"), 'r')
            for line in no_key_file:
                parts = line.strip().split('\t')
                curr_lane_dict[parts[0]] = parts[1]
        return deleted_dict
            
    def get_abundance_info(self):
        base_path = self.BASE_OUTPUT + self.run.run_date
        lane_keys = self.run.run_keys
        abundance_data = {}
        for key in lane_keys:            
            abundance_data[key] = key_dict = {}
            for record in SeqIO.parse(os.path.join(base_path, key) + ".abund.fa", 'fasta'):  
                parts = record.description.split(';')
                id = parts[0]
                count = int(parts[1].split('=')[1])
                seq = record.seq.tostring()
                key_dict[id] = {'size' : count, 'seq' : seq}
        return abundance_data
                
    def get_trimmed_info(self):
        base_path = self.BASE_OUTPUT + self.run.run_date
        lane_keys = self.run.run_keys
        trimmed_data = {}
        for key in lane_keys:            
            trimmed_data[key] = key_dict = {}
            for record in SeqIO.parse(os.path.join(base_path, key) + ".trimmed.fa", 'fasta'):  
                id = record.id
                seq = record.seq.tostring()
                key_dict[id] = seq
        return trimmed_data

    def get_unique_info(self):
        base_path = self.BASE_OUTPUT + self.run.run_date
        lane_keys = self.run.run_keys
        unique_data = {}
        for key in lane_keys:            
            unique_data[key] = key_dict = {}
            for record in SeqIO.parse(os.path.join(base_path, key) + ".unique.fa", 'fasta'):  
                id = record.id
                seq = record.seq.tostring()
                key_dict[id] = seq
        return unique_data


if __name__ == '__main__':
    unittest.main()
