# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

from pipeline.sample import Sample
from pipeline.configurationexception import ConfigurationException
import os
import constants as C
import ast

# read a config file and convert to a dictionary
def configDictionaryFromFile(config_file_path):
    import ConfigParser
    
    configDict = {}
    user_config = ConfigParser.ConfigParser()
    user_config.read(config_file_path)
    
    for section in user_config.sections():
        section_dict = configDict[section] = {}
        for option in user_config.options(section):
            section_dict[option] = user_config.get(section,option)

    return configDict


class RunConfig:
    """Doc string here."""
    def __init__(self, config_info, baseoutputdir, basepythondir):
        self.run_date   = None
        self.platform   = None # enum('454','illumina','ion_torrent','')
        self.input_dir  = None
        self.output_dir = None
        self.sff_files  = []
        self.run_keys = []
        self.run_key_lane_dict = {}
        self.samples = {}
        self.base_python_dir = basepythondir


        # if the config_info was a file path to an .ini file then convert to a dictionary
        # we'll take the info as an ini file or dictionary so we can be called by an api
        config_dict = config_info if (type(config_info)==dict) else configDictionaryFromFile(config_info)
        # now extract it all from the dictionary form
        self.initializeFromDictionary(config_dict)
        if 'vamps_user_upload' in config_info['general'] and config_info['general']['vamps_user_upload'] == True:
            pass 
            
        else:
            
    
            # primers should be in json format in a file and that file should be specified in the general section
            # print "curr dir: " + os.getcwd()
            # print "curr file all: " + os.path.realpath(__file__)
            # print "curr file dir: " + os.path.dirname(os.path.realpath(__file__))
            # primer_file = open(config_dict['general']['primer_file'])
            primer_file = open(os.path.join(self.base_python_dir, "config/mbl_primers.json"))
            ascii_primer_str = primer_file.read()
            self.primer_suites = ast.literal_eval(ascii_primer_str)
            
            # anchors should be similarly specified
            # anchor_json_text = open(config_dict['general']['anchor_file']).read()
            anchor_json_text = open(os.path.join(self.base_python_dir, "config/mbl_anchors.json")).read()
            self.anchors = ast.literal_eval(anchor_json_text)
    
        # this is our default output dir
        self.base_output_dir = baseoutputdir #user supplied or default
        self.output_dir = os.path.join(self.base_output_dir, self.run_date)
        # not sure if this setup should be here or in trim?  here for now
        self.trim_status_file_name = os.path.join(self.output_dir, 'trim_status.txt')
        self.trim_status_file_h = None #handle to file
        # not sure if this setup should be here or in chimera?  here for now
        self.chimera_status_file_name = os.path.join(self.output_dir, 'chimera_status.txt')
        self.chimera_status_file_h = None #handle to file
        
    

    # read a config dictionary and extract the info we want into the objects we use
    def initializeFromDictionary(self, configDict):
        # get the general stuff
        general_config = configDict['general']
        #if general_config['gast_data_source'] != 'database':
        self.run_date       = general_config['run_date']
        self.platform       = general_config.get('platform', "unknown")
        self.input_dir      = general_config.get('input_dir', None)
        self.require_distal = general_config.get('require_distal', True)
        self.minimumLength  = general_config.get('minimumLength', C.minimumLength)
        self.maximumLength  = general_config.get('maximumLength', C.maximumLength)
        self.minAvgQual     = general_config.get('minAvgQual',    C.minAvgQual)
        self.force_runkey   = general_config.get('force_runkey', None)
        
        # added gast_input_source for vamps uploads
        # so when users want to gast at a later time they will
        # look in the database and not the files (which may be missing)
        # see /xraid2-2/vampsweb/vampsdev/vamps_trim.py
        self.gast_input_source = 'files' # for regular gast pipeline
        if 'gast_input_source' in general_config: 
            self.gast_input_source = general_config['gast_input_source']
        
        #print general_config
        # parse out the input file info
        input_file_names  = [input_str.strip() for input_str in general_config['input_file_names'].split(',')]
        input_file_types  = [input_str.strip() for input_str in general_config['input_file_formats'].split(',')]
        
        if len(input_file_names) != len(input_file_types):
            raise Exception("Mismatch between the number of input_file_names(" + str(len(input_file_names)) + ") and input_file_types(" + str(len(input_file_types)) + ") in configuration information")
        
        lane_info = general_config['input_file_lanes'].strip()
        input_file_lanes  = [] if lane_info == '' else [input_str.strip() for input_str in lane_info.split(',')]

        # no lane info? better by our custom fasta-mbl format then
        if len(input_file_lanes) == 0 and len([  type for type in input_file_types if type != 'fasta-mbl' ]) > 0:
            raise Exception("Only fasta-mbl formatted sequence files are allowed to not provide a value for input_file_lanes")

        # if they give any lane information it then needs to either be 1 value (for all files) or match them exactly
        if len(input_file_lanes) > 1 and (len(input_file_names) != len(input_file_lanes)):
            raise Exception("Mismatch between the number of input_file_names(" + str(len(input_file_names)) + ") and lanes(" + str(len(input_file_lanes)) + ") in configuration information")
        
        self.input_file_info = {}
        for idx,input_file in enumerate(input_file_names):
            input_file_format = input_file_types[idx]
            if input_file_format not in ['sff', 'fasta', 'fasta-mbl', 'fastq', 'fastq-illumina', 'fastq-sanger']:
                raise Exception("Invalid sequence input file format: " + self.input_file_format)
            # make up a hash...they are allowed to not put in any input_file_lanes...could be 3 mbl fasta files which would all have lane
            # info encoded on each id/description line of the sequence record
            self.input_file_info[input_file] = {"name" : input_file, "format" : input_file_types[idx], "lane" : input_file_lanes[idx] if idx < len(input_file_lanes) else ""}
        
        # now deal with each lane_runkey combo (Sample) that is misnamed though
        # populate sample information for every run_key
        for lane_run_key in [s for s in configDict.keys() if s != 'general']:
            lane_run_dict = configDict[lane_run_key]
            sample = Sample(lane_run_key)
            # has defaults -not required
            try:
                sample.forward_primers = lane_run_dict['forward_primers'].split(',')
            except:
                sample.forward_primers = []
            try:
                sample.reverse_primers = lane_run_dict['reverse_primers'].split(',')
            except:
                sample.reverse_primers = []
            try:
                sample.stop_sequences = lane_run_dict['stop_sequences'].split(',')
            except:
                sample.stop_sequences = []
            try:
                sample.anchor = lane_run_dict['anchor']
            except:
                sample.anchor = ''
            # should we try to trim with mbl primers as well as custom ones
            try:
                sample.use_mbl_primers = lane_run_dict['use_mbl_primers']
            except:
                sample.use_mbl_primers = 1
            # required
            sample.direction = lane_run_dict['direction']                                                                   
            sample.project = lane_run_dict['project_name']
            sample.dataset = lane_run_dict['dataset_name']
            sample.dna_region = lane_run_dict['dna_region']
            sample.taxonomic_domain = lane_run_dict['taxonomic_domain']

            # a list of run_keys
            # convert: change ':' to '_'
            key = lane_run_key[:1]+'_'+lane_run_key[2:]
            self.run_keys.append(key)
            # a dictionary of samples
            self.samples[key] = sample
        
