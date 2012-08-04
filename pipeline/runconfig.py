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
from pipeline.metadata import MetadataUtils
import sys,os
import glob
import constants as C
import ast
from pipeline.get_ini import readCSV

# read a config file and convert to a dictionary
def configDictionaryFromFile_ini(config_file_path):
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
    def __init__(self, config_info, args, basepythondir):
        self.configPath = config_info
        self.args       = args
        self.run_date   = None
        self.platform   = None # enum('454','illumina','ion_torrent','')
        self.input_dir  = None
        self.output_dir = None
        self.config_file_type = args.config_file_type  # ini, csv or dict
        self.sff_files  = []
        self.run_keys = []
        self.run_key_lane_dict = {}
        self.samples = {}
        self.base_python_dir = os.path.normpath(basepythondir)
        self.configFile = config_info
        #
        # IMPORTANT to get a dictionary here from whatever the input is:
        #     platform and configFile type
        # if the config_info was a file path to an .csv or .ini file then convert to a dictionary
        # 
        # for vamps user uploads: the config info is a dictionary
        v = MetadataUtils(args)
    	if type(config_info)==dict:
            config_dict = config_info
                 
        elif self.args.platform == 'illumina':
            
            config_dict = v.create_dictionary_from_ini()
            config_dict['general'] = v.get_command_line_items(config_dict['general'])
            config_dict['general']['config_file'] = os.path.join(config_dict['general']['output_dir'], config_dict['general']['run']+'.ini')
            config_dict['general']['status_file'] = os.path.join(config_dict['general']['output_dir'], 'STATUS.txt')
            config_dict['general']['files_list'] = config_dict['general']['input_files'].split(',')
            #config_dict = v.check_for_input_files(config_dict) 
            
        elif self.args.platform == '454':
        
            config_dict = v.create_dictionary_from_ini()
            config_dict['general'] = v.get_command_line_items(config_dict['general'])
            config_dict['general']['config_file'] = os.path.join(config_dict['general']['output_dir'], config_dict['general']['run']+'.ini')
            config_dict['general']['status_file'] = os.path.join(config_dict['general']['output_dir'], 'STATUS.txt')
            #config_dict = v.check_for_input_files(config_dict)
            
        elif args.platform == 'ion_torrent':
            sys.exit("3-ConfigFile conversion to dictionary not written yet for platform ("+self.args.platform+") ")
        
        
        else:
            sys.exit("Unknown platform for dictionary conversion")
            
        
        # if the config_info was a file path to an .ini file then convert to a dictionary
        # we'll take the info as an ini file or dictionary so we can be called by an api
        # ie vamps user uploads: the config info is a dictionary

        #
        # now extract it all from the dictionary form
        #
        #
        
        self.initializeFromDictionary(config_dict)
         
         
        if type(config_info)==dict and 'vamps_user_upload' in config_info['general'] and config_info['general']['vamps_user_upload'] == True:
            self.vamps_user_upload = True            
        else:
            self.vamps_user_upload = False
    
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
    
        
        # NOTE
        # the base ouput directory is best gotten from the command line (no default)
        # second best is from the metadata file: ini or csv
        # if neither of these are set then output to current directory
        # and always attach the rundate dir to it
        if args.baseoutputdir:
            self.base_output_dir = os.path.normpath(args.baseoutputdir) #user supplied or default 
        elif config_dict['general']['output_dir']:
            self.base_output_dir = os.path.normpath(config_dict['general']['output_dir'])
        else:
            self.base_output_dir = '.'
        # this is our default output dir -- Always rundate?
        self.output_dir = os.path.join(self.base_output_dir, self.run_date)
        #self.output_dir = os.path.join(config_dict['general']['output_dir'])
        self.run_status_file_name = os.path.join(self.output_dir,"STATUS.txt")
        self.run_status_file_h = None #handle to file
                  
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
        
        if 'files_list' in general_config:
            input_file_names = general_config['files_list']
        else:
            input_file_names  = [input_str.strip() for input_str in general_config['input_file_names'].split(',')]
        
#         
#         # for ini file:  (no plurals)
#         # 1) if input_file_format is a comma sep list then it should match the count of input_file_name
#         #       The same with input_file_lane
#         # 2) if input_file_format is supplied and is a single item it will apply to all the input files
#         #       either in input_dir or the list (or single) of input_file_name
#         # 3) EITHER input_dir OR input_file_name will be supplied (but not both)
#         #
#         if self.platform == '454':
#             
#             if 'input_file_format' in general_config and general_config['input_file_format'] != '':
#                 input_file_types = general_config['input_file_format']
#             elif 'file_formats_list' in general_config:    
#                 input_file_types = general_config['file_formats_list']
#             else:
#                 input_file_types  = [input_str.strip() for input_str in general_config['input_file_formats'].split(',')]
#             
#             print 'input_file_types= ',input_file_types
#             if len(input_file_names) != len(input_file_types):
#                 raise Exception("Mismatch between the number of input_file_names(" + str(len(input_file_names)) + ") and input_file_types(" + str(len(input_file_types)) + ") in configuration information")
#             
#             if 'lanes_list' in general_config: 
#                 input_file_lanes = general_config['lanes_list']
#             else:        
#                 lane_info = general_config['input_file_lanes'].strip()
#                 input_file_lanes  = [] if lane_info == '' else [input_str.strip() for input_str in lane_info.split(',')]
#     
#             # no lane info? better by our custom fasta-mbl format then
#             if len(input_file_lanes) == 0 and len([  type for type in input_file_types if type != 'fasta-mbl' ]) > 0:
#                 raise Exception("Only fasta-mbl formatted sequence files are allowed to not provide a value for input_file_lanes")
#     
#             # if they give any lane information it then needs to either be 1 value (for all files) or match them exactly
#             if len(input_file_lanes) > 1 and (len(input_file_names) != len(input_file_lanes)):
#                 raise Exception("Mismatch between the number of input_file_names(" + str(len(input_file_names)) + ") and lanes(" + str(len(input_file_lanes)) + ") in configuration information")
#         else:
#             input_file_types = []   
#             input_file_lanes = []
#         
#         
#         
#         
 
 
        self.input_file_info = {}
#        print general_config
        for idx,input_file in enumerate(input_file_names):
            
            if "input_file_format" in general_config:
                file_format = general_config['input_file_format']
            else:
                # default
                file_format = 'fasta'
            
            
            if file_format not in C.input_file_formats:
                raise Exception("Invalid sequence input file format: " + self.args.input_file_format)
                
            if "input_file_lane" in general_config:
                file_lane = general_config['input_file_lane']
            else:
                # default
                file_lane = ''    
                
            # make up a hash...they are allowed to not put in any input_file_lanes...could be 3 mbl fasta files which would all have lane
            # info encoded on each id/description line of the sequence record
            
            self.input_file_info[input_file] =  {  "name" : input_file, 
                                                   "format" : file_format, 
                                                   "lane" : file_lane
                                                }
        
        
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
#################################
            try:
                sample.run_key = lane_run_dict['run_key']
            except:
                sample.run_key = ''
            try:
                sample.lane = lane_run_dict['lane']
            except:
                sample.lane = ''
            try:
                sample.adaptor = lane_run_dict['adaptor']
            except:
                sample.adaptor = ''
            try:
                sample.barcode = lane_run_dict['barcode']
            except:
                sample.barcode = ''
            try:
                sample.seq_operator = lane_run_dict['seq_operator']
            except:
                sample.seq_operator = ''
            try:
                sample.amp_operator = lane_run_dict['amp_operator']
            except:
                sample.amp_operator = ''
            try:
                sample.primer_suite = lane_run_dict['primer_suite']
            except:
                sample.primer_suite = ''
            try:
                sample.tubelabel = lane_run_dict['tubelabel']
            except:
                sample.tubelabel = ''
            try:    
                sample.dna_region = lane_run_dict['dna_region'] 
            except:
                sample.dna_region = ''
                
            sample.data_owner           = lane_run_dict['data_owner']
            sample.first_name           = lane_run_dict['first_name']
            sample.last_name            = lane_run_dict['last_name']
            sample.email                = lane_run_dict['email']
            sample.institution          = lane_run_dict['institution']
            sample.project_title        = lane_run_dict['project_title']
            sample.project_description  = lane_run_dict['project_description']
            sample.funding              = lane_run_dict['funding']
            sample.env_sample_source    = lane_run_dict['env_sample_source']
            sample.dataset_description  = lane_run_dict['dataset_description']
                
            if self.platform == 'illumina':
                # req specifically for illumina
                sample.barcode_index = lane_run_dict['barcode_index'] 
                sample.overlap = lane_run_dict['overlap'] 
                sample.read_length = lane_run_dict['read_length'] 
                sample.file_prefix = lane_run_dict['file_prefix'] 
                sample.insert_size = lane_run_dict['insert_size']
                # concatenate: barcode_index and run_key and lane
                key = lane_run_dict['barcode_index'] +'_'+ lane_run_dict['run_key'] +'_'+ lane_run_dict['lane'] 
                #sample.key = key
                self.run_keys.append(key)  
                
            elif self.platform == '454':
                # required for 454
                sample.direction = lane_run_dict['direction'] 
                sample.taxonomic_domain = lane_run_dict['domain']
                # a list of run_keys
                # convert: change ':' to '_'
                key = lane_run_key[:1]+'_'+lane_run_key[2:]
                #sample.key = key
                self.run_keys.append(key)
                
            sample.project = lane_run_dict['project']
            sample.dataset = lane_run_dict['dataset']
                      

            
            # a dictionary of samples
            self.samples[key] = sample
        
