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
from pipeline.validate import CSV_utils
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

def configDictionaryFromFile_csv(config_file_path, args):
    
    
    data = {}
    projects = {}
    megadata = {}
    megadata['general'] = {}
    megadata['general']['input_dir'] = args.input_dir
    megadata['general']['configFile'] = config_file_path
    #megadata['general']['output_dir'] = self.args.output_dir
    megadata['general']['platform'] = args.platform
    megadata['general']['run'] = args.run
    megadata['general']['run_date'] = args.run
    #megadata['general']['run'] = self.args.run
    megadata['general']["input_file_format"] = args.input_file_format
    #input_dir,"/xraid2-2/sequencing/Illumina/20120525_recalled/Project_Sandra_v6/analysis/"
    megadata['general']["input_file_suffix"] = args.input_file_suffix
    
    # HERE make list of input_file_names
    # open dir
    # and list of types and lanes
    if 'input_file_suffix' not in megadata['general']:
            megadata['general']['input_file_suffix'] = ''
    file_count = 0
    files_list = []
    if os.path.isdir(megadata['general']['input_dir']):
        p = megadata['general']['input_dir'], '*'+megadata['general']['input_file_suffix']
        print p
        for infile in glob.glob( os.path.join(megadata['general']['input_dir'], '*'+megadata['general']['input_file_suffix']) ):
            files_list.append(os.path.basename(infile))
            file_count += 1
    else:
        sys.exit("ERROR:no input directory or directory permissions problem")
        
    if not file_count:
        sys.exit("ERROR:No files were found in '"+megadata['general']['input_dir']+"' with a suffix of '"+megadata['general']['input_file_suffix']+"'")
        
    megadata['general']['input_file_names'] = ','.join(files_list)
    megadata['general']['input_file_formats'] = ','.join([megadata['general']['input_file_format'] for i in files_list])
    # assign 1 to lanes -- kludge
    megadata['general']['input_file_lanes'] = ','.join(['1']*file_count)
        
        
    test_datasets = {}
    dataset_counter = {}
    headers = ''
    # must be comma sep
    
    f_in_md = open(config_file_path, 'r')
    # must be comma sep
    lines = f_in_md.readlines()
    for line in lines:
        
        line = line.strip()
            
        if not line:
            continue
        lst = [i.strip('"').replace(" ", "_") for i in line.strip().split(',')]
        
        
        if not lst[0]:
            continue
        temp = {}   
        if not headers:
            headers = [i.strip('"').lower().replace(" ", "_") for i in line.split(',')]

            if sorted(known_header_list) != sorted(headers):
                sys.exit("ERROR : unknown_headers:\nyours: "+ ' '.join(sorted(headers))+"\nours:  "+' '.join(sorted(known_header_list)))
        else:
            for n in range(0,len(headers)):
                #print headers[n], lst[n]
                try:
                    temp[headers[n]] = lst[n]
                except:
                    sys.exit("ERROR:It looks like the header count and the data column count are different.")
        
            temp['file_prefix'] = temp['dataset']+'_'+temp['barcode'].upper().replace('N','')
        
            #data[lst[0]] = temp
            
            
            unique_identifier = temp['barcode_index']+'_'+temp['run_key']+'_'+temp['lane']
            megadata[unique_identifier]={}
            if unique_identifier in test_datasets:
                sys.exit("ERROR: duplicate run_key:barcode_index:lane: "+unique_identifier+" - Exiting")
            else:                     
                test_datasets[unique_identifier] = 1
                
            megadata[unique_identifier]['dataset'] = temp['dataset']
            megadata[unique_identifier]['project'] = temp['project']
            
            if temp['project'] in dataset_counter:
                dataset_counter[temp['project']] += 1
            else:
                dataset_counter[temp['project']] = 1
            
            #megadata[unique_identifier]['ds_count'] = 1
            megadata[unique_identifier]['project']      = temp['project']
            megadata[unique_identifier]['run_key']      = temp['run_key']
            megadata[unique_identifier]['lane']         = temp['lane']
            megadata[unique_identifier]['tubelabel']    = temp['tubelabel']
            megadata[unique_identifier]['barcode']      = temp['barcode']
            megadata[unique_identifier]['adaptor']      = temp['adaptor']
            megadata[unique_identifier]['dna_region']   = temp['dna_region']
            megadata[unique_identifier]['amp_operator'] = temp['amp_operator']
            megadata[unique_identifier]['seq_operator'] = temp['seq_operator']
            megadata[unique_identifier]['barcode_index']= temp['barcode_index']
            megadata[unique_identifier]['overlap']      = temp['overlap']
            megadata[unique_identifier]['insert_size']  = temp['insert_size']
            megadata[unique_identifier]['file_prefix']  = temp['file_prefix']
            megadata[unique_identifier]['read_length']  = temp['read_length']
            megadata[unique_identifier]['primer_suite'] = temp['primer_suite']
    return megadata

class RunConfig:
    """Doc string here."""
    def __init__(self, config_info, args, basepythondir):
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


        # if the config_info was a file path to an .csv file then convert to a dictionary
        # we'll take the info as an ini file or dictionary so we can be called by an api
        # ie vamps user uploads: the config info is a dictionary
    	if type(config_info)==dict:
            config_dict = config_info
    	elif args.platform == 'illumina' and args.config_file_type == 'csv':
            #config_dict = configDictionaryFromFile_csv(config_info, args)
            v = CSV_utils()
            # read the csv config file
            my_csv = readCSV(file_path = args.configPath)
            config_dict = v.create_dictionary_from_csv(args, my_csv)
        else:
            sys.exit("Unknown platform and configFile type for dictionary conversion")
            
#         
#         Validate here and return the dict for both ini and csv
#             
#     	v = Validate(config_info) # either a file path or dict
#     	if self.config_file_type == 'csv':
#         	config_dict = v.validate_csv()        
#         elif self.config_file_type == 'ini':
#         	config_dict = v.validate_ini()
#         elif self.config_file_type == 'dict':
#         	config_dict = v.validate_dict()
#         else:
#              sys.exit("could not determine config type: "+self.config_file_type)
        
        # if the config_info was a file path to an .ini file then convert to a dictionary
        # we'll take the info as an ini file or dictionary so we can be called by an api
        # ie vamps user uploads: the config info is a dictionary

#        config_dict = config_info if (type(config_info)==dict) else configDictionaryFromFile_ini(config_info)
        # now extract it all from the dictionary form
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
    
        # this is our default output dir -- Always rundate?
        self.base_output_dir = os.path.normpath(args.baseoutputdir) #user supplied or default
        self.output_dir = os.path.join(self.base_output_dir, self.run_date)
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
        
        #print general_config
        # parse out the input file info
        if 'files_list' in general_config:
            input_file_names = general_config['files_list']
        else:
            input_file_names  = [input_str.strip() for input_str in general_config['input_file_names'].split(',')]
        
        if 'file_formats_list' in general_config:    
            input_file_types = general_config['file_formats_list']
        else:
            input_file_types  = [input_str.strip() for input_str in general_config['input_file_formats'].split(',')]
        
        if len(input_file_names) != len(input_file_types):
            raise Exception("Mismatch between the number of input_file_names(" + str(len(input_file_names)) + ") and input_file_types(" + str(len(input_file_types)) + ") in configuration information")
        
        if 'lanes_list' in general_config: 
            input_file_lanes = general_config['lanes_list']
        else:        
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
            if input_file_format not in C.input_file_formats:
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
            
            
            
            if self.platform == 'illumina':
                # req specifically for illumina
                sample.barcode_index = lane_run_dict['barcode_index'] 
                sample.overlap = lane_run_dict['overlap'] 
                sample.read_length = lane_run_dict['read_length'] 
                sample.file_prefix = lane_run_dict['file_prefix'] 
                sample.insert_size = lane_run_dict['insert_size']
                # concatenate: barcode_index and run_key and lane
                key = lane_run_dict['barcode_index'] +'_'+ lane_run_dict['run_key'] +'_'+ lane_run_dict['lane'] 
                self.run_keys.append(key)  
                
            elif self.platform == '454':
                # required for 454
                sample.direction = lane_run_dict['direction'] 
                sample.taxonomic_domain = lane_run_dict['taxonomic_domain']
                # a list of run_keys
                # convert: change ':' to '_'
                key = lane_run_key[:1]+'_'+lane_run_key[2:]
                self.run_keys.append(key)
                
            sample.project = lane_run_dict['project']
            sample.dataset = lane_run_dict['dataset']
            sample.dna_region = lane_run_dict['dna_region']           

            
            # a dictionary of samples
            self.samples[key] = sample
        
