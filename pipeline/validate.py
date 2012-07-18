#!/usr/bin/env python

##!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

import subprocess
import sys, os, stat
import glob
import time
import shutil
import random
from pipeline.pipelinelogging import logger
import constants as C
#import pipeline.fastalib


class FastaReader:
    def __init__(self,file_name=None):
        self.file_name = file_name
        self.h = open(self.file_name)
        self.seq = ''
        self.id = None
        self.revcomp_seq = None
        self.base_counts = None
        
    def next(self): 
        def read_id():
            return self.h.readline().strip()[1:]

        def read_seq():
            ret = ''
            while True:
                line = self.h.readline()
                
                while len(line) and not len(line.strip()):
                    # found empty line(s)
                    line = self.h.readline()
                
                if not len(line):
                    # EOF
                    break
                
                if line.startswith('>'):
                    # found new defline: move back to the start
                    self.h.seek(-len(line), os.SEEK_CUR)
                    break
                    
                else:
                    ret += line.strip()
                    
            return ret
        
        self.id = read_id()
        self.seq = read_seq()
        
        
        if self.id:
            return True    

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

def validate_dict():
    """
    This is only used for data that comes in as a dictionary rather than a file
    such as with vamps user uploads
    """
    print "Validating input dictionary -TODO"
    # must be a general section
    # should I create a dict here??? -That would render much code in
    #    runconfig useless.
    # are we going to continue developing ini style config files if
    #   no one uses them?  
    configDict = config_info
    print "Finished Validating"
    return configDict
    
    
class MetadataUtils:
    """
    Class to read metadata files (csv and ini style)
    validate and create a dictionary from them
    """
    Name = "Metadata_utils"
    def __init__(self, run=None):
        self.run = run
        self.known_header_list = C.csv_header_list
        self.primer_suites     = C.primer_suites 
        self.dna_regions       = C.dna_regions
        

    def create_dictionary_from_illumina_csv(self, args, my_csv):
    
        data_object = self.populate_data_object_illumina(args, my_csv)
        
        data_object = self.check_for_input_files(data_object)
        
        return data_object
    def create_dictionary_from_454_csv(args, my_csv):
        data_object = self.populate_data_object_454(args, my_csv)
        
        data_object = self.check_for_input_files(data_object)
        
        return data_object
    
    def validate_454_csv(args, my_csv):
        print "TODO: write validate_454_csv def"
        pass
    def validate_illumina_csv(self, args, my_csv):
    
        data_object = self.populate_data_object_illumina(args, my_csv)
                      
       

#        print 'general:', megadata['general']
        
        
        
        # start error checking here
        # MUST be in list: "Domain","Primer Suite","DNA Region"
        # MUST MATCH: "Domain","Primer Suite","DNA Region"
        #
        #
        # VAMPS project name format:  SLM_GCB_Bv6
        #
        #
        data_object = self.check_for_input_files(data_object)
        
        
        
        for x in data_object['general']:
            print "%s = %s" % (x, data_object['general'][x])
        
        
        
        """
        TODO:
            1) split into methods
            2) we can use "content" variable instead of megadata here
        """       
        #print dataset_counter
        for item in data_object:
            if item != 'general':            
                self.check_for_datasets(data_object[item])
                self.check_for_missing_values(data_object[item])
                self.check_domain_suite_region(data_object[item])
                self.check_project_name(data_object[item])





        """
            TODO:
                 other checks to put in:
                 check for duplicate dataset name:  NO
                 that data == file prefix
                     if we have an input directory that each dataset has a coresponding file - for illumina
                 Missing data is ok for barcode and adaptor (illumina only - for env454 also - A.)
                 
            
        """
                #print item,megadata[item],"\n\n"
        print "SUCCESS: Finished Validating"        

    def populate_data_object_454(self, args, my_csv):
        print "TODO: write populate_data_object_454 def"
        pass
        
        
    def populate_data_object_illumina(self, args, my_csv):
        data = {}
        data['general'] = {}
        test_datasets = {}
        dataset_counter = {}
        headers = ''
        if self.run:
            infile = self.configPath
            data['general']['input_dir'] = self.run.input_dir
            #megadata['general']['output_dir'] = self.args.output_dir
            data['general']['platform'] = self.run.platform
            data['general']['run'] = self.run.run_date
            data['general']['run_date'] = self.run.run_date
            #megadata['general']['run'] = self.args.run
            data['general']["input_file_format"] = self.run.input_file_format
            #input_dir,"/xraid2-2/sequencing/Illumina/20120525_recalled/Project_Sandra_v6/analysis/"
            data['general']["input_file_suffix"] = self.run.input_file_suffix
        else:            
            infile = args.configPath
            data['general']['input_dir'] = args.input_dir
            #megadata['general']['output_dir'] = self.args.output_dir
            data['general']['platform'] = args.platform
            data['general']['run'] = args.run
            data['general']['run_date'] = args.run
            #megadata['general']['run'] = self.args.run
            data['general']["input_file_format"] = args.input_file_format
            #input_dir,"/xraid2-2/sequencing/Illumina/20120525_recalled/Project_Sandra_v6/analysis/"
            data['general']["input_file_suffix"] = args.input_file_suffix
            
        print "Validating csv type ConfigFile"
        
        # changes spaces to '_' and all lowercase

        temp = {}   

        
#        my_read_csv = readCSV(file_path = infile)
#        my_read_csv.put_run_info()
#        print "content[1].keys(): "
#        print content[1].keys()
#        # To see the list of statistics available for each line
#        for k, v in content.items():
#            print k, v['dataset'], v 
        content     = my_csv.read_csv()
        headers     = content[1].keys()
        headers_clean = [x.strip('"').replace(" ", "_").lower() for x in headers]
        if self.check_headers(headers_clean):

#
#                try:
#                    temp[headers[n]] = lst[n]
#                except:
#                    sys.exit("ERROR:It looks like the header count and the data column count are different.")
            for k, v in content.items():
                run_key = v['run_key'].replace('N','').upper()
                temp['file_prefix'] = v['dataset']+'_'+ run_key
#                print "v = %s\n" % v
#                v = {'barcode_index': 'ATCACG', 'project': 'JCR_SPO_Bv6', 'lane': '3', 'run': '20120613', 'dna_region': 'v6', 'adaptor': '', 
#                      'barcode': '', 'seq_operator': 'JV', 'overlap': 'complete', 'dataset': 'H40', 'run_key': 'NNNNACGCA', 'read_length': '101', 
#                       'file_prefix': 'H40', 'data_owner': 'jreveillaud', 'primer_suite': 'Bacterial v6 Suite', 'tubelabel': 'H40', 'amp_operator': 'JR', 'insert_size': '230'}; 
#                        temp['file_prefix'] = H40_
                unique_identifier   = v['barcode_index']+'_'+run_key+'_'+v['lane']
                data[unique_identifier] = {}
                if unique_identifier in test_datasets:
                    sys.exit("ERROR: duplicate run_key:barcode_index:lane: "+unique_identifier+" - Exiting")
                else:
                    test_datasets[unique_identifier] = 1
#                print "test_datasets = %s;\ntemp['file_prefix'] = %s\nunique_identifier = %s" % (test_datasets,temp['file_prefix'], unique_identifier)
                
                data[unique_identifier]['dataset'] = v['dataset']
                data[unique_identifier]['project'] = v['project']
                
                if v['project'] in dataset_counter:
                    dataset_counter[v['project']] += 1
                else:
                    dataset_counter[v['project']] = 1
                
                #megadata[unique_identifier]['ds_count'] = 1
                data[unique_identifier]['project']              = v['project']
                data[unique_identifier]['run_key']              = v['run_key']
                data[unique_identifier]['lane']                 = v['lane']
                data[unique_identifier]['tubelabel']            = v['tubelabel']
                data[unique_identifier]['barcode']              = v['barcode']
                data[unique_identifier]['adaptor']              = v['adaptor']
                data[unique_identifier]['dna_region']           = v['dna_region']
                data[unique_identifier]['amp_operator']         = v['amp_operator']
                data[unique_identifier]['seq_operator']         = v['seq_operator']
                data[unique_identifier]['barcode_index']        = v['barcode_index']
                data[unique_identifier]['overlap']              = v['overlap']
                data[unique_identifier]['insert_size']          = v['insert_size']
                data[unique_identifier]['file_prefix']          = v['file_prefix']
                data[unique_identifier]['read_length']          = v['read_length']
                data[unique_identifier]['primer_suite']         = v['primer_suite']
                data[unique_identifier]['first_name']           = v['first_name']
                data[unique_identifier]['last_name']            = v['last_name']
                data[unique_identifier]['email']                = v['email']
                data[unique_identifier]['institution']          = v['institution']
                data[unique_identifier]['project_title']        = v['project_title']
                data[unique_identifier]['project_description']  = v['project_description']
                data[unique_identifier]['funding']              = v['funding']
                data[unique_identifier]['env_sample_source']    = v['env_sample_source']
                data[unique_identifier]['dataset_description']  = v['dataset_description']
        for item in data:
            if item != 'general':
                data[item]['primer_suite']  = data[item]['primer_suite'].lower().replace(" ", "_")
                data[item]['dna_region']    = data[item]['dna_region'].lower().replace(" ", "_")
                data[item]['barcode']       = data[item]['barcode'].upper()
                data[item]['barcode_index'] = data[item]['barcode_index'].upper()
                data[item]['ds_count']      = str(dataset_counter[data[item]['project']])
        return data
        
    def check_for_input_files(self,data_object):
    
        file_count = 0
        files_list = []
        imports_list = []
        lanes_list = []
        if os.path.isdir(data_object['general']['input_dir']):
            p = data_object['general']['input_dir'], '*'+data_object['general']['input_file_suffix']
            
            for infile in glob.glob( os.path.join(data_object['general']['input_dir'], '*'+data_object['general']['input_file_suffix']) ):
                files_list.append(os.path.basename(infile))
                for x in data_object:
                    if 'file_prefix' in data_object[x]:
                        #print data_object[x]['file_prefix']
                        
                        if os.path.basename(infile).split('-')[0] == data_object[x]['file_prefix']:
                            lanes_list.append(data_object[x]['lane'])
                        
                file_count += 1
        else:
            sys.exit("ERROR: no input directory or directory permissions problem: "+data_object['general']['input_dir'])
            
        if not file_count:
            sys.exit("ERROR: No files were found in '"+data_object['general']['input_dir']+"' with a suffix of '"+data_object['general']['file_suffix']+"'")
        
        data_object['general']['files_list'] = files_list
        
        data_object['general']['file_count'] = file_count
        # all the files in an illumina directory should be the same type
        data_object['general']['file_formats_list'] = [data_object['general']["input_file_format"]] * file_count
        data_object['general']['lanes_list'] = lanes_list
        
        
        
        return data_object
        
        
        
    def check_for_datasets(self,item):
        if not item['dataset']:
            sys.exit("ERROR:Current dataset name is missing or corrupt - Exiting")
            
    def check_for_missing_values(self,item):
        missing_key   = ''
        missing_value = ''
        for k,v in item.iteritems():
            if not k:
                #sys.exit("ERROR: key for: '"+v+"' is missing or corrupt - Exiting")
                logger.warn("ERROR: key for: '"+v+"' is missing or corrupt - Continuing")
                missing_key = v
            if not v:
                if (k != 'barcode' and k != 'adaptor'): #these could be empty
                    logger.warn("ERROR: value of: '"+k+"' is missing or corrupt - Continuing")
                    missing_value = k
        
        if missing_key:
            #sys.exit("ERROR: value of: "+missing_key+" is missing or corrupt - Exiting")
            pass
        if missing_value:
            #sys.exit("ERROR: value of: "+missing_value+" is missing or corrupt - Exiting")
            pass
        
    def check_domain_suite_region(self,item):       
        
        # CHECK MUST MATCH: "Domain","Primer Suite","DNA Region"
        if item['primer_suite'] not in self.primer_suites:
            sys.exit("ERROR: Primer Suite not found: "+item['primer_suite'])
        #if dataset_items['domain'] not in domains:
        #   sys.exit("ERROR: Domain not found: "+dataset_items['domain'])
        if item['dna_region'] not in self.dna_regions:
            sys.exit("ERROR: DNA Region not found: "+item['dna_region'])
        # "Bacterial v6","BacterialV6Suite","v6"
        #if dataset_items['domain'][:6] != dataset_items['primer_suite'][:6]:
        #    sys.exit("ERROR: Domain ("+dataset_items['domain']+") -- Primer Suite ("+dataset_items['primer_suite']+") mismatch.")
        #if dataset_items['domain'][-2:].lower() != dataset_items['dna_region'].lower():
        #    sys.exit("ERROR: DNA Region ("+dataset_items['dna_region']+") -- Domain ("+dataset_items['domain']+") mismatch.")
        if item['dna_region'] not in item['primer_suite']:
            sys.exit("ERROR: DNA Region ("+item['dna_region']+") not found in Primer Suite ("+item['primer_suite']+")")
                
        
        
    def check_project_name(self,item):
        # CHECK: project name format: 3 parts; end with Bv6,Ev9,Av6 or something similar
        try:
            (a,b,c) = item['project'].split('_')
        except:
            sys.exit("ERROR: project not in correct format: "+item['project'])
        (a,b,c) = item['project'].split('_')
        #if c[0] not in [i[0].upper() for i in domains]:
        #    sys.exit("ERROR : Project suffix has incorrect/non-existant domain: "+c)
        if c[1:] not in self.dna_regions:
            sys.exit("ERROR : Project suffix has incorrect DNA region: "+c)
            
            
    def check_headers(self, headers):
        if sorted(self.known_header_list) != sorted(headers):
            sys.exit("ERROR : unknown_headers:\nyours: "+ ' '.join(sorted(headers))+"\nours:  "+' '.join(sorted(self.known_header_list)))
        else:
            return True
            
    def validate_454_ini():
        print "Validating ini type Config File -TODO"
        # must be a general section
        # should I create a dict here??? -That would render much code in
        #    runconfig useless.
        # are we going to continue developing ini style config files if
        #   no one uses them?
        
        import ConfigParser
    
        configDict = {}
        user_config = ConfigParser.ConfigParser()
        try:
            user_config.read(config_info)
        except:
            sys.exit("Failed to read Config file: Are you sure it is in the correct .ini format?")
        for section in user_config.sections():
    
            section_dict = configDict[section] = {}
            for option in user_config.options(section):
                section_dict[option] = user_config.get(section,option)
        print "Finished Validating"
        return configDict
        
        
#     def validate_csv2(self):
#         print "Validating csv type Config File"
#         
#         # changes spaces to '_' and all lowercase
#         known_header_list = ["run_key","lane","dataset","project","tubelabel","barcode","adaptor","dna_region",
#                                 "amp_operator","seq_operator","barcode_index","overlap","insert_size","file_prefix","read_length","primer_suite" ]
#         primer_suites = ["bacterialv6suite","archaealv6suite","eukaryalv9suite"]
#         dna_regions = ["v1","v3","v4","v5","v6","v9"]
#         
#         
#         data = {}
#         projects = {}
#         megadata = {}
#         megadata['general'] = {}
#         test_datasets = {}
#         dataset_counter = {}
#         headers = ''
#         f_in_md = open(self.config_info, 'r')
#         # must be comma sep
#         lines = f_in_md.readlines()
#         #headerLine = lines.pop(0).strip() #removes and returns the first line
#         #headers = [i.strip('"').lower().replace(" ", "_") for i in headerLine.split(',')]
#         
#         
#     
#         #lines.pop(0)
#         
#         
#         #[general] section
#         #run_date = 20120601
#         #platform = illumina
#         #input_dir = /xraid2-2/sequencing/Illumina/20120525_recalled/Project_Sandra_v6/analysis/
#         #output_dir = .
#         #input_file_names = 
#         #input_file_suffix = fa.uniques
# 
#         #input_file_formats = fasta
#         #input_file_lanes = 1 
#         found_data_lines = False
#         
#         
#         
#         for line in lines:
#             line = line.strip()
#             
#             if not line:
#                 continue
#             lst = [i.strip('"').replace(" ", "_") for i in line.strip().split(',')]
#             
#             
#             if not lst[0]:
#                 continue
#             #There can be variable numbers of lines before the data (##DATA##) line
#             if not found_data_lines:
#                 # need to get required general header items: run, platform, vamps_user
#                 if lst[0] == '##DATA##':
#                     found_data_lines = True
#                 else:
#                     try:
#                         megadata['general'][lst[0]] = lst[1]
#                     except:
#                         logger.error("The general item "+lst[0]+" has no value - Continuing without")
#                         megadata['general'][lst[0]] = ''
#                     if lst[0] == 'run':
#                         megadata['general']['run_date'] = lst[1]
#                 
#                 
#             elif found_data_lines == True:
#                 
#                 temp = {}   
#                 if not headers:
#                     headers = [i.strip('"').lower().replace(" ", "_") for i in line.split(',')]
#         
#                     if sorted(known_header_list) != sorted(headers):
#                         sys.exit("ERROR : unknown_headers:\nyours: "+ ' '.join(sorted(headers))+"\nours: "+' '.join(sorted(known_header_list)))
#                 else:
#                     for n in range(0,len(headers)):
#                         #print headers[n], lst[n]
#                         temp[headers[n]] = lst[n]
# 
#                 
#                     temp['file_prefix'] = temp['dataset']+'_'+temp['barcode'].upper().replace('N','')
#                 
#                     #data[lst[0]] = temp
#                     
#                     idx_run_key = temp['barcode_index']+'_'+temp['run_key']
#                     
#                     
#                     #unique = str(random.randrange(1000000, 9999999))
#                     #submit_code=temp['vamps_user']+'_'+unique
#                     megadata[idx_run_key]={}
#                     
#                     if idx_run_key in test_datasets:
#                         sys.exit("ERROR: duplicate index:run_key: "+idx_run_key+" - Exiting")
#                     else:                     
#                         test_datasets[idx_run_key] = 1
#                         
#                     megadata[idx_run_key]['dataset'] = temp['dataset']
#                     #megadata[idx_run_key]['submit_code'] = submit_code
#                     megadata[idx_run_key]['project'] = temp['project']
#                     
#                     if temp['project'] in dataset_counter:
#                         dataset_counter[temp['project']] += 1
#                     else:
#                         dataset_counter[temp['project']] = 1
#                     
#                     #megadata[idx_run_key]['ds_count'] = 1
#                     megadata[idx_run_key]['project'] = temp['project']
#                     megadata[idx_run_key]['run_key'] = temp['run_key']
#                     megadata[idx_run_key]['lane'] = temp['lane']
#                     megadata[idx_run_key]['tubelabel'] = temp['tubelabel']
#                     megadata[idx_run_key]['barcode'] = temp['barcode']
#                     megadata[idx_run_key]['adaptor'] = temp['adaptor']
#                     megadata[idx_run_key]['dna_region'] = temp['dna_region']
#                     megadata[idx_run_key]['amp_operator'] = temp['amp_operator']
#                     megadata[idx_run_key]['seq_operator'] = temp['seq_operator']
#                     megadata[idx_run_key]['barcode_index'] = temp['barcode_index']
#                     megadata[idx_run_key]['overlap'] = temp['overlap']
#                     megadata[idx_run_key]['insert_size'] = temp['insert_size']
#                     megadata[idx_run_key]['file_prefix'] = temp['file_prefix']
#                     megadata[idx_run_key]['read_length'] = temp['read_length']
#                     megadata[idx_run_key]['primer_suite'] = temp['primer_suite']
#                 
#             else:
#                 sys.exit("ERROR: Manditory ##DATA## line not found in Config File - Exiting"+line)
#         if not found_data_lines:
#             sys.exit("No data found: csv config file must have a ##DATA## line to mark the beginning of the data\n Are you sure this is a csv file?")
#         
#         print 'general:',megadata['general']
#         
#         
#         
#         # start error checking here
#         # MUST be in list: "Domain","Primer Suite","DNA Region"
#         # MUST MATCH: "Domain","Primer Suite","DNA Region"
#         #
#         #
#         # VAMPS project name format:  SLM_GCB_Bv6
#         #
#         #
#         if 'file_suffix' not in megadata['general']:
#             megadata['general']['file_suffix'] = ''
#         file_count = 0
#         files_list = []
#         if os.path.isdir(megadata['general']['input_dir']):
#             p = megadata['general']['input_dir'], '*'+megadata['general']['file_suffix']
#             print p
#             for infile in glob.glob( os.path.join(megadata['general']['input_dir'], '*'+megadata['general']['file_suffix']) ):
#                 files_list.append(os.path.basename(infile))
#                 file_count += 1
#         else:
#             sys.exit("no input directory")
#             
#         if not file_count:
#             sys.exit("No files were found in '"+megadata['general']['input_dir']+"' with a suffix of '"+megadata['general']['file_suffix']+"'")
#             
#         megadata['general']['input_file_names'] = ','.join(files_list)
#         megadata['general']['input_file_formats'] = ','.join([megadata['general']['input_file_format'] for i in files_list])
#         # assign 1 to lanes -- kludge
#         megadata['general']['input_file_lanes'] = ','.join(['1']*file_count)
#         #print dataset_counter
#         for item in megadata:
#             if item != 'general':
#                 
#             
#             #for dataset_items in megadata[idx_run_key]['datasets']:
#                 #dataset_items['domain']        = dataset_items['domain'].lower().replace(" ", "_")
#                 megadata[item]['primer_suite']  = megadata[item]['primer_suite'].lower().replace(" ", "_")
#                 megadata[item]['dna_region']    = megadata[item]['dna_region'].lower().replace(" ", "_")
#                 megadata[item]['barcode']        = megadata[item]['barcode'].upper()
#                 megadata[item]['barcode_index']  = megadata[item]['barcode_index'].upper()
#                 megadata[item]['ds_count'] = str(dataset_counter[megadata[item]['project']])
#                 #print    dataset_counter,megadata[item]['project']
#                 
#                 #print project,dataset_items,"\n\n"
#     
#                 if not megadata[item]['dataset']:
#                     sys.exit("ERROR:Current dataset name is missing or corrupt - Exiting")
#                 
#                 for k,v in megadata[item].iteritems():
#                     if not k:
#                         sys.exit("ERROR: key for: '"+v+"' is missing or corrupt - Exiting")
#                     if not v:
#                         sys.exit("ERROR: value of: '"+k+"' is missing or corrupt - Exiting")
#                 
#                 # CHECK MUST MATCH: "Domain","Primer Suite","DNA Region"
#                 if megadata[item]['primer_suite'] not in primer_suites:
#                     sys.exit("ERROR: Primer Suite not found: "+megadata[item]['primer_suite'])
#                 #if dataset_items['domain'] not in domains:
#                  #   sys.exit("ERROR: Domain not found: "+dataset_items['domain'])
#                 if megadata[item]['dna_region'] not in dna_regions:
#                     sys.exit("ERROR: DNA Region not found: "+megadata[item]['dna_region'])
#                 # "Bacterial v6","BacterialV6Suite","v6"
#                 #if dataset_items['domain'][:6] != dataset_items['primer_suite'][:6]:
#                 #    sys.exit("ERROR: Domain ("+dataset_items['domain']+") -- Primer Suite ("+dataset_items['primer_suite']+") mismatch.")
#                 #if dataset_items['domain'][-2:].lower() != dataset_items['dna_region'].lower():
#                 #    sys.exit("ERROR: DNA Region ("+dataset_items['dna_region']+") -- Domain ("+dataset_items['domain']+") mismatch.")
#                 if megadata[item]['dna_region'] not in megadata[item]['primer_suite']:
#                     sys.exit("ERROR: DNA Region ("+megadata[item]['dna_region']+") not found in Primer Suite ("+megadata[item]['primer_suite']+")")
#                 
#             
#                 # CHECK: project name format: 3 parts; end with Bv6,Ev9,Av6 or something similar
#                 try:
#                     (a,b,c) = megadata[item]['project'].split('_')
#                 except:
#                     sys.exit("ERROR project not in correct format: "+megadata[item]['project'])
#                 (a,b,c) = megadata[item]['project'].split('_')
#                 #if c[0] not in [i[0].upper() for i in domains]:
#                 #    sys.exit("ERROR : Project suffix has incorrect/non-existant domain: "+c)
#                 if c[1:] not in dna_regions:
#                     sys.exit("ERROR : Project suffix has incorrect DNA region: "+c)
#                    
#                     
#                 #print item,megadata[item],"\n\n"
#         # other checks to put in:
#         # check for duplicate dataset name:  NO
#         # that data == file prefix
#         # if we have an input directory that each dataset has a coresponding file - for illumina
#         # Missing data is ok for barcode and adaptor (illumina only)
#         # get tube number back:  yes ds_count
#         print "Finished Validating"
#         return megadata
    
# def send_metadata_to_database(data, data_object):
#     cursor = data_object['cursor']
#     cursor_env454 = data_object['cursor_env454']
#     sub_table = data_object['submission_tbl']
#     tubes_table_vamps = data_object['tubes_tbl']
#     seqs_table_env454 = data_object['seqs_454']
#     metadata_table_env454 = data_object['metadata_454']
#     vamps_user = data_object['vamps_user']
#     upload_code = data_object['upload_code']
#     datetime = data_object['datetime']
#     submit_code = upload_code+'_'+vamps_user
#     env_source = '100'  # 100 is unknown
#     platform = 'illumina' # this is ONLY for illumina data - right?
#     # get contact
#     cursor.execute("select last_name,first_name,email,institution from vamps_auth where user='"+vamps_user+"'")
#     (last_name,first_name,email,institution) = cursor.fetchone()
#     for project in data:
#         
#         insert_string = "insert ignore into %s (\
#             upload_code,\
#             temp_project,\
#             user,\
#             last_name,\
#             first_name,\
#             email,\
#             institution,\
#             env_source_id,\
#             num_of_tubes,\
#             date_initial\
#             ) \
#             values('"+upload_code+"',\
#                 '"+project+"',\
#                 '"+vamps_user+"',\
#                 '"+last_name+"',\
#                 '"+first_name+"',\
#                 '"+email+"',\
#                 '"+institution+"',\
#                 '"+env_source+"',\
#                 '"+str(data[project]['ds_count'])+"',\
#                 '"+datetime+"'\
#                 )"
#         print insert_string
#         cursor.execute(insert_string % (sub_table))
#         ds_num = 0
#         #known_header_list = ["run_key","run","lane","dataset","project",
#         #"tubelabel","barcode","adaptor","dna_region","amp_operator",
#         #"seq_operator","barcode_index","overlap","insert_size",
#         #"file_prefix","read_length","primer_suite"]
#         for dataset_items in data[project]['datasets']:
#             ds_num += 1
#             insert_string = "insert into %s (\
#                 upload_code,\
#                 runkey,\
#                 run,\
#                 lane,\
#                 project,\
#                 dataset,\
#                 tubelabel,\
#                 barcode,\
#                 adaptor,\
#                 dna_region,\
#                 amp_operator,\
#                 seq_operator,\
#                 barcode_index,\
#                 overlap,\
#                 insert_size,\
#                 file_prefix,\
#                 read_length,\
#                 primer_suite,\
#                 date_initial\
#                 ) \
#                 values('"+upload_code+"',\
#                     '"+dataset_items['run_key']+"',\
#                     '"+dataset_items['run']+"',\
#                     '"+dataset_items['lane']+"',\
#                     '"+project+"',\
#                     '"+dataset_items['dataset']+"',\
#                     '"+dataset_items['tubelabel']+"',\
#                     '"+dataset_items['barcode']+"',\
#                     '"+dataset_items['adaptor']+"',\
#                     '"+dataset_items['dna_region']+"',\
#                     '"+dataset_items['amp_operator']+"',\
#                     '"+dataset_items['seq_operator']+"',\
#                     '"+dataset_items['barcode_index']+"',\
#                     '"+dataset_items['overlap']+"',\
#                     '"+dataset_items['insert_size']+"',\
#                     '"+dataset_items['file_prefix']+"',\
#                     '"+dataset_items['read_length']+"',\
#                     '"+dataset_items['primer_suite']+"',\
#                     '"+datetime+"')"
#             
#             cursor.execute(insert_string % (tubes_table_vamps))
#             
#             cursor_env454.execute(insert_string % (metadata_table_env454))
    # ls /xraid2-2/sequencing/Illumina/20120525_recalled/Project_Sandra_v6/analysis/*fa.unique
    #barcode_index = 'idx': 'AAGCTA',
    #barcode = run_key = 'inline_barcode': 'NNNNACGCA
    
# def read_sequence_files(data, data_object):
#     # sample directory
#     seqs_table_env454 = data_object['seqs_454']
#     cursor_env454 = data_object['cursor_env454']
#     fasta_dir = "/xraid2-2/sequencing/Illumina/20120525_recalled/Project_Sandra_v6/analysis/"
#     file_suffix = "fa.unique"
#     file_count = 0
#     for filename in os.listdir(fasta_dir): 
#         
#         if filename[-len(file_suffix):] == file_suffix:
#             size = os.path.getsize(fasta_dir +filename)
#             
#             #cat 'filename' | grep ">" | wc -l
#             #print "result",x
#             print os.path.split(filename)[-1]
#             p = subprocess.Popen(["cat "+ fasta_dir +filename+ " | grep '>' | wc -l"], shell=True, stdout=subprocess.PIPE)
#             #p = subprocess.Popen(["cat", fasta_dir +filename, "|", "grep", '>',"|","wc","-l"], shell=True, stdout=subprocess.PIPE)
#             print "  seq_count:",p.communicate()[0].strip()
#             file_count += 1
#             
#             dataset = '_'.join(filename.split('-')[0].split('_')[:-1])
#             project = ''
#             for p in data:
#                 #print type(data[p]),data[p]['datasets'],"\n\n"
#                 data[p]['datasets'][0]['dataset']
#                 for n in range(0,len(data[p]['datasets'])):
#                     if data[p]['datasets'][n]['dataset'] == dataset:
#                         project = p
#                         #domain = data[p]['datasets'][n]['domain']
#                         dna_region = data[p]['datasets'][n]['dna_region']
#                         break
#             f = FastaReader(fasta_dir +filename)
#             while f.next():
#                 id = f.id.split('|')[0]
#                 #print id
#                 cursor_env454.execute("insert ignore into "+seqs_table_env454+" (read_id,sequence,length,project,dataset,source) \
#                     VALUES('"+id+"','"+f.seq+"','"+str(len(f.seq))+"','"+project+"','"+dataset+"','"+dna_region+"')")
#     print "File Count",file_count
        
if __name__ == '__main__':

    m =  Metadata_utils
    
    # import argparse
#     
#     # DEFAULTS
#     user = ''  
#     #project = 'p'+str(random.randrange(100000,999999))
#    
#     unique = str(random.randrange(1000000, 9999999))
#     data_object = {}
#     
# 
#     
#     myusage = """usage: illumina_input.py -m metadatafile -d seqdir [options]
#          
#          
#          
#          where
#             -m, --metadatafile The name of the text file.  [required]
#             
#             -d, --directory    The name of the directory where the sequence files are located.   [required]
#             
#               
#             -site            vamps or vampsdev.
#                                 [default: vampsdev]
#             -r,  --run
#             -u, --vamps_user       Needed for database
#                                 vamps_user is to keep a record of who uploaded
#             
#     
#     
#     """
#     parser = argparse.ArgumentParser(description="Read Illumina directory, scan/validate metadata file and import sequences and metadata" ,usage=myusage)
#     
#     parser.add_argument("-u", "--vamps_user",         required=True,  action="store",   dest = "vamps_user", 
#                                                     help="user name")  
#                                          
#     parser.add_argument("-m", "--metadatafile",  required=True,  action="store",   dest = "metadata_file", 
#                                                     help="Metadata File ") 
#     parser.add_argument("-d", "--directory",     required=False,  action="store",   dest = "seqsdir", 
#                                                     help="Taxonomy File ")
#     parser.add_argument("-s", "--site",     required=False,  action="store",   dest = "site", default='vampsdev',
#                                                     help="Taxonomy File ")                                      
#     parser.add_argument("-r", "--run",     required=False,  action="store",   dest = "runcode",
#                                                     help="Taxonomy File ")
#     print "Starting illumina_imput.py"
#     
#     args = parser.parse_args()
#     
#     data_object['metadata_file'] = args.metadata_file
#     data_object['seqsdir'] = args.seqsdir
#     #data_object['basedir'] = args.basedir
#     data_object['datetime'] = str(datetime.date.today())
#     data_object['vamps_user'] = args.vamps_user
#     if args.runcode:
#         data_object['upload_code'] = args.vamps_user+'_'+args.runcode
#     else:
#         data_object['upload_code'] = args.vamps_user+'_'+unique
#     
#     
#     data_object['submission_tbl'] = 'vamps_submissions_illumina'
#     data_object['tubes_tbl'] = 'vamps_submissions_tubes_illumina'
#     data_object['metadata_454'] = 'illumina_metadata_av'    
#     data_object['seqs_454'] = 'illumina_seqs_av'
#     
#     if args.site:
#         site = args.site
#     
#     if site == 'vamps':
#         db_host = 'vampsdb'
#         db_name = 'vamps'
#         db_home = '/xraid2-2/vampsweb/vamps/'
#     else:
#         db_host = 'vampsdev'
#         db_name = 'vamps'
#         db_home = '/xraid2-2/vampsweb/vampsdev/'
#     
#     db_host_env454 = 'newbpcdb2'
#     db_name_env454 = 'env454'
#     db_home_env454 = os.getenv("HOME")
#     print db_home_env454
#     obj_env454=ConMySQL.New(db_host_env454, db_name_env454, db_home_env454)
#     data_object['cursor_env454'] = obj_env454.get_cursor()
#     obj=ConMySQL.New(db_host, db_name, db_home)
#     data_object['cursor'] = obj.get_cursor()
#     
#     dbuser = obj.get_db_user()
#     data_object['db_user'] = dbuser
#     
#     data = validate_metadata_file(data_object)
#     send_metadata_to_database(data, data_object)
#     read_sequence_files(data, data_object)
#     
#     data_object['cursor'].close()
    

