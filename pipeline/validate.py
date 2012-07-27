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
from pipeline.db_upload import MyConnection
#import pipeline.fastalib


  
    
class MetadataUtils:
    """
    Class to read metadata files (csv and ini style)
    validate and create a dictionary from them
    """
    Name = "MetadataUtils"
    def __init__(self, args):
        self.args = args
        self.known_header_list_illumina = C.csv_header_list_illumina
        self.known_header_list_454 = C.csv_header_list_454
        self.primer_suites     = C.primer_suites 
        self.dna_regions       = C.dna_regions
        self.data_object = {}
        self.data_object['general'] = {}
        self.warn_msg = """\n\tThe config File () seems to be okay. If the items above look correct
        then press 'c' to continue the pipeline\n"""
            
            
    def validate(self):  
        if self.args.platform == 'illumina':
            self.warn_msg = self.validate_illumina_ini()
        elif self.args.platform == '454':
            data = self.validate_454_ini()
        elif self.args.platform == 'ion_torrent':
            pass
        else:
            sys.exit("Unknown platform and configFile type for validation")
            
    def convert_and_save_ini(self):
        # converts csv to ini and saves to output_dir
        if self.args.config_file_type == 'csv':
            ini_file = self.convert_csv_to_ini()
            self.args.configPath  = ini_file
            self.args.config_file_type ='ini'
        elif self.args.config_file_type == 'ini':
            ini_file = self.save_ini_file()
            self.args.configPath  = ini_file
            self.args.config_file_type ='ini'
        else:
            sys.exit("Unknown config file type: "+config_file_type) 
    
    def get_general_data(self):
        """
        """
        return self.data_object['general']
        
    def create_dictionary_from_ini(self):
        """
        # read an ini config file and convert to a dictionary
        """
        import ConfigParser
        if os.path.exists(self.args.configPath):
            data_object = {}
            user_config = ConfigParser.ConfigParser()
            user_config.read(self.args.configPath)
            
            for section in user_config.sections():
                
                section_dict = data_object[section] = {}
                for option in user_config.options(section):
                    section_dict[option] = user_config.get(section,option)
                    
        else:
            print "error could not open config file: ",self.args.configPath
        
        return data_object 

    def get_command_line_items(self, general_data):
    
        # command line items take precedence over ini file items of the same name
        # defaults should be here and NOT in argparse/commandline
        if self.args.input_dir:       
            general_data['input_dir'] = self.args.input_dir
        else:
            if not general_data['input_dir']:
                general_data['input_dir'] = './'
        
        if self.args.run:
            general_data['run'] = self.args.run
            general_data['run_date'] = self.args.run
        else:
            if 'run' in general_data:                
                general_data['run_date'] = general_data['run']
            elif 'run_date' in general_data:
                general_data['run'] = general_data['run_date']
            else:
                sys.exit("Cannot find the run or run_date on command line or in config file - Exiting")
        # make sure RUN is before OUTPUT_DIR        
        try:
            general_data['output_dir'] = os.path.join(self.args.baseoutputdir,self.args.run)
        except:
            if 'output_dir' not in general_data:
                general_data['output_dir'] = os.path.join('.',self.args.run)       
        #getattr(args,'force_runkey', "")
        
        
        if self.args.platform:
            general_data['platform'] = self.args.platform
        else:
            if 'platform' not in general_data:
                sys.exit("Cannot find the platform from command line or in config file - Exiting")
                
        
        if self.args.input_file_format:
            general_data['input_file_format'] = self.args.input_file_format
        else:
            if 'input_file_format' not in general_data:
                general_data['input_file_format'] = ''
        if self.args.input_file_suffix:
            general_data['input_file_suffix'] = self.args.input_file_suffix
        else:
            if 'input_file_suffix' not in general_data:
                general_data['input_file_suffix'] = ''
        
        return general_data
        
#     def validate_454_csv(self, args, my_csv):
#         print "TODO: write validate def for 454/csv"
#         data_object = self.populate_data_object_454(args, my_csv)
        
        
    def validate_454_ini(self):
        print "Validating ini type Config File"
        print "TODO - write validation def for 454/ini"
        self.data_object = self.create_dictionary_from_ini() 
        # 454 ini file requirements:
        
        
        
    def validate_illumina_ini(self):
        """
        The csv headers are checked earlier
        """
        
        print "Validating ini type Config File (may have been converted fron csv)"
        return_code = False
        error_code  = False
        warn_code   = False
        msg = ''
        error=False
        warn=False
        self.data_object = self.create_dictionary_from_ini()
        
        
        (error_code,warn_code) = self.check_for_missing_values(self.data_object)  
        if error_code: error=True
        if warn_code: warn=True
        (error_code,warn_code) = self.check_for_datasets(self.data_object)
        if error_code: error=True
        if warn_code: warn=True
        (error_code,warn_code) = self.check_domain_suite_region(self.data_object)
        if error_code: error=True
        if warn_code: warn=True
        (error_code,warn_code) = self.check_project_name(self.data_object)
        if error_code: error=True
        if warn_code: warn=True
        (error_code,warn_code) = self.check_projects_and_datasets(self.data_object)
        if error_code: error=True
        if warn_code: warn=True
        #print self.data_object['input_dir']
        #print self.data_object['input_files']
        if 'input_dir' in self.data_object['general'] and self.data_object['general']['input_dir'] and ('input_files' not in self.data_object['general'] or not self.data_object['general']['input_files']):
            logger.error("There are no files found in the input directory: "+self.data_object['general']['input_dir'])
            error=True
        elif 'input_dir' not in self.data_object['general'] and 'input_files' not in self.data_object['general']:
            logger.warning("No input directory and no input files")        
            warn=True
        
        if error:
            sys.exit( """\n\tTHERE WERE SEVERE PROBLEMS WITH THE CONFIG FILE - EXITING 
            PLEASE CORRECT THEM AND START OVER.\n
            To view the errors add ' --loglevel info' to the command line.\n""")
        elif warn: 
            msg = """\n\tTHERE WERE NON-FATAL PROBLEMS WITH THE CONFIG FILE THAT MAY OR MAY NOT CAUSE PROBLEMS.\n
                To view the warnings add ' --loglevel warning' to the command line.\n"""
        
        return msg
        
    def validate_dictionary(self, config_info):
        """
        This is only used for data that comes in as a dictionary rather than a file
        such as with vamps user uploads
        """
        print "TODO - Validating input dictionary"
        # must be a general section
        # should I create a dict here??? -That would render much code in
        #    runconfig useless.
        # are we going to continue developing ini style config files if
        #   no one uses them?  
        configDict = config_info

        return configDict   
        

















        
    def populate_data_object_454(self, args):
        data = {}
        data['general'] = {}
        test_datasets = {}
        dataset_counter = {}
        headers = ''
        if self.run:
            infile = self.run.configPath
        else:            
            infile = args.configPath
            data['general']['input_dir'] = args.input_dir
            data['general']['output_dir'] = os.path.join(args.baseoutputdir,args.run)
            data['general']['platform'] = args.platform
            data['general']['run'] = args.run
            data['general']['run_date'] = args.run
            data['general']["input_file_format"] = args.input_file_format
            data['general']["input_file_suffix"] = args.input_file_suffix
    
        return data['general']

    
    def populate_data_object_illumina(self, args, my_csv):
        data = {}
        data['general'] = {}
        test_datasets = {}
        dataset_counter = {}
        headers = ''
        if self.run:
            infile = self.run.configPath
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
            data['general']['output_dir'] = os.path.join(args.baseoutputdir,args.run)
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
        projects = {}
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
    
        
    def get_input_files(self,file_suffix):
    
        files_list = []
        if os.path.isdir(self.args.input_dir):
            
            for infile in glob.glob( os.path.join(self.args.input_dir, '*'+file_suffix) ):
                if os.path.isdir(infile) == True:
                    pass
                else:
                    files_list.append(os.path.basename(infile))
        else:
            logger.warning("No input directory or directory permissions problem: "+self.args.input_dir)
            
        return files_list
        
    def check_for_input_files(self,data_object):
    
        file_count = 0
        files_list = []
        imports_list = []
        lanes_list = []


        #input_dir = os.path.join(data_object['general']['input_dir'],"fasta")
        input_dir = data_object['general']['input_dir']
        if os.path.isdir(input_dir):
            p = data_object['general']['input_dir'], '*'+data_object['general']['input_file_suffix']

            
            for infile in glob.glob( os.path.join(input_dir, '*'+data_object['general']['input_file_suffix']) ):
                files_list.append(os.path.basename(infile))
                for x in data_object:
                    if 'file_prefix' in data_object[x]:
                        pass
                        #print data_object[x]['file_prefix']
                        
                        #if os.path.basename(infile).split('-')[0] == data_object[x]['file_prefix']:
                            #lanes_list.append(data_object[x]['lane'])
                        
                file_count += 1
        else:

            logger.info("No input directory or directory permissions problem: "+input_dir)
            print "No input directory or directory permissions problem: "+input_dir
        if not file_count:
            #sys.exit("ERROR: No files were found in '"+input_dir+"' with a suffix of '"+data_object['general']['input_file_suffix']+"'")
            logger.info("ERROR: No files were found in '"+input_dir+"' with a suffix of '"+data_object['general']['input_file_suffix']+"'")

        data_object['general']['files_list'] = files_list
        data_object['general']['file_count'] = file_count
        # all the files in an illumina directory should be the same type
        #data_object['general']['file_formats_list'] = [data_object['general']["input_file_format"]] * file_count
        #data_object['general']['lanes_list'] = lanes_list
        #print "Files LIST",data_object['general']['files_list']
        
        
        return data_object
        
           
    def check_for_missing_values(self, data):
        missing_key   = ''
        error = False
        warn = False
        for item in data:
            if item == 'general':
                for k,v in data[item].iteritems():
                    if not k:
                        #sys.exit("ERROR: key for: '"+v+"' is missing or corrupt - Exiting")
                        logger.warning("(key: "+item+") key for: '"+v+"' is missing or corrupt - Continuing")
                        warn=True
                    if not v:                        
                        logger.warning("(key: "+item+") value of: '"+k+"' is missing or corrupt - Continuing")
                        warn=True
                            
        for item in data:
            if item != 'general':
                for k,v in data[item].iteritems():
                    if not k:
                        #sys.exit("ERROR: key for: '"+v+"' is missing or corrupt - Exiting")
                        logger.warning("(key: "+item+") key for: '"+v+"' is missing or corrupt - Continuing")
                        warn=True
                    if not v:
                        if (k == 'barcode' or k == 'adaptor'): #these could be empty
                            logger.warning("(key: "+item+") value of: '"+k+"' is missing or corrupt - Continuing")
                        else:
                            logger.error("(key: "+item+") value of: '"+k+"' is missing or corrupt - Continuing")
                            error=True
        return (error,warn)

    def check_for_datasets(self,data):
        error = False
        warn=False
        for item in data:
            if item != 'general':
                #print 'ds',data[item]['dataset']
                if not data[item]['dataset']:
                #if 'dataset' not in data[item]:
                    logger.error("Current dataset name is missing or corrupt - Exiting (key: "+item+")")
                    error=True
        return (error,warn) 
        
    def check_domain_suite_region(self,data):
        error = False
        warn=False
        
        for item in data:
            
            if item != 'general':
                # CHECK MUST MATCH: "Domain","Primer Suite","DNA Region"
                if data[item]['primer_suite'] not in self.primer_suites:
                    logger.error("Primer Suite not found: "+data[item]['primer_suite']+" - Exiting (key: "+item+")")
                    error=True
                #if dataset_items['domain'] not in domains:
                #   sys.exit("ERROR: Domain not found: "+dataset_items['domain'])
                if data[item]['dna_region'] not in self.dna_regions:
                    logger.error("DNA Region not found: "+data[item]['dna_region']+" - Exiting (key: "+item+")")
                    error=True
                # "Bacterial v6","BacterialV6Suite","v6"
                #if dataset_items['domain'][:6] != dataset_items['primer_suite'][:6]:
                #    sys.exit("ERROR: Domain ("+dataset_items['domain']+") -- Primer Suite ("+dataset_items['primer_suite']+") mismatch.")
                #if dataset_items['domain'][-2:].lower() != dataset_items['dna_region'].lower():
                #    sys.exit("ERROR: DNA Region ("+dataset_items['dna_region']+") -- Domain ("+dataset_items['domain']+") mismatch.")
                if data[item]['dna_region'] not in data[item]['primer_suite']:
                    logger.error("DNA Region ("+data[item]['dna_region']+") not found in Primer Suite ("+data[item]['primer_suite']+") - Exiting (key: "+item+")")
                    error=True
        return (error,warn)
        
    def check_project_name(self,data):
        """
        # CHECK: project name format: 3 parts; end with Bv6,Ev9,Av6 or something similar
        """
        error   =False
        warn    =False
        for item in data:
            if item != 'general':
                try:
                    (a,b,c) = data[item]['project'].split('_')
                except:
                    logger.error("project not in correct format: "+data[item]['project']+" - Exiting (key: "+data[item]+")")
                    error=True
                (a,b,c) = data[item]['project'].split('_')
                #if c[0] not in [i[0].upper() for i in domains]:
                #    sys.exit("ERROR : Project suffix has incorrect/non-existant domain: "+c)
                if c[1:] not in self.dna_regions:
                    logger.error("Project suffix has incorrect DNA region: "+c+" - Exiting (key: "+data[item]+")")
                    error = True
        return (error,warn)
            
    def check_projects_and_datasets(self,data):
        self.my_conn     = MyConnection(host='newbpcdb2', db="env454")  
        project_dataset = {}
        projects = {}
        datasets = {}
        error   =False
        warn    =False
        for item in data:
            if item != 'general':
                #project_dataset[data[item]['project']+'--'+data[item]['dataset']] = 1
                datasets[data[item]['dataset']] = data[item]['project']
                projects[data[item]['project']] = 1
        for p in projects:
            #print p 
            my_sql = """SELECT project FROM project WHERE project = '%s'""" % (p)
            res    = self.my_conn.execute_fetch_select(my_sql)
            if res:
                logger.warning("project '"+p+"' already exists in the database - is this okay?")
                warn = True
            else:
                logger.debug("project '"+p+"' is new")
                
            ds_found_count = 0   
            for d in datasets:
                if datasets[d] == p:
                    
                    #print "\t%s" % (d)
                    my_sql = """SELECT dataset FROM dataset WHERE dataset = '%s'""" % (d)
                    res    = self.my_conn.execute_fetch_select(my_sql)
                    if res:
                        ds_found_count += 1
                        if ds_found_count >3:
                            logger.warning("\t\tPossibly more .... - Exiting after just three")
                            break
                        logger.warning("\tdataset '"+d+"' already exists in the database - is this okay?")
                        warn=True
                    else:
                        logger.debug("\tdataset '"+d+"' is new")
            logger.debug("\tDataset Count: "+str(len(datasets)))
        return (error,warn)      
        
    def get_confirmation(self, steps, general_data):
        print "\n"
        for item,value in general_data.iteritems():
            print "%20s = %-20s" % (item,value)
        print "\nStep(s) to be performed: ",steps
        print "\n"+self.warn_msg+"\n"
        if 'validate' in steps.split(','):
            # print we are done
            sys.exit()
#        return raw_input("\nDoes this look okay? (q to quit, v to view configFile, c to continue) ")
        return "c"
        
    def convert_csv_to_ini(self):
        #print self.args
        from pipeline.get_ini import readCSV
        ini_file = os.path.join(self.args.baseoutputdir,self.args.run,self.args.run + '.ini')
        fh = open(ini_file,'w')
        my_csv = readCSV(file_path = self.args.configPath)
        
        content     = my_csv.read_csv()
        headers     = content[1].keys()
        headers_clean = [x.strip('"').replace(" ", "_").lower() for x in headers]
        projects = {}
        
        # get list of keys
        keys_list = []
        if self.check_headers(headers_clean):
        
            for k,values in content.iteritems():
                keys_list.append(values['barcode_index']+"_"+values['run_key']+"_"+values['lane'])
        
        # general section
        fh.write("[general]\n") 
        fh.write("run = "+self.args.run+"\n")
        fh.write("run_date = "+self.args.run+"\n")
        fh.write("config_file = "+os.path.join(self.args.baseoutputdir,self.args.run,self.args.run + '.ini')+"\n")
        fh.write("config_format = ini\n")
        fh.write("config_file_orig = "+self.args.configPath+"\n")
        fh.write("config_format_orig = "+self.args.config_file_type+"\n")
        fh.write("platform = "+self.args.platform+"\n")
        
        fh.write("output_dir = "+os.path.join(self.args.baseoutputdir,self.args.run)+"\n")
        fh.write("input_file_suffix = "  + getattr(self.args,'input_file_suffix', "")+"\n")
        fh.write("input_file_format = " + getattr(self.args,'input_file_format', "")+"\n")
        fh.write("anchor_file = "        + getattr(self.args,'anchor_file', "")+"\n")
        fh.write("primer_file = "        + getattr(self.args,'primer_file', "")+"\n")
        fh.write("require_distal = "     + getattr(self.args,'require_distal', "1")+"\n")
        fh.write("idx_keys = "           +','.join(keys_list)+"\n")
        fh.write("input_dir = "+getattr(self.args,'input_dir', ".")+"\n") 
        if self.args.input_dir:
            file_list = self.get_input_files(self.args.input_file_suffix)
            fh.write("input_files = "     + ','.join(file_list)+"\n") 
        else:
            fh.write("input_files = \n") 
        #fh.write(getattr(args,'force_runkey', ""))        
        
        for k,values in content.iteritems():
            fh.write("\n")
            if self.args.platform == 'illumina':
                fh.write("["+values['barcode_index']+"_"+values['run_key']+"_"+values['lane']+"]\n")
            elif self.args.platform == '454':
                fh.write("["+values['lane']+"_"+values['run_key']+"]\n")
                
            for v in values:
                fh.write(v+" = "+values[v]+"\n")
                
        fh.close()
        
        return ini_file 
        
    def save_ini_file(self):
        # give it a new name
        
        with open(os.path.abspath(self.args.configPath),"r") as fh_from:
            lines = fh_from.readlines()
        fh_from.close()
        ini_file = os.path.join(self.args.baseoutputdir,self.args.run,self.args.run + '.ini')
        with open(ini_file,'w') as fh_to:
            for line in lines:
                fh_to.write(line)              
        fh_to.close()
        
        return ini_file
            
    def check_headers(self,headers):
        if self.args.platform=='illumina':
            known_header_list= C.csv_header_list_illumina
        elif self.args.platform == '454':
            known_header_list = C.csv_header_list_454
        else:
            logger.error("in utils: check_headers - unknown platform")
        if sorted(known_header_list) != sorted(headers):
            sys.exit("ERROR : unknown_headers:\nyours: "+ ' '.join(sorted(headers))+"\nours:  "+' '.join(sorted(known_header_list)))
        else:
            return True
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
        
#if __name__ == '__main__':

#    m =  Metadata_utils
    
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
    

