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
import datetime
from types import *
from pipeline.pipelinelogging import logger
import constants as C
from pipeline.db_upload import MyConnection
import re
#import pipeline.fastalib


  
    
class MetadataUtils:
    """
    Class to read metadata files (csv and ini style)
    validate and create a dictionary from them
    Two parts: 
    1) From pipeline-ui.py to validate the input args
    2) From runconfig.py to write the final ini file and create the dictionary
    that is used to create the run object
    """
    Name = "MetadataUtils"
    def __init__(self, command_line_args = None, configuration_dictionary = None):
        self.args = command_line_args
        self.general_config_dict = configuration_dictionary
        self.known_header_list  = C.csv_header_list
        self.pipeline_run_items = C.pipeline_run_items
        self.primer_suites      = C.primer_suites 
        self.dna_regions        = C.dna_regions
        self.data_object = {}
        self.data_object['general'] = {}
        self.warn_msg = """\n\tThe config File seems to be okay. If the items above look correct
        then press 'c' to continue the pipeline\n"""
        
        
            
    def convert_and_save_ini(self):
        
        new_ini_file = os.path.join(self.general_config_dict['output_dir'],self.general_config_dict['run'] + '.ini')
        #new_ini_file = os.path.join(self.general_config_dict['output_dir'],self.general_config_dict['run'],self.general_config_dict['run'] + '.ini')
        # converts csv to ini and saves to output_dir
        if self.general_config_dict['platform'] == 'vamps':
            self.save_ini_file(new_ini_file)
        else:
            self.convert_csv_to_ini(new_ini_file)
        self.general_config_dict['configPath']
        self.general_config_dict['configPath_original'] = self.general_config_dict['configPath']
        self.general_config_dict['configPath'] = new_ini_file
        
        # change path and type to new ini
        # regardless of what they were before        
    
    
    
    def validate(self): 
        
        if self.general_config_dict['platform'] == 'illumina':
            self.warn_msg = self.validate_illumina_ini()
        elif self.general_config_dict['platform'] == '454':
            data = self.validate_454_ini()
        elif self.general_config_dict['platform'] == 'ion_torrent':
            pass
        elif self.general_config_dict['platform'] == 'vamps':
            data = self.validate_vamps_ini()
        else:
            sys.exit("Unknown platform and configFile type for validation")
            

        return self.data_object
            
    def get_general_data(self):
        """
        """
        return self.data_object['general']
        
#     def create_dictionary_from_ini(self):
#         """
#         # read an ini config file and convert to a dictionary
#         """
#         import ConfigParser
#         if os.path.exists(self.general_config_dict['configPath']):
#             data_object = {}
#             user_config = ConfigParser.ConfigParser()
#             user_config.read(self.general_config_dict['configPath'])
#             
#             for section in user_config.sections():
#                 
#                 section_dict = data_object[section] = {}
#                 for option in user_config.options(section):
#                     section_dict[option] = user_config.get(section,option)
#                     
#         else:
#             print "error could not open config file: ",self.general_config_dict['configPath']
#         
#         return data_object 

#     def get_command_line_items(self, general_data):
#     
#         # command line items take precedence over ini file items of the same name
#         # defaults should be here and NOT in argparse/commandline
#         if self.args.input_dir:       
#             general_data['input_dir'] = self.args.input_dir
#         else:
#             if not general_data['input_dir']:
#                 general_data['input_dir'] = './'
#         
#         if self.args.run:
#             general_data['run'] = self.args.run
#             general_data['run_date'] = self.args.run
#         else:
#             if 'run' in general_data:                
#                 general_data['run_date'] = general_data['run']
#             elif 'run_date' in general_data:
#                 general_data['run'] = general_data['run_date']
#             else:
#                 sys.exit("Cannot find the run or run_date on command line or in config file - Exiting")
#         # make sure RUN is before OUTPUT_DIR        
#         try:
#             general_data['output_dir'] = os.path.join(self.args.baseoutputdir,self.args.run)
#         except:
#             if 'output_dir' not in general_data:
#                 general_data['output_dir'] = os.path.join('.',self.args.run)       
#         #getattr(args,'force_runkey', "")
#         
#         
#         if self.args.platform:
#             general_data['platform'] = self.args.platform
#         else:
#             if 'platform' not in general_data:
#                 sys.exit("Cannot find the platform from command line or in config file - Exiting")
#                 
#         
#         if self.args.input_file_format:
#             general_data['input_file_format'] = self.args.input_file_format
#         else:
#             if 'input_file_format' not in general_data:
#                 general_data['input_file_format'] = ''
#         if self.args.input_file_suffix:
#             general_data['input_file_suffix'] = self.args.input_file_suffix
#         else:
#             if 'input_file_suffix' not in general_data:
#                 general_data['input_file_suffix'] = ''
#         
#         return general_data
        
#     def validate_454_csv(self, args, my_csv):
#         print "TODO: write validate def for 454/csv"
#         data_object = self.populate_data_object_454(args, my_csv)
        
    def validate_vamps_ini(self):
        # configPath is the new configPath
        self.data_object = self.configDictionaryFromFile_ini(self.general_config_dict['configPath'])
        if 'fasta_file' in self.data_object and not os.path.exists(self.data_object['fasta_file']):
            sys.exit("Fasta file path doesn't exist: "+self.data_object['fasta_file'] )
        elif 'fasta_file' in self.data_object['general'] and not os.path.exists(self.data_object['general']['fasta_file']): 
            sys.exit("Fasta file path doesn't exist: "+self.data_object['general']['fasta_file'] )
                          
    def validate_454_ini(self):
        print "TODO - write validation def for 454/ini"
        #self.data_object = self.create_dictionary_from_ini() 
        # 454 ini file requirements:
        
        
        
    def validate_illumina_ini(self):
        """
        The csv headers are checked earlier
        """
        
        print "Validating ini type Config File (may have been converted from csv)"
        new_ini_file = os.path.join(self.general_config_dict['output_dir'],self.general_config_dict['run'] + '.ini')
        print "New ini file location: "+new_ini_file
        return_code = False
        error_code  = False
        warn_code   = False
        msg = ''
        error=False
        warn=False
        #print 'configpath',self.general_config_dict['configPath']
        # configPath here is the new configPath
        self.data_object = self.configDictionaryFromFile_ini(self.general_config_dict['configPath'])

        
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
        (error_code,warn_code) = self.check_dataset_name(self.data_object)
        if error_code: error=True
        if warn_code: warn=True
        (error_code,warn_code) = self.check_projects_and_datasets(self.data_object)
        if error_code: error=True
        if warn_code: warn=True
        #print self.data_object['input_dir']
        #print self.data_object['input_files']
 
 
        if 'input_dir' not in self.data_object['general'] and 'input_files' not in self.data_object['general']:
            logger.warning("No input directory and no input files")        
            warn=True
        elif not os.path.isdir(self.data_object['general']['input_dir']):
            logger.error("That is not a directory: "+self.data_object['general']['input_dir'])
            error=True
        elif self.data_object['general']['input_file_format'] == 'fastq' and self.data_object['general']['platform'] == 'illumina':
                file_exists = False
    #            if 'input_dir' in self.data_object['general'] and self.data_object['general']['input_dir']:
                for dirname, dirnames, filenames in os.walk(self.data_object['general']['input_dir']):
    #                if not filenames:
                    for file_name in filenames:
                        if os.path.isfile(os.path.join(dirname, file_name)):
                            file_exists = True
                            break
                if not file_exists:
                    logger.error("There are no files found in the input directory: "+self.data_object['general']['input_dir'])
                    error=True
        elif 'input_dir' in self.data_object['general'] and self.data_object['general']['input_dir'] and ('input_files' not in self.data_object['general'] or not self.data_object['general']['input_files']):
            logger.error("There are no files found in the input directory: "+self.data_object['general']['input_dir'])
            error=True
                        
        if error:
            sys.exit( """\n\t\033[91mTHERE WERE SEVERE PROBLEMS WITH THE CSV and/or CONFIG FILE - EXITING 
            PLEASE CORRECT THEM AND START OVER.\033[0m\n
            To view the errors add ' --loglevel info' to the command line.\n""")
        elif warn: 
            msg = """\n\t\033[93mTHERE WERE NON-FATAL PROBLEMS WITH THE CSV and/or CONFIG FILE THAT MAY OR MAY NOT CAUSE PROBLEMS.\033[0m\n
                To view the warnings add ' --loglevel warning' to the command line.\n"""
            print "\033[92mCSV File Passed Vaidation! (with warnings)\033[0m"
        else:
            print "\033[92mCSV File Passed Vaidation!\033[0m"
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
        if self.runobj:
            infile = self.runobj.configPath
        else:            
            infile = args.configPath
            data['general']['input_dir'] = args.input_dir
            #data['general']['output_dir'] = os.path.join(args.output_dir,args.run)
            data['general']['output_dir'] = args.output_dir
            data['general']['platform'] = args.platform
            data['general']['run'] = args.run
            #data['general']['run_date'] = args.run
            data['general']["input_file_format"] = args.input_file_format
            data['general']["input_file_suffix"] = args.input_file_suffix
    
        return data['general']

    

        
    def get_input_files(self):
        
        files_list = []
        
        if os.path.isdir(self.general_config_dict['input_dir']):
            
            for infile in glob.glob( os.path.join(self.general_config_dict['input_dir'], '*') ):
                if os.path.isdir(infile) == True:
                    
                    for infile2 in glob.glob( os.path.join( infile,'*') ):
                        if os.path.isdir(infile2) == True:
                            pass
                        else:
                            sub_dir = os.path.basename(infile)
                            
                            files_list.append(os.path.join(sub_dir,os.path.basename(infile2)))
                else:
                    files_list.append(os.path.basename(infile))
        else:
            if fasta_file:
                pass
            logger.warning("No input directory or directory permissions problem: "+self.general_config_dict['input_dir'])
            
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
                    if v == '':                        
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
                if data[item]['primer_suite'].lower() not in self.primer_suites:
                    logger.error("Primer Suite not found: "+data[item]['primer_suite']+" - Exiting (key: "+item+")")
                    error=True
                #if dataset_items['domain'] not in domains:
                #   sys.exit("ERROR: Domain not found: "+dataset_items['domain'])
                if data[item]['dna_region'].lower() not in self.dna_regions:
                    logger.error("DNA Region not found: "+data[item]['dna_region']+" - Exiting (key: "+item+")")
                    error=True
                # "Bacterial v6","BacterialV6Suite","v6"
                #if dataset_items['domain'][:6] != dataset_items['primer_suite'][:6]:
                #    sys.exit("ERROR: Domain ("+dataset_items['domain']+") -- Primer Suite ("+dataset_items['primer_suite']+") mismatch.")
                #if dataset_items['domain'][-2:].lower() != dataset_items['dna_region'].lower():
                #    sys.exit("ERROR: DNA Region ("+dataset_items['dna_region']+") -- Domain ("+dataset_items['domain']+") mismatch.")
                if data[item]['dna_region'].lower() not in data[item]['primer_suite'].lower():
                    logger.error("DNA Region ("+data[item]['dna_region']+") not found in Primer Suite ("+data[item]['primer_suite']+") - Exiting (key: "+item+")")
                    error=True
        return (error, warn)
        
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
        
    def check_dataset_name(self,data):
        """
        # CHECK: dataset name can be ONLY alphanumeric and underscore 
                    and cannot start with a number!
        """
        error   =False
        warn    =False
        for item in data:
            if item != 'general':
                dataset_name = data[item]['dataset']
                if not re.match("^[A-Za-z0-9_]*$", dataset_name):
                    logger.error("Dataset name has illeagal character(s): "+dataset_name+" (must be alphanumeric and underscore only)")
                    error = True
                #if  re.match("^[0-9]", dataset_name):
                 #   logger.error("Dataset name cannot begin with a digit: "+dataset_name)
                  #  error = True
                
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
            #print len(value)
            if type(value) != bool and len(value) > 80:
                tmp = value.split(',')
                print "%-20s = %s .. %s" % (item,tmp[0],tmp[-1])
            else:
                print "%-20s = %-20s" % (item,value)
        print "\nStep(s) to be performed: \033[1;36m",steps,'\033[0m'
        print "\n"+self.warn_msg+"\n"
        if 'validate' in steps.split(','):
            # print we are done
            sys.exit()
        print os.uname()
        print os.uname()[1]
        if os.uname()[1] == 'ashipunova.mbl.edu' or os.uname()[1] == 'as-macbook.local' or os.uname()[1] == 'as-macbook.home':
            return "c"
        else:
            return raw_input("\nDoes this look okay? (q to quit, v to view configFile, c to continue) ")
        
    def convert_csv_to_ini(self,new_ini_file):
        #print self.args
        from pipeline.get_ini import readCSV
        
        print 'CSV path',self.general_config_dict['csvPath']
        my_csv = readCSV(file_path = self.general_config_dict['csvPath'])
        
        content     = my_csv.read_csv()
        headers     = content[1].keys()
        headers_clean = [x.strip('"').replace(" ", "_").lower() for x in headers]
        projects = {}
        #print
        #print content[1]
        #print
        # get list of keys
        keys_list = []
        if self.check_headers(headers_clean):
            logger.info("CSV headers okay")
            for k,values in content.iteritems():
                keys_list.append(values['barcode_index']+"_"+values['run_key']+"_"+values['lane'])
        
        fh = open(new_ini_file,'w')
        # general section
        fh.write("#\n#\tCreated by MBL Pipeline for run: "+self.general_config_dict['run']+" on "+self.general_config_dict['date']+"\n#\n\n")  
        fh.write("[general]\n") 
        fh.write("run = "+self.general_config_dict['run']+"\n")
        fh.write("configPath = "+new_ini_file+"\n")
        
        fh.write("configPath_orig = "+self.general_config_dict['configPath']+"\n")
        fh.write("platform = "+self.general_config_dict['platform']+"\n")
        fh.write("output_dir = "          + self.general_config_dict['output_dir']+"\n")
        #fh.write("output_dir = "+os.path.join(self.general_config_dict['baseoutputdir'],self.general_config_dict['run'])+"\n")
        if self.general_config_dict['platform'] == 'illumina':
            #fh.write("input_file_suffix = "  + self.general_config_dict['input_file_suffix']+"\n")
            fh.write("input_file_format = " + self.general_config_dict['input_file_format']+"\n")
            fh.write("anchor_file = "        + self.general_config_dict['anchor_file']+"\n")
            fh.write("primer_file = "        + self.general_config_dict['primer_file']+"\n")
            fh.write("compressed = "          + str(self.general_config_dict['compressed'])+"\n")
            fh.write("database_host = "          + self.general_config_dict['database_host']+"\n")
            fh.write("database_name = "          + self.general_config_dict['database_name']+"\n")
            
        fh.write("input_dir = "          + self.general_config_dict['input_dir']+"\n")
        fh.write("require_distal = "     + str(self.general_config_dict['require_distal'])+"\n")
        fh.write("use_cluster = "              + str(self.general_config_dict['use_cluster'])+"\n")
        fh.write("date = "              + str(datetime.date.today())+"\n")
        fh.write("site = "              + self.general_config_dict['site']+"\n")
        fh.write("load_vamps_database = " + str(self.general_config_dict['load_vamps_database'])+"\n")
        fh.write("idx_keys = "           +','.join(keys_list)+"\n")
        if 'input_dir' in self.general_config_dict and self.general_config_dict['input_dir'] != '':
            file_list = self.get_input_files()
            fh.write("input_files = "     + ','.join(file_list)+"\n") 
        else:
            fh.write("input_files = \n") 
        #fh.write(getattr(args,'force_runkey', ""))        
 
 
        for k,values in content.iteritems():
            fh.write("\n")
            if self.general_config_dict['platform'] == 'illumina':
                fh.write("["+values['barcode_index']+"_"+values['run_key']+"_"+values['lane']+"]\n")
            elif self.general_config_dict['platform'] == '454':
                fh.write("["+values['lane']+"_"+values['run_key']+"]\n")
                
            for v in values:
                fh.write(v+" = "+values[v]+"\n")
                
        fh.close()
        
        return new_ini_file 
        
    def save_ini_file(self,new_ini_file):
        # give it a new name
        out_fh = open(new_ini_file,'w')
        #for line in open(os.path.abspath(self.general_config_dict['configPath']),"r"):
        #    out_fh.write(line)
        self.general_config_dict['configPath_original'] = self.general_config_dict['configPath']
        self.general_config_dict['configPath'] = new_ini_file
        
        out_fh.write("#\n#\tCreated by MBL Pipeline for run: "+self.general_config_dict['run']+" on "+self.general_config_dict['date']+"\n#\n\n")  
        out_fh.write("[general]\n")   
        for item in self.general_config_dict:
            
            out_fh.write(item+" = "+str(self.general_config_dict[item]) + "\n")
        #out_fh.write("\n["+self.general_config_dict['platform']+"]\n") 
        #for item in self.general_config_dict:
        #    if item not in C.general_run_items:
        #        out_fh.write(item+" = "+str(self.general_config_dict[item]) + "\n")
        
        
        
        if 'fasta_file' in self.general_config_dict and self.general_config_dict['fasta_file'] != '':
            (path,fasta) = os.path.split(self.general_config_dict['fasta_file'])
            if 'input_dir' in self.general_config_dict and self.general_config_dict['input_dir'] != path:
                sys.exit("Your input_dir and fasta_file directory don't agree - Exiting\n\t"+self.general_config_dict['input_dir']+" != "+self.general_config_dict['fasta_file'])
            
            out_fh.write("input_dir = "+path+"\n")
            out_fh.write("input_files = "+fasta+"\n")
            #out_fh.write("input_file_suffix = fasta\n")
        elif 'input_dir' in self.general_config_dict and self.general_config_dict['input_dir'] != '':
            file_list = self.get_input_files()
            out_fh.write("input_files = "     + ','.join(file_list)+"\n") 
        else:
            out_fh.write("input_files = \n") 
        out_fh.close()


            
    def check_headers(self,headers):
        if self.general_config_dict['platform']=='illumina':
            known_header_list= self.known_header_list['illumina']
        elif self.general_config_dict['platform'] == '454':
            known_header_list = self.known_header_list['454']
        else:
            logger.error("in utils: check_headers - unknown platform")
        #print   sorted(known_header_list)
        #print sorted(headers)
        if sorted(known_header_list) != sorted(headers):
            print "="*40
            print "csv file header problem"
            print "%-20s %-20s" % ("REQUIRED", "YOUR CSV")
            for i in sorted(known_header_list):
                if i in headers:
                    print "%-20s%-20s" % (i,i)
                else:
                    print "%-20s%-20s" % (i,"----------- <--- missing")
            for i in headers:
                
                if i not in known_header_list:
                    print "%-20s%-20s" % (" ",i+" <--- extra")
            print "="*40
            sys.exit("ERROR : unknown or missing headers\n")
        else:
            return True
        
    "TODO: add function to convert env_sample_source to env_sample_source_id"

    def configDictionaryFromFile_ini(self,config_file_path):
        import ConfigParser
        
        configDict = {}
        user_config = ConfigParser.ConfigParser()
        user_config.read(config_file_path)
        
        for section in user_config.sections():
            section_dict = configDict[section] = {}
            for option in user_config.options(section):
                section_dict[option] = user_config.get(section,option)
                if section_dict[option] == 'True' or section_dict[option] == 'true':
                    section_dict[option] = True
                elif section_dict[option] == 'False' or section_dict[option] == 'false':
                    section_dict[option] = False
                    
        return configDict
        
    def get_values(self, args, general_config_dict = {} ):
        collector={}

        for item in self.pipeline_run_items[args.platform]:
            
            # set collector[item] to the default first
            collector[item] = self.pipeline_run_items[args.platform][item]
            
            # now look for args (then ini) values to replace
            if item in args and getattr( args, item ) != None:
                collector[item]  = getattr( args, item )
            elif general_config_dict and item in general_config_dict[args.platform] and general_config_dict[args.platform][item] != '':
                collector[item]  = general_config_dict[args.platform][item]
        
        # get all the items from general_config_dict['general']
        if 'general' in general_config_dict:
            for item in general_config_dict['general']:
                collector[item]  = general_config_dict['general'][item]
            
               
        return collector
    
    def validate_args(self):
        """
        # THOUGHTS
        # vamps users
        # single project and dataset
        # Supply an ini file OR commandline (for web interface), but no csv file
        #
        # MBL pipeline
        # REQUIRE a csv file and a ini file
        """
        collector={}
        
        if self.args.configPath:
            general_config_dict = self.configDictionaryFromFile_ini(self.args.configPath) 
            if self.args.platform in general_config_dict and 'general' in general_config_dict:
                collector= self.get_values( self.args, general_config_dict)
            else:
                sys.exit("The ini file needs both a [general] and ["+ self.args.platform +"] section - Exiting.")
        else:
            # no configPath
            collector= self.get_values( self.args )
            
        if self.args.platform == 'illumina':
            print "Starting Illumina Pipeline"
            if not self.args.csvPath:
                sys.exit("illumina requires a csv file - Exiting")
            
        elif self.args.platform == 'vamps':
            print "Starting VAMPS Pipeline:"
            
            if 'project' not in collector or collector['project'] == '':    
                collector['project'] = collector['project'][:1].capitalize() + collector['project'][1:]
            else:
                logger.debug("No project found in vamps pipeline")
            if self.args.fasta_file:
                collector['project'] = self.args.fasta_file
                collector['from_fasta'] = True
        elif self.args.platform == '454':
            print "Starting 454 Pipeline"
            
        elif self.args.platform == 'ion_torrent':
            print "Starting Ion Torrent Pipeline"
            
        else:
            sys.exit("Validate args: Unknown Platform")
        
        if  self.args.configPath:
            collector['configPath'] = self.args.configPath
        else:
            collector['configPath'] = ""
        # these are all the bool items in the collector
        # they need to be converted fron str to bool here
        for i in collector:
            if collector[i] == 'True' or collector[i] == 'true':
                collector[i] = True
            elif collector[i] == 'False' or collector[i] == 'false':
                collector[i] = False
        
        #collector['runcode'] = self.args.run
        collector['run'] = self.args.run
        #collector['run_date'] = self.args.run
        #collector['steps'] = self.args.steps
        collector['platform'] = self.args.platform
        if self.args.input_dir:       
             collector['input_dir'] = self.args.input_dir

        collector['date'] = str(datetime.date.today())
        #print collector
        return collector
            