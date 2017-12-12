# -*- coding: utf-8 -*-
import sys
import os
import traceback
#sys.path.append("/xraid/bioware/linux/seqinfo/bin")
#sys.path.append("/Users/ashipunova/bin/illumina-utils")
#sys.path.append("/Users/ashipunova/bin/illumina-utils/illumina-utils/scripts")
#sys.path.append("/bioware/merens-illumina-utils")
# sys.path.append("/bioware/pythonmodules/illumina-utils/")
import IlluminaUtils.lib.fastqlib as fq
#import fastqlib as fq
import IlluminaUtils.lib.fastalib as fa
#import fastalib as fa
from subprocess import call
import ast
from pipeline.utils import Dirs, PipelneUtils
from collections import defaultdict
import constants as C
import getpass
# import time

# from pipeline.pipelinelogging import logger

"TODO: add tests and test case"
#from collections import defaultdict

class IlluminaFiles:
    """
    0) from run create all dataset_lines names files in output dir
    1) split fastq files from casava into files with dataset_names
    2) create ini files 
    3) process them through Meren's script
    4) result - files dataset_lane-PERFECT_reads.fa.unique with frequencies - to process with env454upload()    
    
    """
    def __init__(self, runobj):
        self.utils = PipelneUtils()
        self.runobj         = runobj
        self.out_files      = {} 
        self.id_dataset_idx = {}
        self.in_file_path   = self.runobj.input_dir
                
        if self.runobj.vamps_user_upload:
            site = self.runobj.site
            dir_prefix=self.runobj.user+'_'+self.runobj.run
        else:
            site = ''
            dir_prefix = self.runobj.run
        if self.runobj.lane_name:
            lane_name = self.runobj.lane_name
        else:
            lane_name = ''
        
        dirs      = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, lane_name = lane_name, site = site) 
        self.dirs = dirs
        self.out_file_path = dirs.check_dir(dirs.analysis_dir)
        self.results_path  = dirs.check_dir(dirs.reads_overlap_dir)
        self.platform = self.runobj.platform
        
    def split_files(self, compressed = False):
        """
        TODO: *) fastq_file_names method to collect all file_names with full path or directories_names (see get_all_files()?)
        """   
#        print "compressed = %s" %       compressed
#        compressed = ast.literal_eval(compressed)     
        (in_files_r1, in_files_r2) = self.get_fastq_file_names(self.in_file_path)
#         correct_file_names = self.get_correct_file_names(in_files_r1)
        if (len(in_files_r1) > 0):
            self.read1(in_files_r1, compressed)
            self.read2(in_files_r2, compressed)
            self.create_inis()
        else:
#             print "ERROR: There is something wrong with fastq file names. Please check if they start with correct indexes."
#             logger.debug("ERROR: There is something wrong with fastq file names. Please check if they start with correct indexes.")
            self.utils.print_both("ERROR: There is something wrong with fastq file names. Please check if they start with correct indexes.")
        self.close_dataset_files()

            

#        self.perfect_reads()
#        self.uniq_fa()

    def open_dataset_files(self):
        file_name_base = [i + "_R1" for i in self.runobj.samples.keys()] + [i + "_R2" for i in self.runobj.samples.keys()]
        for f_name in file_name_base:
            output_file = os.path.join(self.out_file_path, f_name + ".fastq")
            self.out_files[f_name] = fq.FastQOutput(output_file)
        self.out_files["unknown"] = fq.FastQOutput(os.path.join(self.out_file_path, "unknown" + ".fastq"))        

    def close_dataset_files(self):
        [o_file[1].close() for o_file in self.out_files.iteritems()] 
        return
      
#     def perfect_reads(self):
#         self.utils.print_both("Extract perfect V6 reads:")
#         for idx_key in self.runobj.samples.keys():
#             file_name = os.path.join(self.out_file_path, idx_key + ".ini")
#             program_name = C.perfect_overlap_cmd
#             if self.utils.is_local():
#                 program_name = C.perfect_overlap_cmd_local       
#             try:
#                 if self.runobj.samples[idx_key].primer_suite.lower().startswith('archaeal'):
#                     call([program_name, file_name, "--archaea"]) 
#                 else: 
#                     call([program_name, file_name])
#             except:
#                 self.utils.print_both("Problems with program_name = %s, file_name = %s" % (program_name, file_name))
#                 raise  
#     
#     TODO: use from util
#     def call_sh_script(self, script_name_w_path, where_to_run):
#         try:
#             call(['chmod', '0774', script_name_w_path])
#             if self.utils.is_local():
#                 self.utils.print_both("call(['qsub', script_name_w_path], cwd=(where_to_run))")
#                 call(['bash', script_name_w_path], cwd=(where_to_run))                
#             else:
#                 call(['qsub', script_name_w_path], cwd=(where_to_run))
# #             pass
#         except:
#             self.utils.print_both("Problems with script_name = %s or qsub" % (script_name_w_path))
#             raise     
        
#     todo: combine and DRY with partial - it's the same command, different arguments
    def merge_perfect(self):
        self.utils.print_both("merge perfect V6 reads:")
        program_name = C.perfect_overlap_cmd
        if self.utils.is_local():
            program_name = C.perfect_overlap_cmd_local
        add_arg = " --marker-gene-stringent --retain-only-overlap --max-num-mismatches 0"
        command_line          = program_name + add_arg
        file_list             = self.dirs.get_all_files_by_ext(self.out_file_path, "ini")
        script_file_name      = self.create_job_array_script(command_line, self.dirs.analysis_dir, file_list)
        script_file_name_full = os.path.join(self.dirs.analysis_dir, script_file_name)
        self.utils.call_sh_script(script_file_name_full, self.dirs.analysis_dir)  
        return script_file_name    
    
    def trim_primers_perfect(self):
        self.utils.print_both("trim primers from perfect V6 reads:")
        
        merged_file_names = self.dirs.get_all_files_by_ext(self.dirs.reads_overlap_dir, "_MERGED")
        primer_suite = self.get_config_values('primer_suite')
        add_arg = ""
        if any([s.lower().startswith("Archaeal".lower()) for s in primer_suite]):
            add_arg += " --archaea"
        program_name = C.trim_primers_cmd + add_arg
        script_file_name      = self.create_job_array_script(program_name, self.dirs.reads_overlap_dir, merged_file_names)
        script_file_name_full = os.path.join(self.dirs.reads_overlap_dir, script_file_name)
        self.utils.call_sh_script(script_file_name_full, self.dirs.reads_overlap_dir)  
        return script_file_name    

    """    
    def perfect_reads_cluster(self):
        '''
        iu-merge-pairs anna.ini --marker-gene-stringent --retain-only-overlap --max-num-mismatches 0
​            Each flag is critical. ​marker-gene-stringent looks complete overlaps, retain-only-overlap gets rid of adapters, max-num-mismatches retains only perfect overlaps. 
            This generates the test_MERGED file with all complete overlaps without any mismatches. But it has all the primers. 
            Then we process this file with the new and shiny iu-analyze-v6-complete-overlaps script:
        iu-trim-V6-primers test_MERGED

        '''
        self.utils.print_both("Extract perfect V6 reads:")
        script_file_name      = self.merge_perfect()
        trim_script_file_name = self.trim_primers_perfect()

        return (script_file_name, trim_script_file_name)    
    """          
                              
    def partial_overlap_reads_cluster(self):
        self.utils.print_both("Extract partial_overlap reads (from partial_overlap_reads_cluster):")
        program_name = C.partial_overlap_cmd
        if self.utils.is_local():
            program_name = C.partial_overlap_cmd_local       
        dna_region = self.get_config_values('dna_region')
        if set(C.marker_gene_stringent_regions) & set(list(dna_region)):
            add_arg = "--marker-gene-stringent"
        else:
            add_arg = ""
#         TODO: this part is the same in perfect overlap - move into a method    
        command_line          = program_name + " --enforce-Q30-check " + add_arg 
        file_list             = self.dirs.get_all_files_by_ext(self.out_file_path, "ini")
        script_file_name      = self.create_job_array_script(command_line, self.dirs.analysis_dir, file_list)
        script_file_name_full = os.path.join(self.dirs.analysis_dir, script_file_name)
        self.utils.call_sh_script(script_file_name_full, self.dirs.analysis_dir)  
        self.utils.print_both("self.dirs.chmod_all(%s)" % (self.dirs.analysis_dir))
        self.dirs.chmod_all(self.dirs.analysis_dir)        
        
        return script_file_name      
                    
    def partial_overlap_reads(self):
        self.utils.print_both("Extract partial_overlap reads (from partial_overlap_reads):")
        for idx_key in self.runobj.samples.keys():
            ini_file_name = os.path.join(self.out_file_path, idx_key + ".ini")
            program_name = C.partial_overlap_cmd
            if self.utils.is_local():
                program_name = C.partial_overlap_cmd_local        
            try:
                if set(C.marker_gene_stringent_regions) & set(list(self.runobj.samples[idx_key].dna_region)):
                # if (self.runobj.samples[idx_key].dna_region == "ITS1"):
                    call([program_name, "--enforce-Q30-check", "--marker-gene-stringent", ini_file_name])
                else:
                    call([program_name, "--enforce-Q30-check", ini_file_name])
                               
#                 call([program_name, ini_file_name])           
#                 call([program_name, ini_file_name, idx_key])
#                 call([program_name, "--fast-merge", ini_file_name, idx_key])
            except Exception:
#                 except Exception, err:
                message = traceback.format_exc()
                self.utils.print_both(message)
    #or
#     print sys.exc_info()[0]

                self.utils.print_both("Problems with program_name = %s" % (program_name))
                raise  
                
#             print "HERE: program_name = " % (program_name)   
#             call([program_name, "--fast-merge", "--compute-qual-dicts", ini_file_name, idx_key])
            
    def get_config_values(self, key):
        config_path_data = [v for k, v in self.runobj.configPath.items()]
        return set([a[key] for a in config_path_data if key in a.keys()])
        
#     TODO: use from util
    def make_users_email(self):
        username = getpass.getuser() 
        return username + "@mbl.edu"
                
#     TODO: use from util
    def create_job_array_script(self, command_line, dir_to_run, files_list):
        files_string         = " ".join(files_list)
        files_list_size         = len(files_list)
        command_file_name = os.path.basename(command_line.split(" ")[0])
        script_file_name  = command_file_name + "_" + self.runobj.run + "_" + self.runobj.lane_name + ".sh"
        script_file_name_full = os.path.join(dir_to_run, script_file_name)
        log_file_name     = script_file_name + ".sge_script.sh.log"
        email_mbl         = self.make_users_email()
        text = (
                '''#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N %s
# Giving the name of the output log file
#$ -o %s
# Combining output/error messages into one file
#$ -j y
# Send mail to these users
#$ -M %s
# Send mail at job end (e); -m as sends abort, suspend.
#$ -m as
#$ -t 1-%s
# Now the script will iterate %s times.

  file_list=(%s)
  
  i=$(expr $SGE_TASK_ID - 1)
#   echo "i = $i"
  # . /etc/profile.d/modules.sh
  # . /xraid/bioware/bioware-loader.sh
  . /xraid/bioware/Modules/etc/profile.modules
  module load bioware
    
  echo "%s ${file_list[$i]}"  
  %s ${file_list[$i]}  
''' % (script_file_name, log_file_name, email_mbl, files_list_size, files_list_size, files_string, command_line, command_line)
# ''' % (script_file_name, log_file_name, email_mbl, files_list_size, files_list_size, files_string, command_line)
                )
        self.open_write_close(script_file_name_full, text)
        return script_file_name

    def filter_mismatches_cluster(self, max_mismatch = 3):
        self.utils.print_both("Filter mismatches if more then %s" % (max_mismatch))
        command_line = C.filter_mismatch_cmd
        if self.utils.is_local():
            command_line = C.filter_mismatch_cmd_local    
        files_dir = self.dirs.reads_overlap_dir   
                
        file_list             = self.dirs.get_all_files_by_ext(files_dir, "_MERGED")
        script_file_name      = self.create_job_array_script(command_line, files_dir, file_list)
        script_file_name_full = os.path.join(files_dir, script_file_name)
        self.utils.call_sh_script(script_file_name_full, files_dir)
        self.utils.print_both("self.dirs.chmod_all(%s)" % (files_dir))
        self.dirs.chmod_all(files_dir)        
        
        return script_file_name              

    def filter_mismatches(self, max_mismatch = 3):
        self.utils.print_both("Filter mismatches if more then %s" % (max_mismatch))
        n = 0        
        files = self.dirs.get_all_files()
        for full_name in files.keys():    
            if files[full_name][0].endswith('_MERGED'):
                n +=1   
#                print "%s fasta file: %s" % (n, full_name)
                program_name = C.filter_mismatch_cmd
                if self.utils.is_local():
                    program_name = C.filter_mismatch_cmd_local
#                 output_flag = "--output " + full_name + "_FILTERED"
# TODO:    Remove!!!
#                 output_flag = "-o " + full_name + "_FILTERED"           
#                 output_flag = "-o TTAGGC_NNNNTGACT_1_MERGED_FILTERED"           

#                 print "output_flag = %s" % (output_flag)
#                 print "%s %s %s" % (program_name, full_name, output_flag)                
#                 call([program_name, full_name, output_flag])
                call([program_name, full_name])

    def uniq_fa_cluster(self):
        self.utils.print_both("Uniqueing fasta files")
        command_line = C.fastaunique_cmd
        if self.utils.is_local():
            command_line = C.fastaunique_cmd_local   
        files_dir = self.dirs.reads_overlap_dir   
                
        file_list             = self.dirs.get_all_files_by_ext(files_dir, C.filtered_suffix)
        if len(file_list) == 0:
            file_list         = self.dirs.get_all_files_by_ext(files_dir, ".fa")
        if len(file_list) == 0:
            file_list         = self.dirs.get_all_files_by_ext(files_dir, "MERGED_V6_PRIMERS_REMOVED")
        
        script_file_name      = self.create_job_array_script(command_line, files_dir, file_list)
        script_file_name_full = os.path.join(files_dir, script_file_name)
        self.utils.call_sh_script(script_file_name_full, files_dir)  
        self.utils.print_both("self.dirs.chmod_all(%s)" % (files_dir))
        self.dirs.chmod_all(files_dir)        
        return script_file_name                           
                                       
    def uniq_fa(self):
        n = 0        
        self.utils.print_both("Uniqueing fasta files")
        files = self.dirs.get_all_files()
        for full_name in files.keys():    
#             if files[full_name][1] == ".fa" or files[full_name][0].endswith('_MERGED_FILTERED'):
            if files[full_name][1] == ".fa" or files[full_name][0].endswith(C.filtered_suffix):
                n +=1   
#                print "%s fasta file: %s" % (n, full_name)
                program_name = C.fastaunique_cmd
                if self.utils.is_local():
                    program_name = C.fastaunique_cmd_local                
                call([program_name, full_name])

    def get_primers(self):
        proximal_primer = ""
        distal_primer   = ""
        primers         = {}
        for idx_key in self.runobj.samples.keys():
            primer_suite = self.runobj.samples[idx_key].primer_suite.lower()
            
            # print "PPP primer_suite = "
            # print primer_suite

            if primer_suite in C.primers_dict:
                proximal_primer = C.primers_dict[primer_suite]["proximal_primer"]
                distal_primer = C.primers_dict[primer_suite]["distal_primer"]
                # print "RRR proximal_primer: %s. distal_primer: %s" % (proximal_primer, distal_primer)
            else:
                self.utils.print_both("ERROR! Something wrong with the primer suite name: %s. NB: For v6mod it suppose to be 'Archaeal V6mod Suite'\n" % (primer_suite))
            primers[idx_key] = (proximal_primer, distal_primer) 
            
        return primers
        
    def create_inis(self):
        for idx_key in self.runobj.samples.keys():
            run_key = idx_key.split('_')[1].replace("N", ".");
            "todo: check if works w/o NNNN when there is a proper csv"
            email = self.runobj.samples[idx_key].email
#        for dataset in self.dataset_emails.keys():
#            dataset_idx_base = dataset + "_" + self.dataset_index[dataset]
#            print "dataset = %s, self.dataset_emails[dataset] = %s" % (dataset, self.dataset_emails[dataset])
            text = """[general]
project_name = %s
researcher_email = %s
input_directory = %s
output_directory = %s

[files]
pair_1 = %s
pair_2 = %s
""" % (idx_key, email, self.out_file_path, self.results_path, idx_key + "_R1.fastq", idx_key + "_R2.fastq")

            "That's for parital overlap (v4v5 and hapto miseq illumina)" 
            if not self.runobj.do_perfect:
                primers = self.get_primers()    
                # print "run_key = %s, idx_key = %s, primers[idx_key][0], primers[idx_key][1] = %s" (run_key, idx_key, primers[idx_key][0], primers[idx_key][1])
                text += """
# following section is optional
[prefixes]
pair_1_prefix = ^""" + run_key + primers[idx_key][0] + "\npair_2_prefix = ^" + primers[idx_key][1]
                
            ini_file_name = os.path.join(self.out_file_path,  idx_key + ".ini")
            self.open_write_close(ini_file_name, text)

#     TODO: use from utils
    def open_write_close(self, script_file_name, text):
        ini_file = open(script_file_name, "w")
        ini_file.write(text)
        ini_file.close()
 
    def get_fastq_file_names(self, f_input_file_path):
        in_files_r1 = []
        in_files_r2 = []
        "TODO: exclude dir with new created files from the loop"
        for dirname, dirnames, filenames in os.walk(f_input_file_path):
            correct_file_names = self.get_correct_file_names(filenames)

            for filename in sorted(list(correct_file_names)):
                if filename.find('_R1_') > 0:
                    in_files_r1.append(os.path.join(dirname, filename))
                elif filename.find('_R2_') > 0:
                    in_files_r2.append(os.path.join(dirname, filename))
                else:
                    sys.stderr.write("No read number in the file name: %s\n" % filename)
        self.utils.print_both("FFF0: in_files_r1 %s\n, in_files_r2 %s" % (in_files_r1, in_files_r2))                    
        return (in_files_r1, in_files_r2)
    
    def get_correct_file_names(self, filenames):
        correct_file_names = [];
        for file1 in filenames:
            index_sequence = self.get_index(file1)
#             self.runobj.run_keys
#             
            good_run_key_lane_names = [x for x in self.runobj.run_keys if x.startswith(index_sequence)]
            if len(good_run_key_lane_names) > 0:
                correct_file_names.append(file1)
        return set(correct_file_names)
        
        
    def get_run_key(self, e_sequence, has_ns = "True"):
        if has_ns:
            return ("NNNN" + e_sequence[4:9])
        else:
            return e_sequence[0:5]
            
    def get_ini_run_key(self, index_sequence, e):
        has_ns = any("NNNN" in s for s in self.runobj.run_keys)           
        
        lane_number = e.lane_number
        if self.platform == "nextseq":
            lane_number = "1"
        return index_sequence + "_" + self.get_run_key(e.sequence, has_ns) + "_" + lane_number
        
    def read1(self, files_r1, compressed):
        """ loop through the fastq_file_names
            1) e.pair_no = 1, find run_key -> dataset name
            2) collect the relevant part of id
        """
        for file_r1 in files_r1:
            self.utils.print_both("====\nFFF1: file %s" % file_r1)
            f_input  = fq.FastQSource(file_r1, compressed)
            index_sequence = self.get_index(file_r1)
            while f_input.next(trim_to = C.trimming_length):
            # while f_input.next(trim_to = C.trimming_length[self.platform]):
                e = f_input.entry
                # todo: a fork with or without NNNN, add an argument
                #                 ini_run_key  = index_sequence + "_" + "NNNN" + e.sequence[4:9] + "_" + e.lane_number   
                # lane_number = e.lane_number
                # if self.platform == "nextseq":
                #     lane_number = "1"
                # ini_run_key  = index_sequence + "_" + self.get_run_key(e.sequence, has_ns) + "_" + lane_number
                ini_run_key = self.get_ini_run_key(index_sequence, e)
                if int(e.pair_no) == 1:
                    dataset_file_name_base_r1 = ini_run_key + "_R1"
                    if (dataset_file_name_base_r1 in self.out_files.keys()):
                        self.out_files[dataset_file_name_base_r1].store_entry(e)
                        "TODO: make a method:"
                        short_id1 = e.header_line.split()[0]
                        short_id2 = ":".join(e.header_line.split()[1].split(":")[1:])
                        id2 = short_id1 + " 2:" + short_id2
                        self.id_dataset_idx[id2] = ini_run_key
                    else:
                        self.out_files["unknown"].store_entry(e)
                    
    # def truncate_seq(self, seq):
    #     return seq[:C.trimming_length]
    
    def remove_end_ns_strip(self, e_sequence):
        if e_sequence.endswith('N'):
            return e_sequence.rstrip('N')
        else:
            return e_sequence
        
    def read2(self, files_r2, compressed):
        "3) e.pair_no = 2, find id from 2), assign dataset_name"
        for file_r2 in files_r2:
            self.utils.print_both("FFF2: file %s" % file_r2)
            f_input  = fq.FastQSource(file_r2, compressed)
            while f_input.next(trim_to = C.trimming_length):
                e = f_input.entry
                
#                 start = time.time()  
#                 time_before = self.utils.get_time_now()
#                 e.sequence = self.remove_end_ns_strip(e.sequence)
#                 elapsed = (time.time() - start)
#                 print "remove_end_ns_strip with strip is done in: %s" % (elapsed)      
                
                if (int(e.pair_no) == 2) and (e.header_line in self.id_dataset_idx):
                    file_name = self.id_dataset_idx[e.header_line] + "_R2"
                    self.out_files[file_name].store_entry(e)        
                else:
                    self.out_files["unknown"].store_entry(e)

    def get_index(self, file_r1):
        file_name_parts = os.path.basename(file_r1).split("_")
#        if the file name starts with "IDX, then actual idx will be next.
        index = file_name_parts[0]
        if file_name_parts[0].startswith("IDX"):
            index = file_name_parts[1]
        return index
