import sys
import os
import traceback
sys.path.append("/xraid/bioware/linux/seqinfo/bin")
sys.path.append("/Users/ashipunova/bin/illumina-utils")
sys.path.append("/Users/ashipunova/bin/illumina-utils/illumina-utils/scripts")
sys.path.append("/bioware/merens-illumina-utils")
# sys.path.append("/bioware/pythonmodules/illumina-utils/")
import fastqlib as fq
import fastalib as fa
from subprocess import call
import ast
from pipeline.utils import Dirs, PipelneUtils
from collections import defaultdict
import constants as C
import getpass

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
        
    def split_files(self, compressed = False):
        """
        TODO: *) fastq_file_names method to collect all file_names with full path or directories_names (see get_all_files()?)
        """   
#        print "compressed = %s" %       compressed
#        compressed = ast.literal_eval(compressed)     
        (in_files_r1, in_files_r2) = self.get_fastq_file_names(self.in_file_path)
        correct_file_names = self.get_correct_file_names(in_files_r1, compressed)
        self.read1(correct_file_names, compressed)
        self.read2(in_files_r2, compressed)
        self.create_inis()
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
   
    def get_all_files(self):
        files = {}
        for dirname, dirnames, filenames in os.walk(self.out_file_path):
            for file_name in filenames:
                full_name = os.path.join(dirname, file_name)
                (file_base, file_extension) = os.path.splitext(os.path.join(dirname, file_name))
                files[full_name] = (file_base, file_extension)
#        print "len(files) = %s" % len(files)
        return files
    
    def perfect_reads(self):
        print "Extract perfect V6 reads:"
        for idx_key in self.runobj.samples.keys():
            file_name = os.path.join(self.out_file_path, idx_key + ".ini")
            program_name = C.perfect_overlap_cmd
            if self.utils.is_local():
                program_name = C.perfect_overlap_cmd_local                    
            try:
                if self.runobj.samples[idx_key].primer_suite.startswith('Archaeal'):
                    call([program_name, file_name, "--archaea"]) 
                else: 
                    call([program_name, file_name])
            except:
                print "Problems with program_name = %s, file_name = %s" % (program_name, file_name)
                raise  
    
    def call_sh_script(self, script_name_w_path, where_to_run):
        try:
            call(['chmod', '0774', script_name_w_path])
            call(['qsub', script_name_w_path], cwd=(where_to_run))
#             pass
        except:
            print "Problems with script_name = %s" % (script_name_w_path)
            raise  
        
    def perfect_reads_cluster(self):
        print "Extract perfect V6 reads:"
        program_name = C.perfect_overlap_cmd
        if self.utils.is_local():
            program_name = C.perfect_overlap_cmd_local
        primer_suite = self.get_config_values('primer_suite')
        if any("Archaeal" in s for s in primer_suite):
            add_arg = " --archaea"
        else: 
            add_arg = ""
        command_line          = program_name + add_arg
        file_list             = self.dirs.get_all_files_by_ext(self.out_file_path, "ini")
        script_file_name      = self.create_job_array_script(command_line, self.dirs.analysis_dir, file_list)
        script_file_name_full = os.path.join(self.dirs.analysis_dir, script_file_name)
        self.call_sh_script(script_file_name_full, self.dirs.analysis_dir)  
        return script_file_name              
                          
    def partial_overlap_reads_cluster(self):
        print "Extract partial_overlap V4V5 reads:"
        program_name = C.partial_overlap_cmd
        if self.utils.is_local():
            program_name = C.partial_overlap_cmd_local       
        dna_region = self.get_config_values('dna_region')
        if ("ITS1" in list(dna_region)):
            add_arg = "--marker-gene-stringent"
        else:
            add_arg = ""
#         TODO: this part is the same in perfect overlap - move into a method    
        command_line          = program_name + " --enforce-Q30-check " + add_arg
        file_list             = self.dirs.get_all_files_by_ext(self.out_file_path, "ini")
        script_file_name      = self.create_job_array_script(command_line, self.dirs.analysis_dir, file_list)
        script_file_name_full = os.path.join(self.dirs.analysis_dir, script_file_name)
        self.call_sh_script(script_file_name_full, self.dirs.analysis_dir)  
        return script_file_name      
                    
    def partial_overlap_reads(self):
        print "Extract partial_overlap V4V5 reads:"
        for idx_key in self.runobj.samples.keys():
            ini_file_name = os.path.join(self.out_file_path, idx_key + ".ini")
            program_name = C.partial_overlap_cmd
            if self.utils.is_local():
                program_name = C.partial_overlap_cmd_local        
            try:
                if (self.runobj.samples[idx_key].dna_region == "ITS1"):
                    call([program_name, "--enforce-Q30-check", "--marker-gene-stringent", ini_file_name])
                else:
                    call([program_name, "--enforce-Q30-check", ini_file_name])
                               
#                 call([program_name, ini_file_name])           
#                 call([program_name, ini_file_name, idx_key])
#                 call([program_name, "--fast-merge", ini_file_name, idx_key])
            except Exception:
#                 except Exception, err:
                print traceback.format_exc()
    #or
#     print sys.exc_info()[0]

                print "Problems with program_name = %s" % (program_name)
                raise  
                
#             print "HERE: program_name = " % (program_name)   
#             call([program_name, "--fast-merge", "--compute-qual-dicts", ini_file_name, idx_key])
            
    def get_config_values(self, key):
        config_path_data = [v for k, v in self.runobj.configPath.items()]
        return set([a[key] for a in config_path_data if key in a.keys()])
        
    def make_users_email(self):
        username = getpass.getuser() 
        return username + "@mbl.edu"
                
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
# Send mail at job end; -m eas sends on end, abort, suspend.
#$ -m eas
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
        print "Filter mismatches if more then %s" % (max_mismatch)
        command_line = C.filter_mismatch_cmd
        if self.utils.is_local():
            command_line = C.filter_mismatch_cmd_local    
        files_dir = self.dirs.reads_overlap_dir   
                
        file_list             = self.dirs.get_all_files_by_ext(files_dir, "_MERGED")
        script_file_name      = self.create_job_array_script(command_line, files_dir, file_list)
        script_file_name_full = os.path.join(files_dir, script_file_name)
        self.call_sh_script(script_file_name_full, files_dir)  
        return script_file_name              

    def filter_mismatches(self, max_mismatch = 3):
        print "Filter mismatches if more then %s" % (max_mismatch)
        n = 0        
        files = self.get_all_files()
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
        print "Uniqueing fasta files"      
        command_line = C.fastaunique_cmd
        if self.utils.is_local():
            command_line = C.fastaunique_cmd_local   
        files_dir = self.dirs.reads_overlap_dir   
                
        file_list             = self.dirs.get_all_files_by_ext(files_dir, C.filtered_suffix)
        if len(file_list) == 0:
            file_list         = self.dirs.get_all_files_by_ext(files_dir, ".fa")
        script_file_name      = self.create_job_array_script(command_line, files_dir, file_list)
        script_file_name_full = os.path.join(files_dir, script_file_name)
        self.call_sh_script(script_file_name_full, files_dir)  
        return script_file_name                           
                                       
    def uniq_fa(self):
        n = 0        
        print "Uniqueing fasta files"      
        files = self.get_all_files()
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
            if self.runobj.samples[idx_key].primer_suite in C.primers_dict:
                proximal_primer = C.primers_dict[self.runobj.samples[idx_key].primer_suite]["proximal_primer"]
                distal_primer = C.primers_dict[self.runobj.samples[idx_key].primer_suite]["distal_primer"]
#                 print "proximal_primer: %s. distal_primer: %s" % (proximal_primer, distal_primer)
            else:
                print "ERROR! Something wrong with the primer suite name: %s. NB: For v6mod it suppose to be 'Archaeal V6mod Suite'" % (self.runobj.samples[idx_key].primer_suite)
            primers[idx_key] = (proximal_primer, distal_primer) 
            
        return primers
        
    def create_inis(self):
        for idx_key in self.runobj.samples.keys():
            run_key = idx_key.split('_')[1].replace("N", ".");
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

            "That's for parital overlap (v4v5 miseq illumina)" 
            if not self.runobj.do_perfect:
                primers = self.get_primers()    
                text += """
# following section is optional
[prefixes]
pair_1_prefix = ^""" + run_key + primers[idx_key][0] + "\npair_2_prefix = ^" + primers[idx_key][1]
                
            ini_file_name = os.path.join(self.out_file_path,  idx_key + ".ini")
            self.open_write_close(ini_file_name, text)

    def open_write_close(self, script_file_name, text):
        ini_file = open(script_file_name, "w")
        ini_file.write(text)
        ini_file.close()
 
    def get_fastq_file_names(self, f_input_file_path):
        in_files_r1 = []
        in_files_r2 = []
        "TODO: exclude dir with new created files from the loop"
        for dirname, dirnames, filenames in os.walk(f_input_file_path):
            for filename in filenames:
                if filename.find('_R1_') > 0:
                    in_files_r1.append(os.path.join(dirname, filename))
                elif filename.find('_R2_') > 0:
                    in_files_r2.append(os.path.join(dirname, filename))
                else:
                    sys.stderr.write("No read number in the file name: %s\n" % filename)
        return (in_files_r1, in_files_r2)
    
    def get_correct_file_names(self, file_r1, compressed):
        correct_file_names = [];
        for file_r1 in file_r1:
            print "FFF1: file %s" % file_r1
            index_sequence = self.get_index(file_r1)
#             self.runobj.run_keys
#             
            good_run_key_lane_names = [x for x in self.runobj.run_keys if x.startswith(index_sequence)]
            if len(good_run_key_lane_names) > 0:
                correct_file_names.append(file_r1)
        return set(correct_file_names)
        
    def read1(self, files_r1, compressed):
        """ loop through the fastq_file_names
            1) e.pair_no = 1, find run_key -> dataset name
            2) collect the relevant part of id
        """
        for file_r1 in files_r1:
            print "====\nFFF1: file %s" % file_r1
            f_input  = fq.FastQSource(file_r1, compressed)
            while f_input.next():
                e = f_input.entry
                ini_run_key  = e.index_sequence + "_" + "NNNN" + e.sequence[4:9] + "_" + e.lane_number                
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
                    
    def read2(self, files_r2, compressed):
        "3) e.pair_no = 2, find id from 2), assign dataset_name"
        for file_r2 in files_r2:
            print "FFF2: file %s" % file_r2
            f_input  = fq.FastQSource(file_r2, compressed)
            while f_input.next():
                e = f_input.entry
                
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
