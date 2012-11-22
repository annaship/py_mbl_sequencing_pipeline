import sys
import os
sys.path.append("/xraid/bioware/linux/seqinfo/bin")
sys.path.append("/Users/ashipunova/bin/illumina-utils")
import fastqlib as fq
import fastalib as fa
from subprocess import call
import ast

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
        if os.uname()[1] == 'ashipunova.mbl.edu' or os.uname()[1] == "as-macbook.home" or os.uname()[1] == "as-macbook.local" or os.uname()[1] == "Ashipunova.local":
            self.LOCAL = True
        else:
            self.LOCAL = False
        self.runobj         = runobj
        self.out_files      = {} 
        self.id_dataset_idx = {}
        self.dataset_emails = dict((self.runobj.samples[key].dataset, self.runobj.samples[key].email) for key in self.runobj.samples)
        self.dataset_index  = dict((self.runobj.samples[key].dataset + "_" + self.runobj.samples[key].barcode_index, self.runobj.samples[key].dataset) for key in self.runobj.samples)
        self.in_file_path   = self.runobj.input_dir
        self.out_file_path  = self.create_out_dir(os.path.join(self.runobj.output_dir, "analysis"))
#        self.in_file_path  = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/input/illumina_files_test"
#        self.out_file_path = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/output/analysis"
        self.results_path  = self.create_out_dir(os.path.join(self.out_file_path, "perfect_reads"))
        self.open_dataset_files()
        
    def split_files(self, compressed = False):
        """
        TODO: *) fastq_file_names method to collect all file_names with full path or directories_names (see get_all_files()?)
        """   
#        print "compressed = %s" %       compressed
#        compressed = ast.literal_eval(compressed)     
        (in_files_r1, in_files_r2) = self.get_fastq_file_names(self.in_file_path)
        self.read1(in_files_r1, compressed)
        self.read2(in_files_r2, compressed)
        self.create_inis()
        self.close_dataset_files()

#        self.perfect_reads()
#        self.uniq_fa()
            
    def create_out_dir(self, dirname):
        try:
            os.makedirs(dirname)
        except OSError:
            if os.path.isdir(dirname):
                pass
            else:
                # There was an error on creation, so make sure we know about it
                raise        
        return dirname

    def open_dataset_files(self):
        n = 0
        for dataset_idx in self.dataset_index.keys():
            for read_n in ["_R1", "_R2"]:
                n += 1
                dataset_idx_base = dataset_idx + read_n
                output_file = os.path.join(self.out_file_path, dataset_idx_base + ".fastq")
                self.out_files[dataset_idx_base] = fq.FastQOutput(output_file)
        self.out_files["unknown"] = fq.FastQOutput(os.path.join(self.out_file_path, "unknown" + ".fastq"))
        

    def close_dataset_files(self): 
#        map(lambda x: x.close(), [dict[item]['filename'] for item in dict])
#        for dataset in self.dataset_emails.keys():
#            for read_n in ["_R1", "_R2"]:
#                key_d = dataset + read_n
#                self.out_files[key_d].close()
#        self.out_files["unknown"].close()
        for o_file in self.out_files.iteritems():
            o_file[1].close()
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
    
    "TODO: made one method out of the next two"
    def perfect_reads(self):
        n = 0
        print "Extract perfect V6 reads:"
        for dataset_idx in self.dataset_index.keys():
            file_name = os.path.join(self.out_file_path, dataset_idx + ".ini")
            n +=1
#            print "%s ini file: %s" % (n, file_name)
            program_name = "analyze-illumina-v6-overlaps"
            if self.LOCAL:
                program_name = "/Users/ashipunova/bin/illumina-utils/analyze-illumina-v6-overlaps"           
            "TODO: get as parameter!"
#            self.archaea = True
            self.archaea = False
            if self.archaea:
                call([program_name, file_name, "--archaea"]) 
            else: 
                call([program_name, file_name])
    
    def uniq_fa(self):
        n = 0        
        print "Uniqueing fasta files"      
        files = self.get_all_files()
        for full_name in files.keys():    
            if files[full_name][1] == ".fa":
                n +=1   
#                print "%s fasta file: %s" % (n, full_name)
                program_name = "fastaunique"
                if self.LOCAL:
                    program_name = "/Users/ashipunova/bin/illumina-utils/fastaunique"                
                call([program_name, full_name])
#                pass

#                    
#            fasta = fa.SequenceSource(full_name, unique = True) 
#            while fasta.next():
#                e = fasta.entry
#                """TODO: open all files with name -PERFECT_reads.fa.unique
#                files[full_name][0] + '.fa.unique'
#                """
#                self.out_files["W5_4-PERFECT_reads.fa.unique"].store_entry(e)

    def create_inis(self):
        for dataset_idx in self.dataset_index.keys():
            email = self.dataset_emails[self.dataset_index[dataset_idx]]
#        for dataset in self.dataset_emails.keys():
#            dataset_idx_base = dataset + "_" + self.dataset_index[dataset]
#            print "dataset = %s, self.dataset_emails[dataset] = %s" % (dataset, self.dataset_emails[dataset])
            "TODO: one argument"
            text = """[general]
project_name = %s
researcher_email = %s
input_directory = %s
output_directory = %s

[files]
pair_1 = %s
pair_2 = %s
            """ % (dataset_idx, email, self.out_file_path, self.results_path, dataset_idx + "_R1.fastq", dataset_idx + "_R2.fastq")
            ini_file_name = os.path.join(self.out_file_path,  dataset_idx + ".ini")
            self.open_write_close(ini_file_name, text)

    def open_write_close(self, ini_file_name, text):
        ini_file = open(ini_file_name, "w")
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
        
    def read1(self, files_r1, compressed):
        """ loop through the fastq_file_names
            1) e.pair_no = 1, find run_key -> dataset name
            2) collect the relevant part of id
        """
        for file_r1 in files_r1:
            print "FFF1: file %s" % file_r1
            f_input  = fq.FastQSource(file_r1, compressed)
            while f_input.next():
                e = f_input.entry
                ini_run_key  = e.index_sequence + "_" + "NNNN" + e.sequence[4:9] + "_" + e.lane_number                
#                ini_run_key  = e.index_sequence + "_" + e.sequence[4:9] + "_" + e.lane_number
                if ini_run_key in self.runobj.samples.keys() and int(e.pair_no) == 1:
                    sample = self.runobj.samples[ini_run_key] 
                    dataset_file_name_base = sample.dataset + "_" + sample.barcode_index
                    dataset_file_name_base_r1 = dataset_file_name_base + "_R1"
                    if (dataset_file_name_base_r1 in self.out_files.keys()):
                        self.out_files[dataset_file_name_base_r1].store_entry(e)
    #                    print "id = %s,\nseq = %s" % (e.id, e.sequence)
                        "TODO: make a method:"
                        short_id1 = e.header_line.split()[0]
                        short_id2 = ":".join(e.header_line.split()[1].split(":")[1:])
                        id2 = short_id1 + " 2:" + short_id2
                        self.id_dataset_idx[id2] = dataset_file_name_base
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

