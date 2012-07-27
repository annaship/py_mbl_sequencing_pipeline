import sys
import os
sys.path.append("/bioware/pythonmodules/fastqlib")
import fastqlib as fq
#from collections import defaultdict

class IlluminaFiles:
    """
    0) from run create all dataset_lines names files in output dir
    1) split fastq files from casava into files with dataset_names
    2) create in files 
    3) process them through Meren's script
    4) result - files dataset_lane-PERFECT_reads.fa.unique with frequencies - to process with env454upload()
    
    
    """
    def __init__(self, run):
        self.fastq_dir     = os.path.join(run.input_dir, "fastq/")
        self.run           = run
        self.out_files     = {} 
        self.id_dataset    = {}
        self.datasets      = list(set([self.run.samples[key].dataset for key in self.run.samples]))
        self.dataset_emails = dict((self.run.samples[key].dataset, self.run.samples[key].email) for key in self.run.samples)
#        aDict = dict((x, x+1) for x in a if x % 2 ==0)

        self.out_file_path = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/analysis"
        self.open_dataset_files()
        self.total_seq = 0


    def open_dataset_files(self):
        for dataset in self.datasets + ["unknown"]:
            output_file = os.path.join(self.out_file_path, dataset + ".fastq")
            self.out_files[dataset] = fq.FastQOutput(output_file)

    def close_dataset_files(self):
        for dataset in self.datasets + ["unknown"]:
            self.out_files[dataset].close
            pass

    
#    """
#        while f_input.next():
#            self.out_files.store(f_input)
#            
#        for dataset in datasets + ["unknown"]:
#            self.out_files[dataset].close
#    """
    
    
    def split_files(self, f_in_dir_path, f_out_dir_path, compressed = False):
#        f_input_file_path = self.fastq_dir
        """
        TODO: 1) path should be argument, not hard-coded!
              2) loop through directories, until got files recursively
        """
        f_in_dir_path  = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/" 
        "TODO: fastq_file_names method to collect all file_names with full path or directories_names"
        (in_files_r1, in_files_r2) = self.get_fastq_file_names(f_in_dir_path)
        self.read1(in_files_r1, compressed)
        self.read2(in_files_r2, compressed)
        print "TTT: total_seq = %s" % self.total_seq
        self.create_inis(f_in_dir_path, self.out_file_path)
#        return

    def create_inis(self, f_in_dir_path, f_out_dir_path):
        """[general]
        project_name = x4 ???ASK Meren why project, not dataset???
        researcher_email = meren@mbl.edu
        input_directory = /xraid2-2/sequencing/Illumina/20120613/JulieR
        output_directory = /xraid2-2/sequencing/Illumina/20120613/JulieR/analysis
        
        [files]
        lane_1 = x4.fastq
         """
        for dataset in self.dataset_emails.keys():
#            print "dataset = %s, self.dataset_emails[dataset] = %s" % (dataset, self.dataset_emails[dataset])
            "TODO: dataset+\".fastq\" should real name, take from creation, not create here"
            text = """
[general]
project_name = %s
researcher_email = %s
input_directory = %s
output_directory = %s

[files]
lane_1 = %s
            """ % (dataset, self.dataset_emails[dataset], f_in_dir_path, f_out_dir_path, dataset+".fastq")
#            """ % (dataset, self.dataset_emails[dataset], f_in_dir_path, f_out_dir_path, self.out_files[dataset].file_path)
            print "SSS = %s" % text
            ini_file_name = os.path.join(f_out_dir_path,  dataset+".ini")
            ini_file = open(ini_file_name, "w")
            ini_file.write(text)
            ini_file.close()
        
#        self.datasets      = list(set([self.run.samples[key].dataset for key in self.run.samples]))

 
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
            f_input  = fq.FastQSource(file_r1, compressed)
            while f_input.next():
                e = f_input.entry
                ini_run_key  = e.index_sequence + "_" + "NNNN" + e.sequence[4:9] + "_" + e.lane_number
                if ini_run_key in self.run.samples.keys() and int(e.pair_no) == 1:
                    sample       = self.run.samples[ini_run_key] 
                    dataset_name = sample.dataset
                    self.out_files[dataset_name].store_entry(e)
                    self.collect_dataset_id()
                    "TODO: make a method:"
                    short_id1 = e.header_line.split()[0]
                    short_id2 = ":".join(e.header_line.split()[1].split(":")[1:])
                    id2 = short_id1 + " 2:" + short_id2
                    self.id_dataset[id2] = dataset_name
                    self.total_seq +=1           
                else:
                    self.out_files["unknown"].store_entry(e)
                    
    def read2(self, files_r2, compressed):
        "3) e.pair_no = 2, find id from 2), assign dataset_name"
        for file_r2 in files_r2:
            f_input  = fq.FastQSource(file_r2, compressed)

            while f_input.next():
                e = f_input.entry
                self.total_seq +=1                
                
                if int(e.pair_no) == 2:
                    dataset_name = self.id_dataset[e.header_line]
                    self.out_files[dataset_name].store_entry(e)        
                else:
                    self.out_files["unknown"].store_entry(e)

    def collect_dataset_id(self):
        pass
