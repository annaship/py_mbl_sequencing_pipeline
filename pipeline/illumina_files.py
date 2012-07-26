import sys
import MySQLdb
import os
from pipeline.utils import PipelneUtils
sys.path.append("/bioware/pythonmodules/fastqlib")
import fastqlib as fq
from collections import defaultdict

class IlluminaFiles:
    """
    0) from run create all dataset_lines names files in output dir
    1) split fastq files from casava into files with dataset_names
    2) create in files 
    3) process them through Meren's script
    4) result - files dataset_lane-PERFECT_reads.fa.unique with frequencies - to process with env454upload()
    
    
    """
    def __init__(self, run):
       self.fastq_dir        = os.path.join(run.input_dir, "fastq/")
       self.run              = run
       self.out_file_names   = {} 
       self.id_dataset       = {}
       self.datasets         = list(set([self.run.samples[key].dataset for key in self.run.samples]))
       self.output_file_path = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/analysis"
       self.open_dataset_files()


    def open_dataset_files(self):
        for dataset in self.datasets + ["unknown"]:
            output_file = os.path.join(self.output_file_path, dataset + ".fastq")
#            self.out_file_names[dataset] = fq.FastQOutput(output_file)

    def close_dataset_files(self):
        for dataset in self.datasets + ["unknown"]:
#            output_file = os.path.join(self.output_file_path, dataset + ".fastq")
#            self.out_file_names[dataset].close
            pass

    
#    """
#        while input.next():
#            self.out_file_names.store(input)
#            
#        for dataset in datasets + ["unknown"]:
#            self.out_file_names[dataset].close
#    """
    
    
    def split_files(self, input_file_path, output_file_path, compressed = False):
#        input_file_path = self.fastq_dir
        """
        TODO: 1) path should be argument, not hard-coded!
              2) loop through directories, until got files recursively
        """
        input_file_path  = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/" 
        output_file_path = self.output_file_path
        print "input_file_path = %s, output_file_path = %s, compressed = %s\n" % (input_file_path, output_file_path, compressed)
#        a = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/Project_Julie_v6_30/Sample_v6_Amplicon_IDX1"
        "TODO: fastq_file_names method to collect all file_names with full path or directories_names"
        self.collect_file_names(input_file_path)
        self.read1()
        self.read2()
#        (dirpath1, dirname1, files1) = self.get_fastq_file_names(input_file_path)
        
#        if dirname1:
#            for dirname2 in dirname1:
#                print "dirname2 = %s in dirname1 = %s" % (dirname2, dirname1)
##            dirname = [x for x in range(20) if x % 2 == 0]
#                dirname3 = os.path.join(dirpath1, dirname2)
#                (dirpath, dirname, files) = self.get_fastq_file_names(dirname3)
#
#                for file in files:
#                    print "\n=============\nfile = %s" % file
#                    input  = fq.FastQSource(os.path.join(dirpath, file), compressed)
#
##                    output = fq.FileOutput(os.path.join(output_file_path, file), compressed)
#                
##                    oo = open(os.path.join(dirpath, file))
##                    print oo.readline()
#
#                    """TODO: loop through the fastq_file_names
#                        1) e.pair_no = 1, find run_key -> dataset name
#                        2) collect the relevant part of id
#                        3) e.pair_no = 2, find id from 2), assign dataset_name
#                    """
#                    while input.next():
#                        e = input.entry
##                        values = list(set([run_key.split('_')[1] for run_key in self.run.run_keys]))
#                        ini_run_key  = e.index_sequence + "_" + "NNNN" + e.sequence[4:9] + "_" + e.lane_number
#                        if ini_run_key in self.run.samples.keys():
#                            sample       = self.run.samples[ini_run_key] 
#                            dataset_name = sample.dataset
#                            self.out_file_names[dataset_name].store_entry(e)
#                            
#                            "TODO: make a method:"
#                            short_id1 = e.header_line.split()[0]
#                            short_id2 = ":".join(e.header_line.split()[1].split(":")[1:])
#                            id2 = short_id1 + " 2:" + short_id2
##                            print  "\n=========\nid2 = %s\nshort_id1 = %s; short_id2 = %s\n=========\n" % (id2, short_id1, short_id2)                     
#                            self.id_dataset[id2] = dataset_name
#                            print "e.x_coord = %s, e.pair_no = %s, e.lane_number = %s, e.index_sequence = %s\ne.header_line=%s\ne.sequence = %s" % (e.x_coord, e.pair_no, e.lane_number, e.index_sequence, e.header_line, e.sequence)
#                        elif (e.pair_no == 2):
#                            """TODO: how to run all pairs1 first!"""
#                            dataset_name = self.id_dataset[e.header_line]
#                            print "HERE: dataset_name = %s" % (dataset_name)
#                            self.out_file_names[dataset_name].store_entry(e)                            
#                        else:
#                            self.out_file_names["unknown"].store_entry(e)
                            

#        #            output.write('>%s\n%s\n' % ('.'.join([e.machine_name,
#        #                                                  e.run_id,
#        #                                                  e.x_coord,
#        #                                                  e.y_coord,
#        #                                                  e.pair_no]),
#        #                                                 e.sequence))
#            
#        #        sys.stderr.write('\n')
#                        input.close()
#                        output.close()
            
        return

 
    def get_fastq_file_names(self, input_file_path):
        in_files_r1 = []
        in_files_r2 = []
        for dirname, dirnames, filenames in os.walk(input_file_path):
            for filename in filenames:
                print "\n=====\nfilename = %s" % filename
                if filename.find('_R1_') > 0:
                    print "R1: filename.find('_R1_') = %s" % filename.find('_R1_')
                    in_files_r1.append(os.path.join(dirname, filename))
                elif filename.find('_R2_') > 0:
                    print "R2"
                    in_files_r2.append(os.path.join(dirname, filename))
                else:
                    print "No read number in the file name: %s" % filename
        return (in_files_r1, in_files_r2)
    
    def collect_file_names(self, input_file_path):
        (in_files_r1, in_files_r2) = self.get_fastq_file_names(input_file_path)
        print (in_files_r1, in_files_r2)
#        (dirpath1, dirname1, files1) = self.get_fastq_file_names(input_file_path)
#        print "(dirpath1 = %s, dirname1 = %s, files1 = %s)" % (dirpath1, dirname1, files1)
#       
#        if dirname1:
#            self.collect_file_names(dirpath1)
#            for dirname2 in dirname1:
##            dirname = [x for x in range(20) if x % 2 == 0]
#                dirname3 = os.path.join(dirpath1, dirname2)
#                (dirpath, dirname, files) = self.get_fastq_file_names(dirname3)

        
        pass
    def read1(self):
        pass
    def read2(self):
        pass