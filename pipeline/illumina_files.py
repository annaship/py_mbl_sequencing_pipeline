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
       self.datasets         = list(set([self.run.samples[key].dataset for key in self.run.samples]))
       self.output_file_path = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/analysis"
       self.open_dataset_files()


    def open_dataset_files(self):
        for dataset in self.datasets + ["unknown"]:
            output_file = os.path.join(self.output_file_path, dataset + ".fastq")
            self.out_file_names[dataset] = fq.FastQOutput(output_file)

    def close_dataset_files(self):
        for dataset in self.datasets + ["unknown"]:
#            output_file = os.path.join(self.output_file_path, dataset + ".fastq")
            self.out_file_names[dataset].close

    
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
        input_file_path  = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/Project_Julie_v6_30" 
        output_file_path = self.output_file_path
        print "input_file_path = %s, output_file_path = %s, compressed = %s\n" % (input_file_path, output_file_path, compressed)
#        a = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/Project_Julie_v6_30/Sample_v6_Amplicon_IDX1"
        (dirpath1, dirname1, files1) = self.get_fastq_file_names(input_file_path)
        "TODO: fastq_file_names method to collect all file_names with full path or directories_names"
        if dirname1:
            for dirname2 in dirname1:
                print "dirname2 = %s in dirname1 = %s" % (dirname2, dirname1)
#            dirname = [x for x in range(20) if x % 2 == 0]
                dirname3 = os.path.join(dirpath1, dirname2)
                (dirpath, dirname, files) = self.get_fastq_file_names(dirname3)

                for file in files:
                    input  = fq.FastQSource(os.path.join(dirpath, file), compressed)

#                    output = fq.FileOutput(os.path.join(output_file_path, file), compressed)
                
#                    oo = open(os.path.join(dirpath, file))
#                    print oo.readline()

                    """TODO: loop through the fastq_file_names
                        1) e.pair_no = 1, find run_key -> dataset name
                        2) collect the relevant part of id
                        3) e.pair_no = 2, find id from 2), assign dataset_name
                    """
                    while input.next():
                        e = input.entry
#                        values = list(set([run_key.split('_')[1] for run_key in self.run.run_keys]))
                        ini_run_key  = e.index_sequence + "_" + "NNNN" + e.sequence[4:9] + "_" + e.lane_number
                        if ini_run_key in self.run.samples.keys():
                            sample       = self.run.samples[ini_run_key] 
                            dataset_name = sample.dataset
                            self.out_file_names[dataset_name].store_entry(e)
    
    #                        print "real_run_key = %s" % real_run_key
                            print "e.x_coord = %s, e.pair_no = %s, e.lane_number = %s, e.index_sequence = %s\nself.header_line=%s\ne.sequence = %s" % (e.x_coord, e.pair_no, e.lane_number, e.index_sequence, e.header_line, e.sequence)
                        else:
                            self.out_file_names["unknown"].store_entry(e)

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

 
    def get_fastq_file_names(self, fastq_dir):
        for (dirpath, dirname, files) in os.walk(fastq_dir):
            return (dirpath, dirname, files)
    