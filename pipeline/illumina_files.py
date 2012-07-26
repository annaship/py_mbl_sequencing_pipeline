import sys
import MySQLdb
import os
from pipeline.utils import PipelneUtils
sys.path.append("/bioware/pythonmodules/fastqlib")
import fastqlib as fq


class IlluminaFiles:
    """
    1) split fastq files from casava into files with dataset_names
    2) process them through Meren's script
    3) result - files dataset_lane-PERFECT_reads.fa.unique with frequencies - to process with env454upload()
    """
    def __init__(self, run):
       self.fastq_dir   = os.path.join(run.input_dir, "fastq/")
       self.run = run
       
    def split_files(self, input_file_path, output_file_path, compressed = False):
#        input_file_path = self.fastq_dir
        """
        TODO: path should be argument, not hard-coded!
        """
        input_file_path  = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/Project_Julie_v6_30" 
        output_file_path = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/analysis"
        print "input_file_path = %s, output_file_path = %s, compressed = %s\n" % (input_file_path, output_file_path, compressed)
#        a = "/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/data/fastq/illumina_files_test/Project_Julie_v6_30/Sample_v6_Amplicon_IDX1"
        (dirpath1, dirname1, files1) = self.get_fastq_file_names(input_file_path)
        if dirname1:
            for dirname2 in dirname1:
                print "dirname2 = %s in dirname1 = %s" % (dirname2, dirname1)
#            dirname = [x for x in range(20) if x % 2 == 0]
                dirname3 = os.path.join(dirpath1, dirname2)
                (dirpath, dirname, files) = self.get_fastq_file_names(dirname3)

                for file in files:
                    input  = fq.FastQSource(os.path.join(dirpath, file), compressed)
                    output = fq.FileOutput(os.path.join(output_file_path, file), compressed)
                
#                    oo = open(os.path.join(dirpath, file))
#                    print oo.readline()

                    while input.next():
                        e = input.entry
                        print "e.x_coord = %s, e.pair_no = %s, e.lane_number = %s, e.index_sequence = %s\nself.header_line=%s\ne.sequence = %s" % (e.x_coord, e.pair_no, e.lane_number, e.index_sequence, e.header_line, e.sequence)
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
            print (dirpath, dirname, files)
            return (dirpath, dirname, files)