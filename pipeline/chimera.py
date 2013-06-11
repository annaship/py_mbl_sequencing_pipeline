import subprocess
import sys, os
import re
import time
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils
from pprint import pprint
from collections import defaultdict

import pipeline.constants as C

class Chimera:
    """ Define here """
    def __init__(self, runobj = None):
        self.utils      = PipelneUtils()
        self.runobj     = runobj
        self.run_keys   = self.runobj.run_keys
        self.rundate    = self.runobj.run
        self.chg_suffix = ".chg"
        self.ref_suffix = ".chimeras.db"      
        self.denovo_suffix      = ".chimeras.txt"        
        self.nonchimeras_suffix = ".nonchimeras.fa"

        if self.runobj.vamps_user_upload:
            site       = self.runobj.site
            dir_prefix = self.runobj.user + '_' + self.runobj.run
        else:
            site = ''
            dir_prefix = self.runobj.run
        if self.runobj.lane_name:
            lane_name = self.runobj.lane_name
        else:
            lane_name = ''
        
        self.dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, lane_name = lane_name, site = site) 
        self.indir  = self.dirs.check_dir(self.dirs.reads_overlap_dir)
        self.outdir = self.dirs.check_dir(self.dirs.chimera_dir)
        self.usearch_cmd = C.usearch_cmd
        self.abskew      = C.chimera_checking_abskew
        self.refdb       = C.chimera_checking_refdb
        self.its_refdb   = C.chimera_checking_its_refdb
        self.input_file_names  = self.make_chimera_input_illumina_file_names()
#         pprint(self.run_keys)
#         self.output_file_names = self.make_chimera_output_illumina_file_names(self.input_file_names)

    def make_chimera_input_illumina_file_names(self):
        input_file_names = {} 
        
        for idx_key in self.run_keys:
            file_name = idx_key + "_" + C.filtered_suffix + ".unique" 
           
            if os.path.exists(os.path.join(self.indir, file_name)):
                input_file_names[idx_key] = file_name
        
        return input_file_names
            
#     def make_chimera_output_illumina_file_names(self, input_file_names):
#         output_file_names = {} 
#         for idx_key, input_file_name in input_file_names.iteritems():
#             output_file_names[idx_key] = input_file_name
#         return output_file_names

    def get_current_dirname(self, in_or_out = ""):
        if in_or_out == "":
            cur_dirname    = self.indir 
        else:
            cur_dirname    = self.outdir
        return cur_dirname

    def is_chimera_check_file(self, filename):
        return filename.endswith((self.denovo_suffix, self.ref_suffix, self.nonchimeras_suffix))

    def get_current_filenames(self, cur_dirname):
        cur_file_names = []
        if cur_dirname == self.indir:
            cur_file_names = self.input_file_names.values()
        elif cur_dirname == self.outdir:
            for dirname, dirnames, filenames in os.walk(cur_dirname):
                for filename in filenames:
                    if (self.is_chimera_check_file(filename)):
                        print "filename = %s" % filename
                        cur_file_names.append(filename)        
        return cur_file_names

    def illumina_frequency_size(self, in_or_out = "", find = "frequency:", replace = ";size="):
        cur_dirname    = self.get_current_dirname(in_or_out)
        cur_file_names = self.get_current_filenames(cur_dirname)
#         print "cur_file_names: "
#         pprint(cur_file_names)
        change_from_suffix = ""
        change_to_suffix   = self.chg_suffix
#         print "find = %s, replace = %s" % (find, replace)
        regex              = re.compile(r"%s" % find)

        for cur_file_name in cur_file_names:
            file_name = os.path.join(cur_dirname, cur_file_name)
            with open(file_name + change_from_suffix, "r") as sources:
                lines = sources.readlines()
            with open(file_name + change_to_suffix, "w") as target:
                for line in lines:
                        target.write(regex.sub(replace, line))
                    
#     def move_out_chimeric(self):
#         chimeric_ids = self.get_chimeric_ids()
#         pprint(chimeric_ids)
#       
#     def get_chimeric_ids(self):
#       
# # http://drive5.com/uchime/uchime_quickref.pdf
# # The --uchimeout file is a tab-separated file with the following 17 fields.
# # Field Name Description
# # 1 Score Value >= 0.0, high score means more likely to be a chimera.
# # 2 Query Sequence label
# # 3 Parent A Sequence label
# # 4 Parent B Sequence label
# # 5 IdQM %id between query and model made from (A, crossover, B)
# # 6 IdQA %id between query and parent A.
# # 7 IdQB %id between query and parent B
# # 8 IdAB %id between parents (A and B).
# # 9 IdQT %id between query and closest reference sequence / candidate parent.
# # 10 LY Yes votes on left
# # 11 LN No votes on left
# # 12 LA Abstain votes on left
# # 13 RY Yes votes on right
# # 14 RN No votes on right
# # 15 RA Abstain votes on right
# # 16 Div Divergence ratio, i.e. IdQM - IdQT
# # 17 YN Y (yes) or N (no) classification as a chimera. Set to Y if score >= threshold
#         chimeric_ids = defaultdict(list)
# 
#         chimeric_files = self.get_current_filenames(self.outdir)
# #         pprint(chimeric_files)
#         
#         for chimeric_file in chimeric_files:
#             file_name = os.path.join(self.outdir, chimeric_file + self.chg_suffix)
#             with open(file_name, "r") as sources:
#                 lines = sources.readlines()[1:]
#             # make a list of chimera deleted read_ids   
#             for line in lines:
# #                 print "line = %s" % line
#                 line_list_tab = line.strip().split("\t")
#                 chimera_yesno = line_list_tab[-1]
#                 if(chimera_yesno) == 'Y':
#                     id = line_list_tab[1]
#                     chimeric_ids[chimeric_file].append(id)           
# #         print "chimeric_ids:"
#         return chimeric_ids 

#         for idx_key in self.run_keys:
#             # open  deleted file and append chimera to it
#             # open and read both chimeras files: chimeras.db and chimeras.txt
#             
#             # hash to remove dupes
#             chimera_deleted = {}
#             for file in [self.files[idx_key]['chimera_db'], self.files[idx_key]['chimera_txt']]:            
#                 fh = open(file,"r") 
#                 # make a list of chimera deleted read_ids            
#                 for line in fh.readlines():
#                     lst = line.strip().split()
#                     id = lst[1].split(';')[0]
#                     chimera_yesno = lst[-1]
#                     if(chimera_yesno) == 'Y':
#                         chimera_deleted[id] = 'chimera'
#         
#             fh_del = open(self.files[idx_key]['deleted'],"a")
#             for id in chimera_deleted:
#                 fh_del.write(id+"\tchimera\n") 
#                     

    def illumina_rm_size_files(self):
        for idx_key in self.input_file_names:
            file_name = os.path.join(self.indir, self.input_file_names[idx_key] + self.chg_suffix)
            if os.path.exists(file_name):
                os.remove(file_name)
          
    def check_if_cluster_is_done(self):
        cluster_done = False
        check_qstat_cmd_line = "qstat | grep usearch | wc -l"
#         check_qstat_cmd_line = "qstat | grep usearch"

        print "check_qstat_cmd_line = %s" % check_qstat_cmd_line
        
        try:
            p = subprocess.Popen(check_qstat_cmd_line, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            num_proc = int(output)
            print "qstat is running %s 'usearch' processes" % num_proc
    #         pprint(p)
            
            if (num_proc == 0):
                cluster_done = True
    #         print "cluster_done from check_if_cluster_is_done = %s" % cluster_done
        except:
            print "Chimera checking can be dan only on grendel."
            raise

        return cluster_done
        
          
    def create_chimera_cmd(self, input_file_name, output_file_name, ref_or_novo, ref_db = ""):
#         from usearch -help
# Chimera detection (UCHIME ref. db. mode):
#   usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
#     [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]
# 
# Chimera detection (UCHIME de novo mode):
#   usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
#      [-uchimeout results.uch] [-uchimealns results.alns]
#   Input is estimated amplicons with integer abundances specified using ";size=N".

        
        uchime_cmd = C.clusterize_cmd
        uchime_cmd += " "
        uchime_cmd += self.usearch_cmd
        uchime_cmd += " --uchime "
        uchime_cmd += input_file_name
 
        if (ref_or_novo == "denovo"):
            output_file_name = output_file_name + self.denovo_suffix 
            cmd_append = " --abskew " + self.abskew                                   
            
        elif (ref_or_novo == "ref"):
            output_file_name = output_file_name + self.ref_suffix           
            cmd_append = " --db " + ref_db   
        else:
            print "Incorrect method, should be \"denovo\" or \"ref\"" 
                                        
        uchime_cmd += " --uchimeout "
        uchime_cmd += output_file_name
        uchime_cmd += " --nonchimeras "
        uchime_cmd += (output_file_name + self.nonchimeras_suffix)
         
        uchime_cmd += cmd_append
#         print "uchime_cmd FROM create_chimera_cmd = %s" % (uchime_cmd)
        return uchime_cmd
        
    def get_ref_db(self, dna_region):
        ref_db = ''
        if dna_region.upper() == 'ITS':
            logger.debug("got an ITS dna region so using refdb: " + self.its_refdb)
            ref_db = self.its_refdb
        else:
            logger.debug("using standard refdb: " + self.refdb)
            ref_db = self.refdb
        return ref_db       
    
    def chimera_checking(self, ref_or_novo):
        chimera_region_found = False
        output = {}
        
        for idx_key in self.input_file_names:
#             print "idx_key, self.input_file_names[idx_key] = %s, %s" % (idx_key, self.input_file_names)
            input_file_name  = os.path.join(self.indir,  self.input_file_names[idx_key] + self.chg_suffix)        
            output_file_name = os.path.join(self.outdir, self.input_file_names[idx_key])        
            dna_region       = self.runobj.samples[idx_key].dna_region
#             print "dna_region = %s" % dna_region
            if dna_region in C.regions_to_chimera_check:
                chimera_region_found = True
            else:
                logger.debug('region not checked: ' +  dna_region)
                continue
            
#             print "input_file_name = %s \noutput_file_name = %s" % (input_file_name, output_file_name)
            ref_db     = self.get_ref_db(dna_region)
#             print "dna_region = %s; ref_db = %s; ref_or_novo = %s" % (dna_region, ref_db, ref_or_novo)
            
            uchime_cmd = self.create_chimera_cmd(input_file_name, output_file_name, ref_or_novo, ref_db)
            print "\n==================\n%s command: %s" % (ref_or_novo, uchime_cmd)
            
            try:
                logger.info("chimera checking command: " + str(uchime_cmd))
                output[idx_key] = subprocess.Popen(uchime_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            except OSError, e:
                print "Problems with this command: %s" % (uchime_cmd)
                if self.utils.is_local():
                    print >>sys.stderr, "Execution of %s failed: %s" % (uchime_cmd, e)
                else:
                    print >>sys.stderr, "Execution of %s failed: %s" % (uchime_cmd, e)
                    raise                  
                               
# ???
        if not chimera_region_found:            
            return ('NOREGION', 'No regions found that need checking', '')
        else:
            return ("The usearch commands were created")
    
    """ For 454.
        not tested 
    """
    def chimera_denovo(self):
        chimera_region_found = False
        output = {}
        cluster_id_list = []
         
        for idx_key in self.input_file_names:
#             print "idx_key, self.input_file_names[idx_key] = %s, %s" % (idx_key, self.input_file_names)
            input_file_name  = os.path.join(self.indir,  self.input_file_names[idx_key] + self.chg_suffix)        
            output_file_name = os.path.join(self.outdir, self.output_file_names[idx_key])        
            dna_region       = self.runobj.samples[idx_key].dna_region
#             print "dna_region = %s" % dna_region
            if dna_region in C.regions_to_chimera_check:
                chimera_region_found = True
            else:
                logger.debug('region not checked: ' +  dna_region)
                continue
             
 
            print "input_file_name = %s \noutput_file_name = %s" % (input_file_name, output_file_name)
 
 
            uchime_cmd = C.clusterize_cmd
            uchime_cmd += " "
            uchime_cmd += self.usearch_cmd
            uchime_cmd += " --uchime "
            uchime_cmd += input_file_name
            uchime_cmd += " --uchimeout "
            uchime_cmd += output_file_name
            uchime_cmd += " --abskew "
            uchime_cmd += self.abskew
             
            print "uchime_cmd = %s" % (uchime_cmd)
             
            try:
                logger.info("chimera denovo command: " + str(uchime_cmd))
#                 subprocess.Popen(uchime_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 
                output[idx_key] = subprocess.Popen(uchime_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#                 print "output[idx_key] = %s" % output[idx_key]
#                 print output[idx_key].split()[2]
#                 cluster_id_list.append(output[idx_key].split()[2])
#                 print 'Have %d bytes in output' % len(output)
#                 print 'denovo', idx_key, output, len(output)
                # len(output) is normally = 47
#                 if len(output[idx_key]) < 50 and len(output[idx_key]) > 40:
#                     logger.debug(idx_key + " uchime denovo seems to have been submitted successfully")
#                 else:
#                     logger.debug("uchime denovo may have broken")  
                 
                # proc.communicate will block - probably not what we want
                #(stdout, stderr) = proc.communicate() #block the last onehere
                #print stderr, stdout
 
#             else:
#                 subprocess.call(uchime_cmd, shell=True)
#                 print uchime_cmd            
#             
#             try:
#                 logger.info("chimera denovo command: " + str(uchime_cmd))
#                 output[idx_key] = subprocess.check_output(uchime_cmd)
#                 print output[idx_key]
#                 print output[idx_key].split()[2]
#                 cluster_id_list.append(output[idx_key].split()[2])
#                 print 'Have %d bytes in output' % len(output)
#                 print 'denovo',idx_key,output,len(output)
#                 # len(output) is normally = 47
#                 if len(output[idx_key]) < 50 and len(output[idx_key]) > 40:
#                     logger.debug(idx_key + " uchime denovo seems to have been submitted successfully")
#                 else:
#                     logger.debug("uchime denovo may have broken")                    
 
            except OSError, e:
                print "Problems with this command: %s" % (uchime_cmd)
                if self.utils.is_local():
                    print >>sys.stderr, "Execution of %s failed: %s" % (uchime_cmd, e)
                else:
                    print >>sys.stderr, "Execution of %s failed: %s" % (uchime_cmd, e)
                    raise                  
                                
# ???
        if not chimera_region_found:            
            return ('NOREGION', 'No regions found that need checking', '')
         
        # ???
#         for idx_key in output:
#             if len(output[idx_key]) > 50 or len(output[idx_key]) < 40:
#                 return ('ERROR','uchime ref may have broken or empty', idx_key)  
         
        # finally
        if cluster_id_list: 
            return ('SUCCESS', 'uchime ref seems to have been submitted successfully', cluster_id_list)
         
    def chimera_reference(self):
     
        chimera_region_found = False
        output = {}
        cluster_id_list = []
        for idx_key in self.run_keys:
             
            dna_region  = self.runobj.samples[idx_key].dna_region
            if dna_region in C.regions_to_chimera_check:
                chimera_region_found = True
            else:
                logger.debug('region not checked: ' + dna_region)                    
                continue
             
            out_file_name = self.prefix[idx_key] + ".chimeras.db"      
             
            # which ref db to use?
            ref_db = ''
            if dna_region.upper() == 'ITS':
                logger.debug("got an ITS dna region so using refdb: " + self.its_refdb)
                ref_db = self.its_refdb
            else:
                logger.debug("using standard refdb: " + self.refdb)
                ref_db = self.refdb
                 
            uchime_cmd = ["clusterize"]
            uchime_cmd.append(self.usearch_cmd)
            uchime_cmd.append("--uchime")
            uchime_cmd.append(self.files[idx_key]['abund'])
            uchime_cmd.append("--uchimeout")
            uchime_cmd.append(out_file_name)
            uchime_cmd.append("--db")
            uchime_cmd.append(ref_db)
                          
            try:
                logger.info("chimera reference command: " + str(uchime_cmd))
                output[idx_key] = subprocess.check_output(uchime_cmd)
                #print 'outsplit',output[idx_key].split()[2]
                cluster_id_list.append(output[idx_key].split()[2])
                #print 'Have %d bytes in output' % len(output)
                #print 'ref',idx_key,output,len(output)
                if len(output[idx_key]) < 50 and len(output[idx_key]) > 40:
                    logger.debug(idx_key + " uchime ref seems to have been submitted successfully")                    
                else:
                    print >>sys.stderr, "uchime ref may be broke"
                
            except OSError, e:
                print >>sys.stderr, "Execution of chimera_reference failed: %s" % (uchime_cmd, e)
                raise
 
        if not chimera_region_found:            
            return ('NOREGION','No regions found that need checking','')
               
        for idx_key in output:
            if len(output[idx_key]) > 50 or len(output[idx_key]) < 40:
                return ('ERROR','uchime ref may have broken or empty',idx_key)  
         
        return ('SUCCESS','uchime ref seems to have been submitted successfully',cluster_id_list)
        
            
    def write_chimeras_to_deleted_file(self): 
    
        for idx_key in self.run_keys:
            # open  deleted file and append chimera to it
            # open and read both chimeras files: chimeras.db and chimeras.txt
            
            # hash to remove dupes
            chimera_deleted = {}
            for file in [self.files[idx_key]['chimera_db'], self.files[idx_key]['chimera_txt']]:            
                fh = open(file,"r") 
                # make a list of chimera deleted read_ids            
                for line in fh.readlines():
                    lst = line.strip().split()
                    id = lst[1].split(';')[0]
                    chimera_yesno = lst[-1]
                    if(chimera_yesno) == 'Y':
                        chimera_deleted[id] = 'chimera'
        
            fh_del = open(self.files[idx_key]['deleted'],"a")
            for id in chimera_deleted:
                fh_del.write(id+"\tchimera\n") 
            
