import subprocess
import sys, os
import re
import time
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils
from pprint import pprint
from collections import defaultdict
sys.path.append("/xraid/bioware/linux/seqinfo/bin")
sys.path.append("/Users/ashipunova/bin/illumina-utils")
sys.path.append("/Users/ashipunova/bin/illumina-utils/illumina-utils/scripts")
sys.path.append("/bioware/merens-illumina-utils")

import fastalib as fa
import pipeline.constants as C

class Chimera:
    """ Define here """
    def __init__(self, runobj = None):
        self.utils      = PipelneUtils()
        self.runobj     = runobj
        self.run_keys   = self.runobj.run_keys
        self.rundate    = self.runobj.run
        
        self.chg_suffix         = ".chg"
        self.chimeras_suffix    = ".chimeras"      
        self.ref_suffix         = ".db"      
        self.denovo_suffix      = ".txt"        
        self.nonchimeric_suffix = ".nonchimeric.fa"
        self.chimeric_suffix    = ".chimeric.fa"
        self.base_suffix        = "unique" + self.chimeras_suffix


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
#         self.usearch_cmd = C.usearch_cmd
        self.usearch_cmd = C.usearch6_cmd        
#         self.abskew      = C.chimera_checking_abskew
        self.refdb       = C.chimera_checking_refdb_6
        self.its_refdb   = C.chimera_checking_its_refdb_6
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
        return filename.endswith((self.chimeras_suffix + self.denovo_suffix, self.chimeras_suffix + self.ref_suffix, self.chimeric_suffix, self.nonchimeric_suffix))

    def get_current_filenames(self, cur_dirname):
        cur_file_names = []
        if cur_dirname == self.indir:
            cur_file_names = self.input_file_names.values()
        elif cur_dirname == self.outdir:
            cur_file_names = self.get_chimera_file_names(self.outdir)
        return cur_file_names

    def get_chimera_file_names(self, cur_dirname):
        cur_file_names = []        
        for dirname, dirnames, filenames in os.walk(cur_dirname):
            for filename in filenames:
                if (self.is_chimera_check_file(filename)):
#                     print "filename = %s" % filename
                    cur_file_names.append(filename)
        return cur_file_names

#     def illumina_frequency_size(self, in_or_out = "", find = "frequency:", replace = ";size="):
#         cur_dirname    = self.get_current_dirname(in_or_out)
#         cur_file_names = self.get_current_filenames(cur_dirname)
# #         print "cur_file_names: "
# #         pprint(cur_file_names)
#         change_from_suffix = ""
#         change_to_suffix   = self.chg_suffix
# #         print "find = %s, replace = %s" % (find, replace)
#         regex              = re.compile(r"%s" % find)
# 
#         for cur_file_name in cur_file_names:
#             file_name = os.path.join(cur_dirname, cur_file_name)
#             with open(file_name + change_from_suffix, "r") as sources:
#                 lines = sources.readlines()
#             with open(file_name + change_to_suffix, "w") as target:
#                 for line in lines:
#                         target.write(regex.sub(replace, line))

    def read_file(self, source_name):
        with open(source_name, "r") as sources:
            return sources.readlines()

    def illumina_sed(self, cur_file_name, lines, target_name, regex, replace):
        with open(target_name, "w") as target:
            for line in lines:
                if line.startswith(">"):
                    line1 = regex.sub(replace, line)
                else:
                    line1 = line.upper()
                target.write(line1)  


    def call_illumina_sed(self):
        find    = "frequency:"
        replace = ";size="
        regex   = re.compile(r"%s" % find)        
        
        cur_dirname        = self.indir
        cur_file_names     = self.get_current_filenames(cur_dirname)
        change_from_suffix = ""
        change_to_suffix   = self.chg_suffix
#         print "find = %s, replace = %s" % (find, replace)
 
        for cur_file_name in cur_file_names:
            source_name = cur_file_name + change_from_suffix
            target_name  = cur_file_name + change_to_suffix 
            lines = self.read_file(source_name)
            self.illumina_sed(self, cur_file_name, lines, target_name, regex, replace)

    def illumina_freq_to_size_in_chg(self):
#         TODO: refactor
        find1    = "frequency:"
        replace1 = ";size="
        regex1   = re.compile(r"%s" % find1)        
        
#         print "cur_file_names: "
#         pprint(cur_file_names)
        cur_dirname        = self.get_current_dirname()
        cur_file_names     = self.get_current_filenames(cur_dirname)
        change_from_suffix = ""
        change_to_suffix   = self.chg_suffix
#         print "find = %s, replace = %s" % (find, replace)
 
        for cur_file_name in cur_file_names:
            file_name = os.path.join(cur_dirname, cur_file_name)
            with open(file_name + change_from_suffix, "r") as sources:
                lines = sources.readlines()
            with open(file_name + change_to_suffix, "w") as target:
                for line in lines:
                    if line.startswith(">"):
                        line1 = regex1.sub(replace1, line)
                    else:
                        line1 = line.upper()
                    target.write(line1)  


    def illumina_size_to_freq_in_chimer(self):
        find1           = ";size="
        replace1        = "frequency:"
        regex1          = re.compile(r"%s" % find1)        
 
        cur_file_names = self.get_chimera_file_names(self.outdir)
                    
        for file_chim in cur_file_names:
            file_chim_path = os.path.join(self.outdir, file_chim)
            with open(file_chim_path, "r") as sources:
                lines = sources.readlines()
            with open(file_chim_path, "w") as target:
                for line in lines:
                    line1 = regex1.sub(replace1, line)
                    target.write(line1)                    
              
    def illumina_rm_size_files(self):
        for idx_key in self.input_file_names:
            file_name = os.path.join(self.indir, self.input_file_names[idx_key] + self.chg_suffix)
            if os.path.exists(file_name):
                os.remove(file_name)
    
#     def illumina_chimera_size_files(self):
#     
#     import os
# [os.rename(f, f.replace('_', '-')) for f in os.listdir('.') if not f.startswith('.')]
    def get_time_now(self):
        """date and hour only!"""
        return time.strftime("%m/%d/%Y %H", time.localtime())
# '2009-01-05 22'
        
          
    def check_if_cluster_is_done(self, time_before):
        cluster_done = False
        check_qstat_cmd_line = "qstat | grep \"%s\" | grep usearch | wc -l" % time_before
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
        """
        http://www.drive5.com/usearch/manual/uchime_denovo.html
        from usearch -help
        Chimera detection (UCHIME ref. db. mode):
          usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
            [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]
         
        Chimera detection (UCHIME de novo mode):
          usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
             [-uchimeout results.uch] [-uchimealns results.alns]
          Input is estimated amplicons with integer abundances specified using ";size=N".
        usearch -uchime_denovo amplicons.fasta -uchimeout results.uchime
        """        

        uchime_cmd_append = ""
        db_cmd_append     = ""
        dir_cmd_append    = ""

        if (ref_or_novo == "denovo"):
            uchime_cmd_append = " -uchime_denovo "           
            output_file_name  = output_file_name + self.chimeras_suffix + self.denovo_suffix 
        elif (ref_or_novo == "ref"):
            uchime_cmd_append = " -uchime_ref "
            output_file_name  = output_file_name + self.chimeras_suffix + self.ref_suffix           
            db_cmd_append     = " -db " + ref_db   
            dir_cmd_append    = " -strand plus"
        else:
            print "Incorrect method, should be \"denovo\" or \"ref\"" 
        print "output_file_name = %s" % output_file_name 


        uchime_cmd = C.clusterize_cmd
        uchime_cmd += " "
        uchime_cmd += self.usearch_cmd
        uchime_cmd += uchime_cmd_append + input_file_name
        uchime_cmd += db_cmd_append
        uchime_cmd += " -uchimeout " + output_file_name
        """if we need nonchimeric for denovo and db separate we might create them here
#         uchime_cmd += " -nonchimeras "
#         uchime_cmd += (output_file_name + self.nonchimeric_suffix)
"""
        uchime_cmd += " -chimeras " + (output_file_name + self.chimeric_suffix)         
        uchime_cmd += dir_cmd_append
        uchime_cmd += " -notrunclabels"
        
        
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
        
    def get_chimeric_ids(self):
        ids = set()
        chimera_file_names = self.get_chimera_file_names(self.outdir)
        file_ratio = self.check_chimeric_stats()
        
        for file_name in chimera_file_names:
#             print "from get_chimeric_ids: file_name = %s" % file_name
            if file_name.endswith(self.chimeric_suffix):
                both_or_denovo = self.get_chimeras_suffix(file_ratio, file_name)
#                 TODO: run ones for each file_base = ".".join(file_name.split(".")[0:3]) (for txt and db)
                if file_name.endswith(both_or_denovo):                    
                    file_name_path = os.path.join(self.outdir, file_name)        
                    print "Get ids from %s" % file_name_path
                    read_fasta     = fa.ReadFasta(file_name_path)
                    ids.update(set(read_fasta.ids))
        return ids
    
    
    def get_chimeras_suffix(self, file_ratio, file_name):
        """ use only de-novo (.txt) chimeric if
            check_chimeric_stats shows
            ref > 15% and ratio ref to de-novo > 2
            e.g.
            if denovo_only:
                chimeric_suffix = self.chimeras_suffix + self.denovo_suffix + self.chimeric_suffix
            if no: 
                chimeras_suffix = self.chimeric_suffix
                
            if file_name.endswith(chimeric_suffix):
            ...        
                #     first_name, last_name = get_name()

        """         
#         for file_basename in file_ratio:
        (percent_ref, ratio) = file_ratio[".".join(file_name.split(".")[0:3])]   

        chimeric_fa_suffix = ""
#         print "percent_ref = %s, ratio = %s" % (percent_ref, ratio)
        if (percent_ref > 15) and (ratio > 2):
            chimeric_fa_suffix = self.chimeras_suffix + self.denovo_suffix + self.chimeric_suffix
        else: 
            chimeric_fa_suffix = self.chimeric_suffix   
        return chimeric_fa_suffix 
    
    def move_out_chimeric(self):
        chimeric_ids = self.get_chimeric_ids()
        for idx_key in self.input_file_names:
            fasta_file_path    = os.path.join(self.indir, self.input_file_names[idx_key])   
            read_fasta         = fa.ReadFasta(fasta_file_path)
            read_fasta.close()
            
            non_chimeric_file  = fasta_file_path + self.nonchimeric_suffix
            non_chimeric_fasta = fa.FastaOutput(non_chimeric_file)

            fasta              = fa.SequenceSource(fasta_file_path, lazy_init = False) 
            while fasta.next():
                if not fasta.id in chimeric_ids:
                    non_chimeric_fasta.store(fasta, store_frequencies = False)
            non_chimeric_fasta.close()


    def check_chimeric_stats(self):
        all_lines_suffix      = self.denovo_suffix # ".txt" or ".db, doesn't matter"
        chimera_ref_suffix    = self.ref_suffix + self.chimeric_suffix #".db.chimeric.fa"
        chimera_denovo_suffix = self.denovo_suffix + self.chimeric_suffix # ".txt.chimeric.fa"
        filenames             = self.get_basenames(self.get_current_filenames(self.outdir))
        file_ratio            = {}
        for file_basename in filenames:
            # print file_basename
            all_lines      = 0
            ref_lines      = 0
            denovo_lines   = 0
            ratio          = 0
            percent_ref    = 0
            percent_denovo = 0
        
            all_lines_file_name    = os.path.join(self.outdir, file_basename + all_lines_suffix)
            ref_lines_file_name    = os.path.join(self.outdir, file_basename + chimera_ref_suffix)
            denovo_lines_file_name = os.path.join(self.outdir, file_basename + chimera_denovo_suffix)
        
            all_lines    = int(self.wccount(all_lines_file_name) or 0)
            ref_lines    = int(self.get_fa_lines_count(ref_lines_file_name) or 0)
            denovo_lines = int(self.get_fa_lines_count(denovo_lines_file_name) or 0)
        
            # denovo_lines = int(denovo_lines or 0)
            if (ref_lines == 0) or (all_lines == 0):
                file_ratio[file_basename] = (0, 0)
                continue
            else:
                percent_ref = self.percent_count(all_lines, ref_lines)
                
            if (denovo_lines == 0):
                file_ratio[file_basename] = (percent_ref, percent_ref) #use ref instead of ratio, because we are actually looking for a huge difference between ref and denovo (ref > 15 and denovo = 0)
                continue
        
            if (denovo_lines > 0):            
                ratio          = self.count_ratio(ref_lines, denovo_lines)        
                percent_denovo = self.percent_count(all_lines, denovo_lines)
            file_ratio[file_basename] = (percent_ref, ratio)
            # percent_ref = int(percent_ref or 0)
            if (percent_ref > 1):
                print "=" * 50
            
                print file_basename
                # print "all_lines_file_name = %s, ref_lines_file_name = %s, denovo_lines_file_name = %s" % (all_lines_file_name, ref_lines_file_name, denovo_lines_file_name)
                print "all_lines = %s, ref_lines = %s, denovo_lines = %s" % (all_lines, ref_lines, denovo_lines)
                print "ratio = %s" % ratio 
                print "percent_ref = %s, percent_denovo = %s" % (percent_ref, percent_denovo)
        return file_ratio
        
        
    def get_basenames(self, filenames):
        file_basenames = set()
        for f in filenames:
            file_basename = ".".join(f.split(".")[0:3])
            if file_basename.endswith(self.base_suffix):
                file_basenames.add(file_basename)
    
        return file_basenames
    
    def wccount(self, filename):
        return subprocess.check_output(['wc', '-l', filename]).split()[0]
    
    def count_ratio(self, ref_num, denovo_num):
        try:
            return float(ref_num or 0) / float(denovo_num or 0) 
        except ZeroDivisionError:
            # print "There is no denovo chimeras to count ratio."
            pass
        else:
            raise
    
    def get_fa_lines_count(self, file_name):
        # return fa.SequenceSource(file_name, lazy_init = False).total_seq
        try:
            file_open = open(file_name)
            return len([l for l in file_open.readlines() if l.startswith('>')])
        except IOError, e:
            print e
            return 0
            # print "%s\nThere is no such file: %s" % (e, file_name)
        else:
            raise
    
    def percent_count(self, all_lines, chimeric_count):
        try:
            return float(chimeric_count or 0) * 100 / float(all_lines or 0)
        except ZeroDivisionError:
            # print "There is no denovo chimeras to count ratio."
            pass
        else:
            raise
        
    
 
    """ 
    -----------------------------------------------------------------------------
        For 454.
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