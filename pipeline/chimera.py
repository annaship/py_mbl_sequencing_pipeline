import subprocess
import sys, os
import re
import time
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils
from pprint import pprint

import pipeline.constants as C

class Chimera:
    """ Define here """
    def __init__(self, runobj = None):
        self.runobj     = runobj
        self.run_keys   = self.runobj.run_keys
        self.rundate    = self.runobj.run
        self.chg_suffix = ".chg"
        self.utils      = PipelneUtils()

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
        
        dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, lane_name = lane_name, site = site) 
        self.indir  = dirs.check_dir(dirs.reads_overlap_dir)
        self.outdir = dirs.check_dir(dirs.chimera_dir)
        self.usearch_cmd = C.usearch_cmd
        self.abskew      = C.chimera_checking_abskew
        self.refdb       = C.chimera_checking_refdb
        self.its_refdb   = C.chimera_checking_its_refdb
        self.input_file_names  = self.make_chimera_input_illumina_file_names()
#         pprint(self.run_keys)
        self.output_file_names = self.make_chimera_output_illumina_file_names(self.input_file_names)

    def make_chimera_input_illumina_file_names(self):
        input_file_names = {} 
        
        for idx_key in self.run_keys:
            file_name = idx_key + "_" + C.filtered_suffix + ".unique" 
           
            if os.path.exists(os.path.join(self.indir, file_name)):
                input_file_names[idx_key] = file_name
        
        return input_file_names
            
    def make_chimera_output_illumina_file_names(self, input_file_names):
        output_file_names = {} 
        for idx_key, input_file_name in input_file_names.iteritems():
            output_file_names[idx_key] = input_file_name + ".chimeras"
        return output_file_names

    def illumina_frequency_size(self, cur_dirname = "", find = "frequency:", replace = ";size="):
        if cur_dirname == "":
            cur_dirname    = self.indir 
            cur_file_names = self.input_file_names
        else:
            cur_dirname    = self.outdir
            cur_file_names = self.output_file_names

        change_from_suffix = ""
        change_to_suffix   = self.chg_suffix
        print "find = %s, replace = %s" % (find, replace)
        regex = re.compile(r"%s" % find)
            
        for idx_key in cur_file_names:
            file_name = os.path.join(cur_dirname, cur_file_names[idx_key])
            with open(file_name + change_from_suffix, "r") as sources:
                lines = sources.readlines()
            with open(file_name + change_to_suffix, "w") as target:
                for line in lines:
                    target.write(regex.sub(replace, line))

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
        
        p = subprocess.Popen(check_qstat_cmd_line, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        print "qstat is running %s 'usearch' processes", output
        pprint(p)
        
        if (int(output) == 0):
            cluster_done = True
        print "cluster_done from check_if_cluster_is_done = %s" % cluster_done
        return cluster_done
          
    def create_chimera_cmd(self, input_file_name, output_file_name, ref_or_novo, ref_db = ""):
        uchime_cmd = C.clusterize_cmd
        uchime_cmd += " "
        uchime_cmd += self.usearch_cmd
        uchime_cmd += " --uchime "
        uchime_cmd += input_file_name
 
        if (ref_or_novo == "denovo"):
            output_file_name = output_file_name + ".txt"
            cmd_append = " --abskew " + self.abskew                                   
            
        elif (ref_or_novo == "ref"):
            output_file_name = output_file_name + ".db"                        
            cmd_append = " --db " + ref_db                                   

        uchime_cmd += " --uchimeout "
        uchime_cmd += output_file_name
        uchime_cmd += cmd_append
        print "uchime_cmd FROM create_chimera_cmd = %s" % (uchime_cmd)
        return uchime_cmd
        
    def get_ref_db(self, dna_region):
        ref_db = ''
        if dna_region.upper() == 'ITS':
            logger.debug("got an ITS dna region so using refdb: " + self.its_refdb)
            ref_db = self.its_refdb
        else:
            logger.debug("using standard refdb: " + self.refdb)
            ref_db = self.refdb
        print  "ref_db = %s" % ref_db  
        return ref_db       
    
    def chimera_checking(self, ref_or_novo):
        chimera_region_found = False
        output = {}
        
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
            
#             print "input_file_name = %s \noutput_file_name = %s" % (input_file_name, output_file_name)
            ref_db     = self.get_ref_db(dna_region)
            print "dna_region = %s; ref_db = %s; ref_or_novo = %s" % (dna_region, ref_db, ref_or_novo)
            
            uchime_cmd = self.create_chimera_cmd(input_file_name, output_file_name, ref_or_novo, ref_db)
            print "uchime_cmd = %s" % (uchime_cmd)
            
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
                
            uchime_cmd.append("--db")
            uchime_cmd.append(ref_db)                
                
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
            
