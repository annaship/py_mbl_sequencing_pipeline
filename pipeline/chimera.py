import subprocess
import sys, os
import re
import time
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs

import pipeline.constants as C

class Chimera:
    """ Define here """
    def __init__(self, runobj = None):
        self.runobj   = runobj
        self.run_keys = self.runobj.run_keys
        self.rundate  = self.runobj.run

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
        #         self.outdir    = run.output_dir
        self.indir  = dirs.check_dir(dirs.reads_overlap_dir)
        self.outdir = dirs.check_dir(dirs.chimera_dir)

        
#         self.usearch_cmd = 'usearch'
        self.usearch_cmd = C.usearch_cmd
        self.abskew      = C.chimera_checking_abskew
        self.refdb       = C.chimera_checking_refdb
        self.its_refdb   = C.chimera_checking_its_refdb
        self.input_file_names  = self.make_chimera_input_illumina_file_names()
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
            output_file_names[idx_key] = input_file_name + ".chimeras.txt"
        return output_file_names

#             
#             allTheLists = [[] for x in range(int(L))]
# #...
# allTheLists[j].append(jListItem)
# printAllTheLiats[listIndex][itemIndex]
            
    def illumina_frequncy_size(self):
        for idx_key in self.input_file_names:
            file_name = os.path.join(self.indir, self.input_file_names[idx_key])
            with open(file_name, "r") as sources:
                lines = sources.readlines()
            with open(file_name + ".size", "w") as target:
                for line in lines:
                    target.write(re.sub(r'frequency:', ';size=', line))
          
       
    def chimera_denovo(self):
        
        chimera_region_found = False
        output = {}
        cluster_id_list = []
        for idx_key in self.input_file_names:
            file_name  = os.path.join(self.indir, self.input_file_names[idx_key])        
            dna_region = self.runobj.samples[idx_key].dna_region
            if dna_region in C.regions_to_chimera_check:
                chimera_region_found = True
            else:
                logger.debug('region not checked: ' +  dna_region)
                continue
            

            print file_name

            uchime_cmd = "clusterize "
            uchime_cmd += self.usearch_cmd
            uchime_cmd += " --uchime "
            uchime_cmd += file_name + ".size"
            uchime_cmd += " --uchimeout "
            uchime_cmd += self.output_file_names[idx_key]
            uchime_cmd += " --abskew "
            uchime_cmd += self.abskew
            
            print "uchime_cmd = %s" % (uchime_cmd)
            try:
                logger.info("chimera denovo command: " + str(uchime_cmd))
                output[idx_key] = subprocess.check_output(uchime_cmd)
                #print output[idx_key]
                #print output[idx_key].split()[2]
                cluster_id_list.append(output[idx_key].split()[2])
                #print 'Have %d bytes in output' % len(output)
                #print 'denovo',idx_key,output,len(output)
                # len(output) is normally = 47
                if len(output[idx_key]) < 50 and len(output[idx_key]) > 40:
                    logger.debug(idx_key + " uchime denovo seems to have been submitted successfully")
                else:
                    logger.debug("uchime denovo may have broken")                    

            except OSError, e:
                print >>sys.stderr, "Execution failed:", e               
        
        if not chimera_region_found:            
            return ('NOREGION','No regions found that need checking','')
        
        for idx_key in output:
            if len(output[idx_key]) > 50 or len(output[idx_key]) < 40:
                return ('ERROR','uchime ref may have broken or empty',idx_key)  
        
        # finally
        return ('SUCCESS','uchime ref seems to have been submitted successfully',cluster_id_list)
        
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
                logger.info("chimera referenc command: " + str(uchime_cmd))
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
                print >>sys.stderr, "Execution failed:", e 
        
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
            
