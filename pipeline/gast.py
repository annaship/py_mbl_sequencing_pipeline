import subprocess
import sys, os, stat
import time
import shutil
from pipeline.pipelinelogging import logger
#import logging
import constants as C
import re
import json  
from types import *
from pipeline.utils import Dirs, PipelneUtils


class Gast:
    """Doc string here.."""
    Name = "GAST"
    def __init__(self, run_object = None, idx_keys=[]):                       
        self.runobj = run_object       
        self.test   = True
        self.utils  = PipelneUtils()
        
        # Create and move all from below
        # self.set_directories()
       
        dirs = Dirs(self.runobj.vamps_user_upload, self.runobj.run, self.runobj.platform) 
#
#        self.out_file_path = dirs.check_dir(dirs.analysis_dir)
#        self.results_path  = dirs.check_dir(dirs.reads_overlap_dir)
        
        self.basedir = self.runobj.output_dir
        self.reads_dir = dirs.check_dir(dirs.reads_overlap_dir)
#        if self.runobj.platform == 'illumina' and not self.runobj.vamps_user_upload:
#            self.basedir = os.path.join(self.runobj.input_dir, self.runobj.run)
        
        self.use_cluster = self.runobj.use_cluster
        if self.runobj.vamps_user_upload:        
            self.idx_keys  = [self.runobj.user+self.runobj.run]
            self.refdb_dir = C.vamps_ref_database_dir
            self.iterator  = self.runobj.datasets
        else:
            self.idx_keys  = idx_keys
            self.iterator  = self.idx_keys
            self.refdb_dir = C.ref_database_dir

            if self.utils.is_local():
#                program_name = "/Users/ashipunova/bin/illumina-utils/analyze-illumina-v6-overlaps"        
                self.refdb_dir = "/Users/ashipunova/bin/illumina-utils/"

        os.environ['SGE_ROOT'] ='/usr/local/sge'
        os.environ['SGE_CELL'] ='grendel'
        path                   = os.environ['PATH']
        os.environ['PATH']     = '/usr/local/sge/bin/lx24-amd64:'+path
        
        # for testing
        self.limit = 400
       
        # If we are here from a vamps gast process
        # then there should be just one dataset to gast
        # but if MBL/illumina pipe then many datasets are probably involved.
        self.analysis_dir = dirs.check_dir(dirs.analysis_dir)
#        if not os.path.exists(self.analysis_dir):
#            os.mkdir(self.analysis_dir)

        self.global_gast_dir = dirs.check_dir(dirs.gast_dir)
#        if self.runobj.vamps_user_upload:
#            self.global_gast_dir = self.basedir
#        else:
#            self.global_gast_dir = os.path.join(self.basedir, C.gast_dir)
#            if os.path.exists(self.global_gast_dir):
#                # delete gast directory and recreate
#                shutil.rmtree(self.global_gast_dir, True)
#                os.mkdir(self.global_gast_dir)
#            else:
#                os.mkdir(self.global_gast_dir)
                
#        if self.runobj.platform == 'illumina' and not self.runobj.vamps_user_upload:
#            reads_dir = dirs.check_dir(dirs.reads_overlap_dir)
#            if os.path.exists(reads_dir):
#                self.input_dir = reads_dir
#            else:
#                self.input_dir = self.runobj.input_dir
            
        # create our directories for each key
        dirs.create_gast_name_dirs(self.iterator)
#        for key in self.iterator:
#            output_dir = os.path.join(self.global_gast_dir, key)
#            if output_dir and not os.path.exists(output_dir):
#                os.mkdir(output_dir)
#                 gast_dir = output_dir
#                 if os.path.exists(gast_dir):
#                     # empty then recreate directory
#                     shutil.rmtree(gast_dir)
#                     os.mkdir(gast_dir)
#                 else:
#                     os.mkdir(gast_dir)
        
    def clustergast(self):
        """
        clustergast - runs the GAST pipeline on the cluster.
               GAST uses UClust to identify the best matches of a read sequence
               to references sequences in a reference database.
               VAMPS: The uniques and names files have previously been created in trim_run.py.
               Illumina: The uniques and names files have been created by illumina_files.py
        """
        logger.info("Starting Clustergast")
        
        self.runobj.run_status_file_h.write(json.dumps({'status':'STARTING_CLUSTERGAST'})+"\n")
        # Step1: create empty gast table in database: gast_<rundate>
        # Step2: Count the number of sequences so the job can be split for nodes
        # $facount = `grep -c \">\" $fasta_uniqs_filename`;
        # $calcs = `/bioware/seqinfo/bin/calcnodes -t $facount -n $nodes -f 1`;

        #   /bioware/seqinfo/bin/fastasampler -n $start, $end ${gastDir}/${fasta_uniqs_filename} $tmp_fasta_filename
        #   $usearch_binary --global --query $tmp_fasta_filename --iddef 3 --gapopen 6I/1E --db $refhvr_fa --uc $tmp_usearch_filename --maxaccepts $max_accepts --maxrejects $max_rejects --id $pctid_threshold
        #   # sort the results for valid hits, saving only the ids and pct identity
        #   grep -P \"^H\\t\" $tmp_usearch_filename | sed -e 's/|.*\$//' | awk '{print \$9 \"\\t\" \$4 \"\\t\" \$10 \"\\t\" \$8}' | sort -k1,1b -k2,2gr | clustergast_tophit > $gast_filename
        #   Submit the script
        #   /usr/local/sge/bin/lx24-amd64/qsub $qsub_priority $script_filename
        
        calcnodes = C.calcnodes_cmd
        if self.utils.is_local():
            calcnodes = "/Users/ashipunova/bin/illumina-utils/calcnodes"           
        
        sqlImportCommand = C.mysqlimport_cmd
        if self.utils.is_local():
            sqlImportCommand = "/usr/local/mysql/bin/mysqlimport"           
        
        #qsub = '/usr/local/sge/bin/lx24-amd64/qsub'
        clusterize = C.clusterize_cmd

        ###################################################################
        # use fasta.uniques file
        # split into smaller files
        # usearch --cluster each
        #######################################
        #
        # Split the uniques fasta and run UClust per node
        #
        #######################################
        qsub_prefix = 'clustergast_sub_'
        gast_prefix = 'gast_'
        if self.use_cluster:
            logger.info("Using cluster for clustergast")
        else:
            logger.info("Not using cluster")
            
            
        counter=0
        gast_file_list = []
        cluster_nodes = C.cluster_nodes
        logger.info("Cluster nodes set to: "+str(cluster_nodes))
        for key in self.iterator:
            counter += 1
            if self.runobj.vamps_user_upload:
                output_dir  = os.path.join(self.global_gast_dir, key)
                gast_dir    = os.path.join(self.global_gast_dir, key)
                fasta_file  = os.path.join(output_dir, 'fasta.fa')
                unique_file = os.path.join(output_dir, 'unique.fa')
                names_file  = os.path.join(output_dir, 'names')
                #datasets_file = os.path.join(self.global_gast_dir, 'datasets')
                print 'gast_dir:', gast_dir
                print 'unique_file:', unique_file
            else:
                if self.runobj.platform == 'illumina':
                    output_dir  = os.path.join(self.global_gast_dir, key)
                    gast_dir    = os.path.join(self.global_gast_dir, key)
                    file_prefix = key
#                    file_prefix = self.runobj.samples[key].file_prefix
                    unique_file = os.path.join(self.reads_dir, file_prefix + "-PERFECT_reads.fa.unique")
                    names_file  = os.path.join(self.reads_dir, file_prefix + "-PERFECT_reads.fa.unique.names")
                elif self.runobj.platform == '454':
                    pass
                else:
                    sys.exit("clustergast: no platform")
            
            if counter >= self.limit:
                pass
                   
            #print 'samples', key, self.runobj.samples
            if key in self.runobj.samples:
                dna_region = self.runobj.samples[key].dna_region
            else:            
                dna_region = self.runobj.dna_region
            if not dna_region:
                logger.error("clustergast: We have no DNA Region: Setting dna_region to 'unknown'")
                dna_region = 'unknown'
                
            (refdb, taxdb) = self.get_reference_databases(dna_region)
            
            
            print "\nFile", str(counter), key
            print 'use_cluster:', self.use_cluster
            if os.path.exists(unique_file) and (os.path.getsize(unique_file) > 0):
                print "cluster nodes: "+str(cluster_nodes)
                i = 0
                if cluster_nodes:
                    grep_cmd = ['grep', '-c', '>', unique_file]
                    logger.debug( ' '.join(grep_cmd) )
                    facount = subprocess.check_output(grep_cmd).strip()
                    logger.debug( key+' count '+facount)
                    calcnode_cmd = [calcnodes, '-t', str(facount), '-n', str(cluster_nodes), '-f', '1']
                    
                    calcout = subprocess.check_output(calcnode_cmd).strip()
                    logger.debug("calcout:\n"+calcout)
                    print "facount: "+str(facount)
                    #calcout:
                    # node=1 start=1 end=1 rows=1
                    # node=2 start=2 end=2 rows=1
                    # node=3 start=3 end=3 rows=1           
                    lines = calcout.split("\n")
                    #gast_file_list = []
                    for line in lines:
                        i += 1
                        if i >= cluster_nodes:
                            continue
                        script_filename      = os.path.join(gast_dir, qsub_prefix + str(i))
                        gast_filename        = os.path.join(gast_dir, gast_prefix + str(i))
                        fastasamp_filename   = os.path.join(gast_dir, 'samp_' + str(i))
                        clustergast_filename = os.path.join(gast_dir, key+".gast_" + str(i))
                        gast_file_list.append(clustergast_filename)
                        usearch_filename     = os.path.join(gast_dir, "uc_" + str(i))
                        log_file             = os.path.join(gast_dir, 'clustergast.log_' + str(i))
                        
                        data = line.split()
                        
                        if len(data) < 2:
                            continue
                        start = data[1].split('=')[1]
                        end  = data[2].split('=')[1]
                        
                        if self.use_cluster:
                            fh = open(script_filename, 'w')
                            qstat_name = "gast" + key + '_' + self.runobj.run + "_" + str(i)
                            fh.write("#!/bin/sh\n\n")
                            
                            # don't need these commands unless running qsub directly
                            #fh.write("#$ -j y\n" )
                            #fh.write("#$ -o " + log_file + "\n")
                            #fh.write("#$ -N " + qstat_name + "\n\n")
        
                            # setup environment
                            fh.write("source /xraid/bioware/Modules/etc/profile.modules\n")
                            fh.write("module load bioware\n\n")
                        
                        fs_cmd = self.get_fastasampler_cmd(unique_file, fastasamp_filename, start, end)
                        
    
                        logger.debug("fastasampler command: "+fs_cmd)
                        
                        if self.use_cluster:
                            fh.write(fs_cmd + "\n")
                        else:
                            subprocess.call(fs_cmd, shell=True)
                        
                        us_cmd = self.get_usearch_cmd(fastasamp_filename, refdb, usearch_filename)
    
                        logger.debug("usearch command: "+us_cmd)
                        #print 'usearch', cmd2
                        if self.use_cluster:
                            fh.write(us_cmd + "\n")
                        else:
                            subprocess.call(us_cmd, shell=True)
                        
                        grep_cmd = self.get_grep_cmd(usearch_filename, clustergast_filename)
    
                        logger.debug("grep command: "+grep_cmd)
                        if self.use_cluster:  
                            fh.write(grep_cmd + "\n")
                            fh.close()
                            
                            # make script executable and run it
                            os.chmod(script_filename, stat.S_IRWXU)
                            qsub_cmd = clusterize + " " + script_filename
                            
                            # on vamps and vampsdev qsub cannot be run - unless you call it from the
                            # cluster aware directories /xraid2-2/vampsweb/vamps and /xraid2-2/vampsweb/vampsdev
                            #qsub_cmd = C.qsub_cmd + " " + script_filename
                            logger.debug("qsub command: "+qsub_cmd)
                            print "qsub command: "+qsub_cmd
                            #subprocess.call(qsub_cmd, shell=True)
                            proc = subprocess.Popen(qsub_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            # proc.communicate will block - probably not what we want
                            #(stdout, stderr) = proc.communicate() #block the last onehere
                            #print stderr, stdout
        
                        else:
                            subprocess.call(grep_cmd, shell=True)
                            print grep_cmd
                
                else:
                    #fastasamp_filename = os.path.join(gast_dir, 'samp')
                    # no nodes means that just one file will be run by clusterize
                    usearch_filename= os.path.join(gast_dir, "uc")
                    clustergast_filename_single   = os.path.join(gast_dir, "gast"+dna_region)
                    gast_file_list = [clustergast_filename_single]
                    print usearch_filename, clustergast_filename_single
                    
                    us_cmd = self.get_usearch_cmd(unique_file, refdb, usearch_filename)
                    print us_cmd
                    subprocess.call(us_cmd, shell=True)
                    grep_cmd = self.get_grep_cmd(usearch_filename, clustergast_filename_single)
                    print grep_cmd
                    subprocess.call(grep_cmd, shell=True)
            else:
                logger.warning( "unique_file not found or zero size: "+unique_file)
                
        if self.use_cluster:
            # wait here for all the clustergast scripts to finish
            temp_file_list = gast_file_list
        
            c = False
            maxwaittime = C.maxwaittime  # seconds
            sleeptime   = C.sleeptime    # seconds
            counter2 = 0
            while c == False:
                counter2 += 1
                if counter2 >= maxwaittime / sleeptime:
                    raise Exception("Max wait time exceeded in gast.py: "+maxwaittime+" seconds")
                for index, file in enumerate(temp_file_list):
                    #print temp_file_list
                    if os.path.exists(file):
                        # remove from tmp list
                        logger.debug("Found file now removing from list: "+file)
                        temp_file_list = temp_file_list[:index] + temp_file_list[index+1:]
                
                if temp_file_list:
                    logger.info("waiting for clustergast files to fill...")
                    logger.debug(' '.join(temp_file_list))
                    logger.info("\ttime: "+str(counter2 * sleeptime)+" | files left: "+str(len(temp_file_list)))
                    time.sleep(sleeptime)
                else:
                    c = True
                  
                  
                  
                  
        for key in self.iterator: 
            if self.runobj.vamps_user_upload:
                gast_dir = os.path.join(self.global_gast_dir, key)
            else:
                if self.runobj.platform == 'illumina':
                    gast_dir = os.path.join(self.global_gast_dir, key)
                elif self.runobj.platform == '454':
                    pass
                else:
                    sys.exit("clustergast: no platform")
      
            # now concatenate all the clustergast_files into one file (if they were split)
            if cluster_nodes:
                # gast file
                clustergast_filename_single   = os.path.join(gast_dir, "gast"+dna_region)
                clustergast_fh = open(clustergast_filename_single, 'w')
                # have to turn off cluster above to be able to 'find' these files for concatenation
                for n in range(1, i-1):
                    #cmd = "cat "+ gast_dir + key+".gast_" + str(n) + " >> " + gast_dir + key+".gast"
                    file = os.path.join(gast_dir, key+".gast_" + str(n))
                    if(os.path.exists(file)):                    
                        shutil.copyfileobj(open(file, 'rb'), clustergast_fh)
                    else:
                        logger.info( "Could not find file: "+os.path.basename(file)+" Skipping")

                clustergast_fh.flush()
                clustergast_fh.close()
            
                
        if not self.test:    
            # remove tmp files
            for n in range(i+1):
                #print "Trying to remove "+os.path.join(gast_dir, "uc_"+str(n))
                if os.path.exists(os.path.join(gast_dir, "uc_"+str(n))):
                    os.remove(os.path.join(gast_dir, "uc_"+str(n)))
                    pass
                #print "Trying to remove "+os.path.join(gast_dir, "samp_"+str(n))
                if os.path.exists(os.path.join(gast_dir, "samp_"+str(n))):    
                    os.remove(os.path.join(gast_dir, "samp_"+str(n)))
                    pass
                #print "Trying to remove "+os.path.join(self.gast_dir, key+".gast_"+str(n))
                if os.path.exists(os.path.join(gast_dir, key+".gast_"+str(n))):    
                    os.remove(os.path.join(gast_dir, key+".gast_"+str(n)))
                    pass
                    
                    
        
        print "Finished clustergast"
        logger.info("Finished clustergast")
        return {'status':"SUCCESS", 'message':"Clustergast Finished"}    
        
            
    def gast_cleanup(self):
        """
        gast_cleanup - follows clustergast, explodes the data and copies to gast_concat and gast files
        """
        self.runobj.run_status_file_h.write(json.dumps({'status':'STARTING_GAST_CLEANUP'})+"\n")
        for key in self.iterator:
            if self.runobj.vamps_user_upload:
                output_dir = os.path.join(self.global_gast_dir, key)
                gast_dir = os.path.join(self.global_gast_dir, key)
                fasta_file = os.path.join(output_dir, 'fasta.fa')
                unique_file = os.path.join(output_dir, 'unique.fa')
                names_file = os.path.join(output_dir, 'names')
                #datasets_file = os.path.join(self.global_gast_dir, 'datasets')
            else:
                if self.runobj.platform == 'illumina':
                    output_dir = os.path.join(self.global_gast_dir, key)
                    gast_dir = os.path.join(self.global_gast_dir, key)
                    file_prefix = self.runobj.samples[key].file_prefix
                    unique_file = os.path.join(self.input_dir, file_prefix+"-PERFECT_reads.fa.unique")
                    names_file = os.path.join(self.input_dir, file_prefix+"-PERFECT_reads.fa.unique.names")
                elif self.runobj.platform == '454':
                    pass
                else:
                    sys.exit("gast_cleanup: no platform")
                
            if key in self.runobj.samples:
                dna_region = self.runobj.samples[key].dna_region
            else:            
                dna_region = self.runobj.dna_region
            if not dna_region:
                logger.error("gast_cleanup: We have no DNA Region: Setting dna_region to 'unknown'")
                self.runobj.run_status_file_h.write(json.dumps({'status':'WARNING', 'message':"gast_cleanup: We have no DNA Region: Setting dna_region to 'unknown'"})+"\n")
                dna_region = 'unknown'
            # find gast_dir
            
            
            if not os.path.exists(gast_dir):
                logger.error("Could not find gast directory: "+gast_dir+" Exiting")
                sys.exit()
                
            clustergast_filename_single   = os.path.join(gast_dir, "gast"+dna_region)
            try:
                logger.debug('gast filesize:'+str(os.path.getsize(clustergast_filename_single)))
            except:
                logger.debug('gast filesize: zero')
                
            gast_filename          = os.path.join(gast_dir, "gast")
            gastconcat_filename    = os.path.join(gast_dir, "gast_concat")  
            #dupes_filename    = os.path.join(gast_dir, "dupes") 
            #nonhits_filename    = os.path.join(gast_dir, "nonhits")   
            copies = {}
            nonhits = {}
            # open and read names file
            
            if os.path.exists(names_file) and os.path.getsize(names_file) > 0:
                names_fh = open(names_file, 'r')
                for line in names_fh:
                    s = line.strip().split("\t")
                    
                    index_read = s[0]                
                    copies[index_read] = s[1].split(',')
                    
                    if index_read in nonhits:
                        nonhits[index_read] += 1
                    else:
                        nonhits[index_read] = 1
                        
                    
                    
                names_fh.close()            
            #print nonhits
            #print copies
            
            #######################################
            # 
            #  Insert records with valid gast hits into gast_file
            # 
            #######################################   
            # read the .gast file from clustergast            
            concat = {}
            
            gast_fh     = open(gast_filename, 'w')
            if os.path.exists(clustergast_filename_single):
                in_gast_fh  = open(clustergast_filename_single, 'r')
                          
                for line in in_gast_fh:
                    
                    s = line.strip().split("\t")
                    if len(s) == 4:
                        read_id     = s[0]
                        refhvr_id   = s[1].split('|')[0]
                        distance    = s[2]
                        alignment   = s[3]
                        frequency   = 0
                    elif len(s) == 5:
                        read_id     = s[0]
                        refhvr_id   = s[1].split('|')[0]
                        distance    = s[2]
                        alignment   = s[3]
                        frequency   = s[4]
                    else:
                        logger.debug("gast_cleanup: wrong field count")
                    #print read_id, refhvr_id
                    # if this was in the gast table zero it out because it had a valid hit
                    # so we don't insert them as non-hits later
                    if read_id in nonhits:
                        del nonhits[read_id]
                        #print 'deleling', read_id
                    #print 'nonhits', nonhits
                    if read_id not in copies:
                        logger.info(read_id+' not in names file: Skipping')
                        continue
                        
                    # give the same ref and dist for each duplicate
                    for id in copies[read_id]:
                        
                        if id != read_id:
                            #print id, read_id, distance, refhvr_id  
                            gast_fh.write( id + "\t" + refhvr_id + "\t" + distance + "\t" + alignment +"\t"+frequency+"\n" )
                            
                                                   
                in_gast_fh.close()
                 
                #######################################
                # 
                #  Insert a record for any valid sequence that had no blast hit and therefore no gast result
                #       into gast_filename
                # 
                #######################################   
                for read in sorted(nonhits.iterkeys()):                
                    for d in copies[read]: 
                        gast_fh.write( d+"\t0\t1\t0\t0\n")
                        
                        
                gast_fh.close()
                
                # concatenate the two gast files
                clustergast_fh = open(clustergast_filename_single, 'a')            
                shutil.copyfileobj(open(gast_filename, 'rb'), clustergast_fh)
                clustergast_fh.close()
                #then open again and get data for gast concat
                concat = {}
                #print clustergast_filename_single
                for line in open(clustergast_filename_single, 'r'):
                    data = line.strip().split("\t")
                    id = data[0]
                    try:
                        refhvr_id = data[1].split('|')[0]
                    except:
                        refhvr_id = data[1]
                    distance = data[2]
                    #print 'data', data
                    if id in concat:
                        concat[id]['refhvrs'].append(refhvr_id)                        
                    else:
                        concat[id] = {}
                        concat[id]['refhvrs'] = [refhvr_id]
                    concat[id]['distance'] = distance     
                    
                
                
                #######################################
                #
                # Insert records into gast_concat_filename
                #
                #######################################             
                # first we need to open the gast_filename
                gastconcat_fh     = open(gastconcat_filename, 'w')
                for id, value in concat.iteritems():
                    #print 'trying gastconcat', id, value
                    gastconcat_fh.write( id + "\t" + concat[id]['distance'] + "\t" + ' '.join(concat[id]['refhvrs']) + "\n" )
                gastconcat_fh.close()
           
            else:
                logger.warning("No clustergast file found:"+clustergast_filename_single+"\nContinuing on ...")
                self.runobj.run_status_file_h.write(json.dumps({'status':'WARNING', 'message':"No clustergast file found: "+clustergast_filename_single+" Continuing"})+"\n")
  
            
        print "Finished gast_cleanup"   
        logger.info("Finished gast_cleanup")
        return {'status':"SUCCESS", 'message':"gast_cleanup finished"}

    def gast2tax(self):
        """
        Creates taxtax files
        """
        self.runobj.run_status_file_h.write(json.dumps({'status':'STARTING_GAST2TAX'})+"\n")
        counter = 0
        qsub_prefix = 'gast2tax_sub_'
        #print tax_file
        tax_files = []
        for key in self.iterator:
            counter += 1
            if self.runobj.vamps_user_upload:
                output_dir = os.path.join(self.global_gast_dir, key)
                gast_dir = os.path.join(self.global_gast_dir, key)
                fasta_file = os.path.join(output_dir, 'fasta.fa')
                unique_file = os.path.join(output_dir, 'unique.fa')
                names_file = os.path.join(output_dir, 'names')
                tagtax_file = os.path.join(output_dir, 'tagtax_terse')
                if os.path.exists(unique_file) and os.path.getsize(unique_file)>0:
                    tax_files.append(tagtax_file)
                #datasets_file = os.path.join(self.global_gast_dir, 'datasets')
            else:
                if self.runobj.platform == 'illumina':
                    output_dir = os.path.join(self.global_gast_dir, key)
                    gast_dir = os.path.join(self.global_gast_dir, key)
                    file_prefix = self.runobj.samples[key].file_prefix
                    unique_file = os.path.join(self.input_dir, file_prefix+"-PERFECT_reads.fa.unique")
                    names_file = os.path.join(self.input_dir, file_prefix+"-PERFECT_reads.fa.unique.names")
                elif self.runobj.platform == '454':
                    pass
                else:
                    sys.exit("gast2tax: no platform")
                
            if key in self.runobj.samples:
                dna_region = self.runobj.samples[key].dna_region
            else:            
                dna_region = self.runobj.dna_region
            if not dna_region:
                logger.error("gast2tax: We have no DNA Region: Setting dna_region to 'unknown'")
                self.runobj.run_status_file_h.write(json.dumps({'status':'WARNING', 'message':"gast2tax: We have no DNA Region: Setting dna_region to 'unknown'"})+"\n")
                dna_region = 'unknown'
            max_distance = C.max_distance['default']
            if dna_region in C.max_distance:
                max_distance = C.max_distance[dna_region] 
            if self.use_cluster:
                
                clusterize = C.clusterize_cmd
                # create script - each file gets script
                script_filename = os.path.join(gast_dir, qsub_prefix + str(counter))
                fh = open(script_filename, 'w')
                qstat_name = "gast2tax" + key + '_' + self.runobj.run + "_" + str(counter)
                log_file = os.path.join(gast_dir, 'gast2tax.log_' + str(counter))
                fh.write("#!/bin/sh\n\n")
                #fh.write("#$ -j y\n" )
                #fh.write("#$ -o " + log_file + "\n")
                #fh.write("#$ -N " + qstat_name + "\n\n")

                # setup environment
                fh.write("source /xraid/bioware/Modules/etc/profile.modules\n")
                fh.write("module load bioware\n")
                fh.write("export PYTHONPATH=/xraid2-2/vampsweb/"+self.runobj.site+"/:$PYTHONPATH\n\n")
                fh.write("/xraid2-2/vampsweb/"+self.runobj.site+"/pipeline/gast2tax.py -dna "+dna_region+" -max "+str(max_distance)+" -o "+gast_dir+" -n "+names_file+" -site "+self.runobj.site+"\n" )
                fh.close()
                            
                # make script executable and run it
                os.chmod(script_filename, stat.S_IRWXU)
                qsub_cmd = clusterize + " -log /xraid2-2/vampsweb/"+self.runobj.site+"/clusterize.log " + script_filename
                print qsub_cmd
                # on vamps and vampsdev qsub cannot be run - unless you call it from the
                # cluster aware directories /xraid2-2/vampsweb/vamps and /xraid2-2/vampsweb/vampsdev
                #qsub_cmd = C.qsub_cmd + " " + script_filename
                #qsub_cmd = C.clusterize_cmd + " " + script_filename
                logger.debug("qsub command: "+qsub_cmd)
                
                subprocess.call(qsub_cmd, shell=True)
                #proc = subprocess.Popen(qsub_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                
            else:
                if os.path.exists(names_file) and os.path.getsize(names_file) > 0: 
                    (refdb, taxdb) = self.get_reference_databases(dna_region)
                    
                    ref_taxa = self.load_reftaxa(taxdb)
            
                    self.assign_taxonomy(key, gast_dir, dna_region, names_file, ref_taxa);
                    
                    
        if self.use_cluster:
            # wait here for tagtax files to finish
            temp_file_list = tax_files
            tagtax_terse_filename     = os.path.join(gast_dir, "tagtax_terse")
            tagtax_long_filename     = os.path.join(gast_dir, "tagtax_long")
            c = False
            maxwaittime = C.maxwaittime  # seconds
            sleeptime   = C.sleeptime    # seconds
            counter3 = 0
            while c == False:
                print 'temp_file_list:', temp_file_list
                counter3 += 1
                if counter3 >= maxwaittime / sleeptime:
                    raise Exception("Max wait time exceeded in gast.py")
                for index, file in enumerate(temp_file_list):
                    #print temp_file_list
                    if os.path.exists(file):
                        # remove from tmp list
                        logger.debug("Found file now removing from list: "+file)
                        temp_file_list = temp_file_list[:index] + temp_file_list[index+1:]
                
                if temp_file_list:
                    logger.info("waiting for tagtax files to fill...")
                    logger.debug(' '.join(temp_file_list))
                    logger.info("\ttime: "+str(counter3 * sleeptime)+" | files left: "+str(len(temp_file_list)))
                    time.sleep(sleeptime)
                else:
                    c = True
            
            
            
#                 counter3 += 1
#                 if counter3 >= maxwaittime / sleeptime:
#                     raise Exception("Max wait time exceeded in gast.py: gast2tax")
#                 
#                 #print temp_file_list
#                 if os.path.exists(tagtax_long_filename) and os.path.getsize(tagtax_long_filename) > 100:
#                     # remove from tmp list
#                     logger.debug("Found file: "+tagtax_long_filename+" - Continuing")
#                     c = True
#                 else:
#                     logger.info("waiting for tagtax files to fill...")
#                     logger.info("\ttime: "+str(counter3 * sleeptime))
#                     time.sleep(sleeptime)
                    
                    
        print "Finished gast2tax" 
        return {'status':"SUCCESS", 'message':"gast2tax finished"} 
        
    
    def get_reference_databases(self, dna_region):
        
        #if dna region == v6v4(a) change it to v4v6
        # other reverse regions? 
        if dna_region == 'v6v4':
            dna_region = 'v4v6'
        if dna_region == 'v6v4a':
            dna_region = 'v4v6a'
        print  'dna_region ', dna_region 
        # try this first:
        if os.path.exists(os.path.join(self.refdb_dir, 'ref'+dna_region+'.udb')):
            refdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.udb')
            taxdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.tax')
        elif os.path.exists(os.path.join(self.refdb_dir, 'ref'+dna_region+'.fa')):
            refdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.fa')
            taxdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.tax')
        elif C.use_full_length or dna_region == 'unknown' or dna_region not in C.refdbs:
            if os.path.exists(os.path.join(self.refdb_dir, 'refssu.udb')):
                refdb = os.path.join(self.refdb_dir, 'refssu.udb')
                taxdb = os.path.join(self.refdb_dir, 'refssu.tax')
            elif os.path.exists(os.path.join(self.refdb_dir, 'refssu.fa')):
                refdb = os.path.join(self.refdb_dir, 'refssu.fa')
                taxdb = os.path.join(self.refdb_dir, 'refssu.tax')
        else:
            
            # try udb first
            if dna_region in C.refdbs:
                if os.path.exists(os.path.join(self.refdb_dir, C.refdbs[dna_region]+".udb")):
                    refdb = os.path.join(self.refdb_dir, C.refdbs[dna_region]+".udb")
                    taxdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.tax')
                elif os.path.exists(os.path.join(self.refdb_dir, C.refdbs[dna_region])):
                    refdb = os.path.join(self.refdb_dir, C.refdbs[dna_region])
                    taxdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.tax')
                else:
                    print 'could not find refdb '+os.path.join(self.refdb_dir, C.refdbs[dna_region])+".udb - Using full length"
                    refdb = os.path.join(self.refdb_dir, 'refssu.fa')
                    taxdb = os.path.join(self.refdb_dir, 'refssu.tax')
            elif os.path.exists(os.path.join(self.refdb_dir, 'ref'+dna_region+'.fa')):
                refdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.fa')
                taxdb = os.path.join(self.refdb_dir, 'ref'+dna_region+'.tax')
                
            elif os.path.exists(os.path.join(self.refdb_dir, 'refssu.udb')):
                refdb = os.path.join(self.refdb_dir, 'refssu.udb')
                taxdb = os.path.join(self.refdb_dir, 'refssu.tax')
                
            elif os.path.exists(os.path.join(self.refdb_dir, 'refssu.fa')):
            
                refdb = os.path.join(self.refdb_dir, 'refssu.fa')
                taxdb = os.path.join(self.refdb_dir, 'refssu.tax')
                
            else:
                refdb = os.path.join(self.refdb_dir, 'refssu.fa')
                taxdb = os.path.join(self.refdb_dir, 'refssu.tax')
                logger.error("Could not find reference database in "+self.refdb_dir+" - Using full length")
                
            
                
        logger.info('tax_file: '+taxdb)
        logger.info('ref_file: '+refdb) 
        return (refdb, taxdb)   
          
    def get_fastasampler_cmd(self, unique_file, fastasamp_filename, start, end):
        fastasampler = C.fastasampler_cmd
        if self.utils.is_local():
            fastasampler = "/Users/ashipunova/bin/illumina-utils/fastasampler"        
        fastasampler_cmd = fastasampler
        fastasampler_cmd += ' -n '+ str(start)+','+ str(end)
        fastasampler_cmd += ' ' + unique_file
        fastasampler_cmd += ' ' + fastasamp_filename        
        return fastasampler_cmd
        
    def get_usearch_cmd(self, fastasamp_filename, refdb, usearch_filename  ):    
        
#         if C.use_full_length:            
#             usearch6 = C.usearch64
#         else:
#             usearch6 = C.usearch6_cmd
#             
        usearch_cmd = C.usearch6_cmd
        if self.utils.is_local():
            usearch_cmd = "/Users/ashipunova/bin/illumina-utils/usearch"
        
        usearch_cmd += ' -usearch_global ' + fastasamp_filename
        usearch_cmd += ' -gapopen 6I/1E'
        usearch_cmd += ' -uc_allhits'
        usearch_cmd += ' -db ' + refdb  
        usearch_cmd += ' -strand plus'   
        usearch_cmd += ' -uc ' + usearch_filename 
        usearch_cmd += ' -maxaccepts ' + str(C.max_accepts)
        usearch_cmd += ' -maxrejects ' + str(C.max_rejects)
        usearch_cmd += ' -id ' + str(C.pctid_threshold)
        
        
#         usearch_cmd = C.usearch5_cmd
#         usearch_cmd += ' --global '
#         usearch_cmd += ' --query ' + fastasamp_filename
#         usearch_cmd += ' --iddef 3'
#         usearch_cmd += ' --gapopen 6I/1E'
#         usearch_cmd += ' --db ' + refdb               
#         usearch_cmd += ' --uc ' + usearch_filename 
#         usearch_cmd += ' --maxaccepts ' + str(C.max_accepts)
#         usearch_cmd += ' --maxrejects ' + str(C.max_rejects)
#         usearch_cmd += ' --id ' + str(C.pctid_threshold)
        
        
        
        return usearch_cmd
        
    def get_grep_cmd(self, usearch_filename, clustergast_filename ):
        
        use_full_length = ''
        if C.use_full_length:
            use_full_length = "-use_full_length"
        if hasattr(self.runobj, 'use_full_length') and self.runobj.use_full_length: 
            use_full_length = "-use_full_length"
            
        grep_cmd = "grep"
        grep_cmd += " -P \"^H\\t\" " + usearch_filename + " |"
        #grep_cmd += " sed -e 's/|.*\$//' |"
        # changed grep command to remove frequency from read_id
        #grep_cmd += " sed -e 's/|frequency:[0-9]*//' |"
        #grep_cmd += " awk -F\"\\t\" '{print $9 \"\\t\" $4 \"\\t\" $10 \"\\t\" $8}' |"
        
        grep_cmd += " sed -e 's/|.*\$//' |"
        grep_cmd += " awk -F\"\\t\" '{print $9 \"\\t\" $4 \"\\t\" $10 \"\\t\" $8}' |"
        
        grep_cmd += " sort -k1,1b -k2,2gr |"
        # append to clustergast file:
        # split_defline adds frequency to gastv6 file (last field)
        tophit_cmd = "/bioware/linux/seqinfo/bin/python_pipeline/py_mbl_sequencing_pipeline/pipeline/clustergast_tophit"
        grep_cmd += " "+tophit_cmd + " -split_defline_frequency "+use_full_length+" >> " + clustergast_filename
    
        return grep_cmd  
          
    def load_reftaxa(self, tax_file):
    
        
        taxa ={}
        #open(TAX, "<$tax_file") || die ("Unable to open reference taxonomy file: $tax_file.  Exiting\n");
        #while (my $line = <TAX>) 
        n=1
        for line in  open(tax_file, 'r'):
        
            # 0=ref_id, 1 = taxa, 2 = count
            data=line.strip().split("\t")
    
            copies = []
    
            # foreach instance of that taxa
            for i in range(0, int(data[2])):
            
                # add that taxonomy to an array
                copies.append(data[1])
            
            # add that array to the array of all taxa for that ref, stored in the taxa hash
            if data[0] in taxa:
                taxa[data[0]].append(copies)
            elif copies:         
                taxa[data[0]] =[copies]
            else:
                taxa[data[0]] =[]
            n += 1
        return taxa
    
    def assign_taxonomy(self, key, gast_dir, dna_region, names_file, ref_taxa):
        from pipeline.taxonomy import Taxonomy, consensus
        #results = uc_results
        results = {}
 
 
        self.runobj.run_status_file_h.write(json.dumps({'status':"STARTING_ASSIGN_TAXONOMY: "+key})+"\n")
        #test_read='FI1U8LC02GEF7N'
        # open gast_file to get results
        tagtax_terse_filename     = os.path.join(gast_dir, "tagtax_terse")
        tagtax_long_filename     = os.path.join(gast_dir, "tagtax_long")
        tagtax_terse_fh = open(tagtax_terse_filename, 'w')
        tagtax_long_fh = open(tagtax_long_filename, 'w')
        tagtax_long_fh.write("\t".join(["read_id", "taxonomy", "distance", "rank", "refssu_count", "vote", "minrank", "taxa_counts", "max_pcts", "na_pcts", "refhvr_ids"])+"\n")
        gast_file          = os.path.join(gast_dir, "gast"+dna_region)
        if not os.path.exists(gast_file):
            logger.info("gast:assign_taxonomy: Could not find gast file: "+gast_file+". Returning")
            return results
        for line in  open(gast_file, 'r'): 
            # must split on tab because last field may be empty and must be maintained as blank
            data=line.strip().split("\t")
            if len(data) == 3:
                data.append("")
            # 0=id, 1=ref, 2=dist, 3=align 4=frequency
            #if data[0]==test_read:
            #    print 'found test in gastv6 ', data[1].split('|')[0], data[2], data[3]
                
            read = data[0]
            if read in results:
                results[read].append( [data[1].split('|')[0], data[2], data[3], data[4]] )
            else:            
                results[read]=[ [data[1].split('|')[0], data[2], data[3], data[4]] ]
            
        
        for line in open(names_file, 'r'):
            data=line.strip().split("\t")
            dupes = data[1].split(",")
            read  = data[0]
            taxObjects  = []
            distance    = 0
            frequency   = 0
            refs_for    = {}
            
            #print 'read', read
            if read not in results:
                results[read]=["Unknown", '1', "NA", '0', '0', "NA", "0;0;0;0;0;0;0;0", "0;0;0;0;0;0;0;0", "100;100;100;100;100;100;100;100"]
                refs_for[read] = [ "NA" ]
            else:
                #print 'read in res', read, results[read]
                #if read == test_read:
                #    print 'found ', test_read, results[test_read]
                for i in range( 0, len(results[read])):
                    #for resultread in results[read]:
                    #print 'resread', results[read]
                    ref = results[read][i][0]
                    if ref in ref_taxa:
                        for tax in ref_taxa[ref]:
                            for t in tax:
                                taxObjects.append(Taxonomy(t))
                    else:
                        pass
                        
                    if read in refs_for:
                        #if read ==test_read:
                        #    print '2', read, refs_for[test_read]
                        if results[read][i][0] not in refs_for[read]:
                            refs_for[read].append(results[read][i][0])  
                    else:
                        #if read == test_read:
                        #    print '1', read, results[read][i][0]
                        refs_for[read] = [results[read][i][0]]                   
                     
                    # should all be the same distance
                    distance = results[read][i][1]
                    frequency = results[read][i][3]
                #Lookup the consensus taxonomy for the array
                taxReturn = consensus(taxObjects, C.majority)
                
                # 0=taxObj, 1=winning vote, 2=minrank, 3=rankCounts, 4=maxPcts, 5=naPcts;
                taxon = taxReturn[0].taxstring()
                #if taxon[-3:] = ';NA':
                #    taxon = taxon[:-3]
                #tax_counter[taxon]
                rank = taxReturn[0].depth()
                #print read, taxon, rank, taxReturn[0], taxReturn[1]
                if not taxon: taxon = "Unknown"
            
                # (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts)
                results[read] = [ taxon, str(distance), rank, str(len(taxObjects)), str(taxReturn[1]), taxReturn[2], taxReturn[3], taxReturn[4], taxReturn[5] ] 
                #print "\t".join([read, taxon, str(distance), rank, str(len(taxObjects)), str(taxReturn[1]), taxReturn[2], taxReturn[3], taxReturn[4], taxReturn[5]]) + "\n"
#read_id taxonomy        distance        rank    refssu_count    vote    minrank taxa_counts     max_pcts        na_pcts refhvr_ids
#D4ZHLFP1:25:B022DACXX:3:1101:12919:40734 1:N:0:TGACCA|frequency:162     Bacteria;Proteobacteria;Gammaproteobacteria     0.117   class   2       100     genus   1;1;1;2;2;2;0;0 100;100;100;50;50;50;0;0        0;0;0;0;0;0;100;100     v6_CI671
#D4ZHLFP1:25:B022DACXX:3:1101:10432:76870 1:N:0:TGACCA|frequency:105     Bacteria;Proteobacteria;Gammaproteobacteria     0.017   class   1       100     class   1;1;1;0;0;0;0;0 100;100;100;0;0;0;0;0   0;0;0;100;100;100;100;100       v6_BW306
                
            # Replace hash with final taxonomy results, for each copy of the sequence
            for d in dupes:
               # print OUT join("\t", $d, @{$results{$read}}, join(", ", sort @{$refs_for{$read}})) . "\n";
                d = d.strip()
                tagtax_long_fh.write( d+"\t"+"\t".join(results[read])+"\t"+', '.join(sorted(refs_for[read]))  + "\n")
                tagtax_terse_fh.write(d+"\t"+results[read][0]+"\t"+results[read][2]+"\t"+results[read][3]+"\t"+', '.join(sorted(refs_for[read]))+"\t"+results[read][1]+"\t"+str(frequency)+"\n")
               
        tagtax_terse_fh.close()
        tagtax_long_fh.close()
        return results

#    def create_uniques_from_fasta(self, fasta_file, key):
#         
#         mothur_cmd = C.mothur_cmd+" \"#unique.seqs(fasta="+fasta_file+", outputdir="+os.path.join(self.basedir, key)+"/);\""; 
#         
#         #mothur_cmd = site_base+"/clusterize_vamps -site vampsdev -rd "+user+"_"+runcode+"_gast -rc "+runcode+" -u "+user+" /bioware/mothur/mothur \"#unique.seqs(fasta="+fasta_file+");\"";    
#         subprocess.call(mothur_cmd, shell=True)
    def check_for_unique_files(self, keys):
        logger.info("Checking for uniques file")
        for key in keys:
            if self.runobj.vamps_user_upload:
                # one fasta file or (one project and dataset from db)
                # if self.runobj.fasta_file is not None then we should have multiple datasets
                # which have already been uniqued
                if self.runobj.fasta_file and os.path.exists(self.runobj.fasta_file):
                    #output_dir = os.path.join(self.basedir, keys[0])
                    unique_file = os.path.join(self.global_gast_dir, key, 'unique.fa')
                    names_file = os.path.join(self.global_gast_dir, key, 'names')
                    
                    # the -x means do not store frequency data in defline of fasta file
                    fastaunique_cmd = C.fastaunique_cmd +" -x -i "+self.runobj.fasta_file+" -o "+unique_file+" -n "+names_file 
                    print fastaunique_cmd
                    
                    subprocess.call(fastaunique_cmd, shell=True)
                    
                    #shutil.move('a.txt', 'b.kml')
                    #os.rename(filename, filename[7:])
                    #os.rename(filename, filename[7:])
                else:
                    if self.runobj.project and self.runobj.dataset:
                        pass
                    else:
                        pass
                #get from database
            else:
                if self.runobj.platform == 'illumina':
#                    reads_dir = dirs.check_dir(dirs.reads_overlap_dir)
#                    os.path.join(self.analysis_dir, 'reads_overlap')
                    
                    file_prefix = key
                    unique_file = os.path.join(self.reads_dir, file_prefix+"-PERFECT_reads.fa.unique")
                    if os.path.exists(unique_file):
                        if os.path.getsize(unique_file) > 0:
                            logger.debug( "GAST: Found uniques file: "+unique_file)
                        else:
                            logger.warning( "GAST: Found uniques file BUT zero size "+unique_file)
                            continue
                    else:
                        logger.error( "GAST: NO uniques file found "+unique_file)
                        
                        
                if self.runobj.platform == '454':
                    pass
                else:
                    pass
                
            
#         for key in keys:
#             fasta_file = ""
#             output_dir = os.path.join(self.basedir, key)
#             unique_file = os.path.join(output_dir, key+'.unique.fa')
#             if not os.path.exists(unique_file):
#                 mothur_cmd = C.mothur_cmd+" \"#unique.seqs(fasta="+fasta_file+", outputdir="+os.path.join(self.basedir, key)+"/);\""; 
#         
#                 #mothur_cmd = site_base+"/clusterize_vamps -site vampsdev -rd "+user+"_"+runcode+"_gast -rc "+runcode+" -u "+user+" /bioware/mothur/mothur \"#unique.seqs(fasta="+fasta_file+");\"";    
#                 subprocess.call(mothur_cmd, shell=True)
        return {"status":"SUCCESS", "message":"checking for uniques"}       
                
    def get_fasta_from_database(self):
        pass