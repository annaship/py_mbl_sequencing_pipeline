import subprocess
import sys, os,stat
import time
from pipeline.pipelinelogging import logger
#import logging
import constants as C

class Gast:
    """Doc string here.."""
    Name = "GAST"
    def __init__(self, run = None):

        self.run 	 = run
        self.outdir  = run.output_dir
        self.rundate = self.run.run_date
        self.use_cluster = 1
        
        
        os.environ['SGE_ROOT']='/usr/local/sge'
        os.environ['SGE_CELL']='grendel'
        path = os.environ['PATH']
        os.environ['PATH'] = path + ':/usr/local/sge/bin/lx24-amd64'
        #First step is to check for (or create via mothur)
        # a uniques fasta file and a names file 
        # one for each dataset.
        # If we are here from a vamps gast process
        # then there should be just one dataset to gast
        # but if MBL pipe then many datasets are prbably involved.
        
        # 1) clustergast
        # 2) gast cleanup
        # 3) gast2tax
        
        
    def clustergast(self, lane_keys):
    	"""
    	clustergast - runs the GAST pipeline on the cluster.
               GAST uses UClust to identify the best matches of a read sequence
               to references sequences in a reference database.
        """
    	# Step1: create empty gast table in database: gast_<rundate>
    	# Step2: Count the number of sequences so the job can be split for nodes
    	# $facount = `grep -c \">\" $fasta_uniqs_filename`;
    	# $calcs = `/bioware/seqinfo/bin/calcnodes -t $facount -n $nodes -f 1`;
    
    	#	/bioware/seqinfo/bin/fastasampler -n $start,$end ${gastDir}/${fasta_uniqs_filename} $tmp_fasta_filename
    	#	$usearch_binary --global --query $tmp_fasta_filename --iddef 3 --gapopen 6I/1E --db $refhvr_fa --uc $tmp_usearch_filename --maxaccepts $max_accepts --maxrejects $max_rejects --id $pctid_threshold
    	#	# sort the results for valid hits, saving only the ids and pct identity
        #	grep -P \"^H\\t\" $tmp_usearch_filename | sed -e 's/|.*\$//' | awk '{print \$9 \"\\t\" \$4 \"\\t\" \$10 \"\\t\" \$8}' | sort -k1,1b -k2,2gr | clustergast_tophit > $gast_filename
    	#	Submit the script
    	#	/usr/local/sge/bin/lx24-amd64/qsub $qsub_priority $script_filename
    	
    	usearch = '/bioware/uclust/usearch'
    	fastasampler = '/bioware/seqinfo/bin/fastasampler'
    	calcnodes = '/bioware/seqinfo/bin/calcnodes'
    	sqlImportCommand = '/usr/bin/mysqlimport'
    	#qsub = '/usr/local/sge/bin/lx24-amd64/qsub'
    	qsub = '/bioware/seqinfo/bin/clusterize'
    	refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
    	

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
        for lane_key in lane_keys:
            nodes = 100
            dna_region = self.run.samples[lane_key].dna_region
            #direction  = self.run.samples[lane_key].direction
            
            refdb = refdb_dir + 'ref' + dna_region
            if not os.path.exists(refdb):
            	if dna_region == 'v6v4':
            		refdb = refdb_dir + 'refv4v6'
            
            unique_file = self.outdir +'/'+lane_key+'.unique.fa'
            
            gast_dir = self.outdir +'/'+lane_key+'_gast/'
            if not os.path.exists(gast_dir):	
            	os.mkdir(gast_dir)
            	
            clustergast_filename   = gast_dir + lane_key+".clustergast"
            
            grep_cmd = ['grep','-c','>',unique_file]
            logger.debug( grep_cmd )
            facount = subprocess.check_output(grep_cmd)
            logger.debug( lane_key+' count '+facount)
            calcnode_cmd = [calcnodes,'-t',str(facount),'-n',str(nodes),'-f','1']
            
            calcout = subprocess.check_output(calcnode_cmd)
            logger.debug("calcout:\n"+calcout)
            #calcout:
            # node=1 start=1 end=1 rows=1
            # node=2 start=2 end=2 rows=1
            # node=3 start=3 end=3 rows=1           
            lines = calcout.split("\n")
            i = 1
            for line in lines:
            	if i == nodes + 1:
            		break
                script_filename = gast_dir + qsub_prefix + str(i)
                gast_filename   = gast_dir + gast_prefix + str(i)
                fastasamp_filename = gast_dir + 'samp_' + str(i)
                usearch_filename= gast_dir + "uc_" + str(i)
                log_file = gast_dir + 'clustergast.log_' + str(i)
                
                if self.use_cluster:
                    fh = open(script_filename,'w')
                    qstat_name = "gast" + lane_key + '_' + self.rundate + "_" + str(i)
                    fh.write("#!/bin/csh\n")
                    fh.write("#$ -j y\n" )
                    fh.write("#$ -o " + log_file + "\n")
                    fh.write("#$ -N " + qstat_name + "\n\n")
                    fh.write("source /xraid/bioware/Modules/etc/profile.modules\n");
                    fh.write("module load bioware\n\n");

                data = line.split()
                
                if len(data) < 2:
                	continue
                start = data[1].split('=')[1]
                end  = data[2].split('=')[1]

                fastasampler_cmd = fastasampler
                fastasampler_cmd += ' -n '+ str(start)+','+ str(end)
                fastasampler_cmd += ' ' + unique_file
                fastasampler_cmd += ' ' + fastasamp_filename
                logger.debug("fastasampler command: "+fastasampler_cmd)
                if self.use_cluster:
                	fh.write(fastasampler_cmd + "\n")
                else:
                	#fastasampler_cmd="fastasampler -n 1,73 20100917/1_GATGA.unique.fa 20100917/1_GATGA_gast/samp_1"
                	subprocess.call(fastasampler_cmd,shell=True)
                	
                usearch_cmd = usearch
                usearch_cmd += ' --global'
                usearch_cmd += ' --query ' + fastasamp_filename
                usearch_cmd += ' --iddef 3'
                usearch_cmd += ' --gapopen 6I/1E'
                usearch_cmd += ' --db ' + refdb                
                usearch_cmd += ' --uc ' + usearch_filename 
                usearch_cmd += ' --maxaccepts ' + str(C.max_accepts)
                usearch_cmd += ' --maxrejects ' + str(C.max_rejects)
                usearch_cmd += ' --id ' + str(C.pctid_threshold)
                logger.debug("usearch command: "+usearch_cmd)
                if self.use_cluster:
                	fh.write(usearch_cmd + "\n")
                	
                else:
                	subprocess.call(usearch_cmd,shell=True)
                
                
                grep_cmd = "grep"
                grep_cmd += " -P \"^H\\t\" " + usearch_filename + " |"
                grep_cmd += " sed -e 's/|.*\$//' |"
                grep_cmd += " awk '{print $9 \"\\t\" $4 \"\\t\" $10 \"\\t\" $8}' |"
                grep_cmd += " sort -k1,1b -k2,2gr |"
                # append to file:
                grep_cmd += " clustergast_tophit >> " + clustergast_filename
                logger.debug("grep command: "+grep_cmd)
                if self.use_cluster:                
                    fh.write(grep_cmd + "\n")
                    fh.close()
                    
                    # make script executable and run it
                    os.chmod(script_filename, stat.S_IRWXU)
                    qsub_cmd = qsub + " " + script_filename
                    logger.debug("qsub command: "+qsub_cmd)
                    subprocess.call(qsub_cmd, shell=True)
                else:
                	subprocess.call(grep_cmd,shell=True)
                	
                # We'll aim to not use the db until done
                # so write theses records to a file instead
                #mysqlimport_cmd = sqlImportCommand 
                #mysqlimport_cmd += " -u " + db_user 
                #mysqlimport_cmd += " -p"+db_password 
                #mysqlimport_cmd += " -C -v -L"
                #mysqlimport_cmd += " -h " + db_hostname+" "+dbName 
                #mysqlimport_cmd += " " + gast_filename
                #fh.write(mysqlimport_cmd + "\n")
                              
                i += 1
            
                  
            logger.info("Finished clustergast")

     
    def gast_cleanup(self, lane_keys):
    	"""
    	gast_cleanup - follows clustergast, explodes the data and copies to gast_concat and gast
    	"""
    	for lane_key in lane_keys:
            unique_file = self.outdir +'/'+lane_key+'.unique.fa'
            names_file = self.outdir +'/'+lane_key+'.names'
            gast_dir = self.outdir +'/'+lane_key+'_gast/'
            clustergast_filename   = gast_dir + lane_key+".clustergast"
            gast_filename = gast_dir+lane_key+".gast";
            gast_file = gast_dir+'gast'+self.dna_region+'_'+lane_key
            gast_concat_file = gast_dir+'gast_concat'+self.dna_region+'_'+lane_key

            copies = {}
            nonhits = {}
    		# open and read names file
            names_fh = open(names_file,'r')
            for line in names_fh:
                (index_read,dupes) = line.strip().split("\t")
                # initially nonhits is all unique readids
                if index_read in nonhits:
                	nonhits[index_read] += 1
                else:
                	nonhits[index_read] = 1
                copies[index_read] = dupes.split(',')
                
                print index_read
            names_fh.close()
            
            # open and read the gast file from clustergast
            in_gast_fh         = open(clustergast_filename,'r')
            
            
            #out_filename_nonhits = gast_dir+lane_key+".nonhits";
            out_gast_fh        = open(gast_filename,'w')
            #out_nohits_fh = open(out_filename_nonhits,'w')
            ########################################
            #
    		# Step through the gast records and insert the dupes
    		#
    		########################################
            logger.debug("length of nonhits: Before:"+len(nonhits))
            logger.info("Inserting BLAST hits into Gast File");
            for line in in_gast_fh:
            	(read_id, refhvr_id, distance, alignment) = line.strip().split()
            	
            	# very important!
            	# if this was in the gast table zero it out because it had a valid hit
        		# so we don't insert them as non-hits later
            	if read_id in nonhits:
#                    delete nonhits[read_id]
                    pass
            	    
                if not read_id in copies:
                    logger.info(read_id+' not in names file: Skipping')
                    continue
                # give the same ref and dist for each duplicate
                for id in copies[read_id]:
                	if id != read_id:
                		out_gast_fh.write(id+"\t"+refhvr_id+"\t"+distance+"\t"+alignment+"\n")
            
            logger.debug("length of nonhits: After:"+len(nonhits))

            #######################################
    		#
    		# Insert a record for any valid sequence that had no blast hit and therefore no gast result
    		#
    		#######################################
            logger.info("Inserting non-BLAST hits into Gast File");
        	# for the list of remaining reads add them in with their duplicates
            for read in nonhits.sort():
                for d in copies[read]:
                    out_gast_fh.write(d+"\t0\t1\t\n")
        	
        	out_gast_fh.close()
        	
        # so gast_filename is the gast file equivalent to the gastTable
        # now the gast_concat table needs to be calculated from this
        # if the gast table was in the database we'd need to do a group by
        # insert $ignore INTO gastConcatTable SELECT read_id, distance, 
    	# 		GROUP_CONCAT(DISTINCT refhvr_id ORDER BY refhvr_id ASC SEPARATOR ' ') 
    	# 		as refhvr_ids FROM gastTable GROUP BY read_id
		# how do we do this from a file????
		fh = open()
            
        logger.info("Finished gast_cleanup")
        
    def gast2tax(self, lane_keys):
        """
        Follows gast_cleanup
        """
        logger.info("Finished gast2tax")