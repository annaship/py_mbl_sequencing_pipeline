import subprocess
import sys, os,stat
import time
import shutil
from pipeline.pipelinelogging import logger
#import logging
import constants as C

class Gast:
    """Doc string here.."""
    Name = "GAST"
    def __init__(self, run = None):

        self.run     = run
        self.outdir  = run.output_dir
        try:
            self.basedir = run.basedir
        except:
            self.basedir = self.outdir
        self.rundate = self.run.run_date
        self.use_cluster = 1
        self.vamps_user_upload = run.vamps_user_upload
        
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
        if self.vamps_user_upload:
            self.refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
        else:
            self.refdb_dir = '/xraid2-2/g454/blastdbs/'
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

        #   /bioware/seqinfo/bin/fastasampler -n $start,$end ${gastDir}/${fasta_uniqs_filename} $tmp_fasta_filename
        #   $usearch_binary --global --query $tmp_fasta_filename --iddef 3 --gapopen 6I/1E --db $refhvr_fa --uc $tmp_usearch_filename --maxaccepts $max_accepts --maxrejects $max_rejects --id $pctid_threshold
        #   # sort the results for valid hits, saving only the ids and pct identity
        #   grep -P \"^H\\t\" $tmp_usearch_filename | sed -e 's/|.*\$//' | awk '{print \$9 \"\\t\" \$4 \"\\t\" \$10 \"\\t\" \$8}' | sort -k1,1b -k2,2gr | clustergast_tophit > $gast_filename
        #   Submit the script
        #   /usr/local/sge/bin/lx24-amd64/qsub $qsub_priority $script_filename
 
        usearch = '/bioware/uclust/usearch'
        fastasampler = '/bioware/seqinfo/bin/fastasampler'
        calcnodes = '/bioware/seqinfo/bin/calcnodes'
        sqlImportCommand = '/usr/bin/mysqlimport'
        #qsub = '/usr/local/sge/bin/lx24-amd64/qsub'
        qsub = '/bioware/seqinfo/bin/clusterize'



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
            if lane_key in self.run.samples:
                dna_region = self.run.samples[lane_key].dna_region
            else:            
                dna_region = self.run.dna_region
            if not dna_region:
                logger.error("clustergast: We have no DNA Region: Setting dna_region to 'unknown'")
                dna_region = 'unknown'

            # if no dna_region OR no refdb can be found then use
            # refssu
            refdb = self.refdb_dir + C.refdbs[dna_region]
            if not os.path.exists(refdb):
                refdb = self.refdb_dir + 'refssu_all'
                    
                    
            unique_file = self.basedir +'/'+lane_key+'.unique.fa'
            
            gast_dir = self.outdir +'/'+lane_key+'_gast/'
            if not os.path.exists(gast_dir):	
                os.mkdir(gast_dir)
            else:
                # empty then recreate directory
                shutil.rmtree(gast_dir)
                os.mkdir(gast_dir)


            grep_cmd = ['grep','-c','>',unique_file]
            logger.debug( ' '.join(grep_cmd) )
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
            file_list = []
            i = 0
            for line in lines:
                i += 1
                if i >= nodes:
                    continue
                script_filename = gast_dir + qsub_prefix + str(i)
                gast_filename   = gast_dir + gast_prefix + str(i)
                fastasamp_filename = gast_dir + 'samp_' + str(i)
                clustergast_filename   = gast_dir + lane_key+".gast_" + str(i)
                file_list.append(clustergast_filename)
                usearch_filename= gast_dir + "uc_" + str(i)
                log_file = gast_dir + 'clustergast.log_' + str(i)
                
                data = line.split()
                
                if len(data) < 2:
                    continue
                start = data[1].split('=')[1]
                end  = data[2].split('=')[1]
                
                if self.use_cluster:
                    fh = open(script_filename,'w')
                    qstat_name = "gast" + lane_key + '_' + self.rundate + "_" + str(i)
                    fh.write("#!/bin/csh\n")
                    fh.write("#$ -j y\n" )
                    fh.write("#$ -o " + log_file + "\n")
                    fh.write("#$ -N " + qstat_name + "\n\n")
                    #fh.write("source /xraid/bioware/Modules/etc/profile.modules\n");
                    #fh.write("module load bioware\n\n");



                fastasampler_cmd = fastasampler
                fastasampler_cmd += ' -n '+ str(start)+','+ str(end)
                fastasampler_cmd += ' ' + unique_file
                fastasampler_cmd += ' ' + fastasamp_filename
                logger.debug("fastasampler command: "+fastasampler_cmd)
                if self.use_cluster:
                    fh.write(fastasampler_cmd + "\n")
                else:
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
                    #subprocess.call(qsub_cmd, shell=True)
                    proc = subprocess.Popen(qsub_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    # proc.communicate will block - probably not what we want
                    #(stdout, stderr) = proc.communicate() #block the last onehere
                    #print stderr,stdout

                else:
                    subprocess.call(grep_cmd,shell=True)


            # wait here for all the clustergast scripts to finish
            temp_file_list = file_list
            c = False
            maxwaittime = C.maxwaittime  # seconds
            sleeptime   = C.sleeptime    # seconds
            counter = 0
            while c == False:
                counter += 1
                if counter >= maxwaittime / sleeptime:
                    raise Exception("Max wait time exceeded in gast.py")
                for index, file in enumerate(temp_file_list):
                    #print temp_file_list
                    if os.path.exists(file) and os.path.getsize(file) > 0:
                        # remove from tmp list
                        #logger.debug("Found file now removing from list: "+file)
                        temp_file_list = temp_file_list[:index] + temp_file_list[index+1:]
                
                if temp_file_list:
                    print "waiting for clustergast files to fill..."
                    print "    time:",counter * sleeptime,"| files left:",len(temp_file_list)
                    time.sleep(sleeptime)
                else:
                    c = True
                    
            # now concatenate all the clustergast_files into one file
            clustergast_filename_single   = gast_dir + lane_key+".gast"
            clustergast_fh = open(clustergast_filename_single,'w')
            # have to turn off cluster above to be able to 'find' these files for concatenation
            for n in range(1,i-1):
                #cmd = "cat "+ gast_dir + lane_key+".gast_" + str(n) + " >> " + gast_dir + lane_key+".gast"
                file = gast_dir + lane_key+".gast_" + str(n)
                if(os.path.exists(file)):                    
                    shutil.copyfileobj(open(file,'rb'), clustergast_fh)
                else:
                    print "Could not find file:",os.path.basename(file)," Skipping"
            print "Finished clustergast concatenation"
            clustergast_fh.flush()
            clustergast_fh.close()
            
            # remove tmp files
            for n in range(1,i-1):
                if os.path.exists(gast_dir+"uc_"+str(n)):
                    os.remove(gast_dir+"uc_"+str(n))
                    pass
                if os.path.exists(gast_dir+"samp_"+str(n)):    
                    os.remove(gast_dir+"samp_"+str(n))
                    pass


            logger.info("Finished clustergast")

     
    def gast_cleanup(self, lane_keys):
        """
        gast_cleanup - follows clustergast, explodes the data and copies to gast_concat and gast files
        """
        for lane_key in lane_keys:
            if lane_key in self.run.samples:
                dna_region = self.run.samples[lane_key].dna_region
            else:            
                dna_region = self.run.dna_region
            if not dna_region:
                logger.error("gast_cleanup: We have no DNA Region: Setting dna_region to 'unknown'")
                dna_region = 'unknown'
            # for vamps user upload
            # basedir is like avoorhis_3453211
            # and outdir is like avoorhis_3453211/2012-06-25
            # for MBL pipeline
            # basedir is like 1_AGTCG
            # and outdir is like 1_AGTCG/2012-06-25
            
            unique_file = self.basedir +'/'+lane_key+'.unique.fa'
            names_file = self.basedir +'/'+lane_key+'.names'
            #print 'names file',names_file
            gast_dir = self.outdir +'/'+lane_key+'_gast/'
            if not os.path.exists(gast_dir):
                logger.error("Could not find gast directory: "+gast_dir+" Exiting")
                sys.exit()
            clustergast_filename   = gast_dir + lane_key+".gast"
            logger.debug('gast filesize:'+str(os.path.getsize(clustergast_filename)))
            
            gast_filename          = gast_dir + lane_key+".gast_table";
            gastconcat_filename    = gast_dir + lane_key+".gast_concat_table";        

            copies = {}
            nonhits = {}
            # open and read names file
            names_fh = open(names_file,'r')
            for line in names_fh:
                s = line.strip().split("\t")
                if len(s) == 2:
                    index_read = s[0]
                    dupes = s[1]
                
                if index_read in nonhits:
                    nonhits[index_read] += 1
                else:
                    nonhits[index_read] = 1
                copies[index_read] = dupes.split(',')
                
            names_fh.close()            
            
            
            #######################################
            # 
            #  Insert records with valid gast hits into gast_file
            # 
            #######################################   
            # read the .gast file from clustergast            
            concat = {}
            gast_fh     = open(gast_filename,'w')
            if(os.path.exists(clustergast_filename)):
                in_gast_fh  = open(clustergast_filename,'r')
            else:
                print "No clustergast file found:",clustergast_filename,"\nExiting"
                sys.exit()
            for line in in_gast_fh:
                #print line.strip()
                s = line.strip().split()
                if len(s) == 4:
                    read_id     = s[0]
                    refhvr_id   = s[1]
                    distance    = s[2]
                    alignment   = s[3]
               
                # if this was in the gast table zero it out because it had a valid hit
                # so we don't insert them as non-hits later
                if read_id in nonhits:
                    del nonhits[read_id]
                if not read_id in copies:
                    logger.info(read_id+' not in names file: Skipping')
                    continue
                # give the same ref and dist for each duplicate
                for id in copies[read_id]:                    
                    if id != read_id:
                        if id not in concat:
                            concat[id] = {}
                            concat[id]['refhvrs'] = []
                        gast_fh.write( id + "\t" + refhvr_id + "\t" + distance + "\t" + alignment + "\n" )
                        concat[id]['refhvrs'].append(refhvr_id)
                        concat[id]['distance'] = distance                        
            in_gast_fh.close()
             
            #######################################
            # 
            #  Insert a record for any valid sequence that had no blast hit and therefore no gast result
            #       into gast_filename
            # 
            #######################################   
            for read in sorted(nonhits.iterkeys()):                
                for d in copies[read]: 
                    if d not in concat:
                        concat[d] = {}
                        concat[d]['refhvrs'] = []
                    gast_fh.write( d+"\t0\t1\t\n")
                    concat[d]['refhvrs'].append('0')
                    concat[d]['distance'] = '1'
                    
            gast_fh.close()
            
            #######################################
            #
            # Insert records into gast_concat_filename
            #
            #######################################             
            # first we need to open the gast_filename
            gastconcat_fh     = open(gastconcat_filename,'w')
            for id, value in concat.iteritems():
                #print 'trying gastconcat', id,value
                gastconcat_fh.write( id + "\t" + concat[id]['distance'] + "\t" + ' '.join(concat[id]['refhvrs']) + "\n" )
            gastconcat_fh.close()
            
            
        print "Finished gast_cleanup"   
        logger.info("Finished gast_cleanup")
        

        
    def gast2tax(self, lane_keys):
        """
        gast2tax : Follows gast_cleanup; assign taxonomy to read_ids in a gast table or file
        """
        from pipeline.taxonomy import Taxonomy,consensus
        
        for lane_key in lane_keys:
            if lane_key in self.run.samples:
                dna_region = self.run.samples[lane_key].dna_region
            else:            
                dna_region = self.run.dna_region
            if not dna_region:
                logger.error("gast2tax: We have no DNA Region: Setting dna_region to 'unknown'")
                dna_region = 'unknown'
                #sys.exit()
            
            # if no dna_region OR no refdb can be found then use
            # refssu
            refdb = self.refdb_dir + C.refdbs[dna_region]
            if not os.path.exists(refdb):
                refdb = self.refdb_dir + 'refssu_all'
            
            max_distance = C.max_distance['default']
            if dna_region in C.max_distance:
                max_distance = C.max_distance[dna_region]
            # for vamps user upload
            # basedir is like avoorhis_3453211
            # and outdir is like avoorhis_3453211/2012-06-25
            # for MBL pipeline
            # basedir is like 1_AGTCG
            # and outdir is like 1_AGTCG/2012-06-25
            
            unique_file = self.basedir +'/'+lane_key+'.unique.fa'
            names_file  = self.basedir +'/'+lane_key+'.names'
            gast_dir    = self.outdir  +'/'+lane_key+'_gast/'
            # 
            if not os.path.exists(gast_dir):
                logger.error("Could not find gast directory: "+gast_dir+" Exiting")
                sys.exit()
            tagtax_filename     = gast_dir + lane_key+".tagtax_table"
            gast_filename       = gast_dir + lane_key+".gast";
            if not os.path.exists(gast_filename):
                logger.error("Could not find gast file from gast_cleanup: "+os.path.basename(gast_filename)+" Exiting")
                sys.exit()
            gastconcat_filename = gast_dir + lane_key+".gast_concat_table";  
            if not os.path.exists(gastconcat_filename):
                logger.error("Could not find gast_concat file from gast_cleanup: "+os.path.basename(gastconcat_filename)+" Exiting")
                sys.exit()
            #######################################
            #
            # Load up the gast references
            #
            #######################################
            # from gast  file
            gast_fh = open(gast_filename,'r')
            distances = {}
            ref_hits = {}
            for line in gast_fh:
                items       = line.strip().split("\t")
                id          = items[0]
                refhvr_id   = items[1].split('|')[0]
                distance    = items[2]
                distances[id] = distance
                if id in ref_hits:                
                    ref_hits[id].append(refhvr_id)
                else:
                    ref_hits[id] = [refhvr_id]
            gast_fh.close()
            #for id in ref_hits:
            #    print id,ref_hits[id]
            #print "\n\nref_hits",ref_hits
            #######################################
            #
            # Load up the ref database
            #
            #######################################
            # now get stuff from ref table
            # if this is from vamps we will look on the vamps or vampsdev db
            # but normally we will look in the env454 db for the table
            # could we do this from a ref file????
            # SQL call to db (env454, vamps or vampsdev - depending....)
            refhvr_taxa = {}
            refdb_fh = open(refdb,'r')
            for line in refdb_fh:
                if line.startswith(">"):
                    items = line.strip().split("|")
                    refhvr_id = items[0][1:]
                    tax       = items[-1]
                    if refhvr_id in refhvr_taxa:
                        refhvr_taxa[refhvr_id].append(tax)
                    else:
                       refhvr_taxa[refhvr_id] = [tax]
            refdb_fh.close()
            #print "\n\n"
            #for refid in refhvr_taxa:
            #    print refid, refhvr_taxa[refid]
            #print "\n\ntaxa",refhvr_taxa
            # open names file and loop through the ids
            # foreach id find the taxa (from reference data)
            # and get consensus etc from taxonomy.py
            taxObjects = []
            names_fh = open(names_file,'r')
            tagtax_fh = open(tagtax_filename,'w')
            for line in names_fh:
                data  = line.strip().split("\t")
                id    = data[0]
                reads = data[1].split(",")
                #HKD9NCB01ANWGP_0	v6_AD520|AB086227|carbazole-degrading	0.017	60M
                #if id == 'HKD9NCB01CH4T3_0':
                #    print 'Found id HKD9NCB01CH4T3_0 '+ref_hits['HKD9NCB01CH4T3_0']
                #if id == 'HKD9NCB01C3J1P_0':
                #    print 'Found id HKD9NCB01C3J1P_0 '+ref_hits['HKD9NCB01C3J1P_0']
                        
                #if ( (exists $ref_hits{$id}) && ($ref_hits{$id}[0] ne 0) ) # 0 is a non-hit
                #if id in ref_hits and ref_hits[id][0] != 0:
                if id in ref_hits and ref_hits[id][0] != '0':
                    logger.debug('Found id from ref_hits: '+id+" "+' '.join(ref_hits[id]))
                    for ref in ref_hits[id]:
                        logger.debug('ref from ref_hits: '+ref)
                        for tax in refhvr_taxa[ref]:
                            
                            # get taxonomy object from taxonomy.py
                            taxObj = Taxonomy(tax)
                            taxObjects.append( taxObj )
                            #taxon = taxObj.taxstring()
                            #rank =  taxObj.depth()
                    taxReturn = consensus( taxObjects, C.majority );
                    # consensus rank and taxon
                    
                    # 0=taxObj, 1=winning vote, 2=minrank, 3=rankCounts, 4=maxPcts, 5=naPcts;
                    taxon = taxReturn[0].taxstring()
                    rank =  taxReturn[0].depth()
                    #print taxReturn[0], taxReturn[0].data,taxReturn[0].taxstring()
                    #while ($taxon =~ /;Unassigned$/) {$taxon =~ s/;Unassigned$//;}
                    while taxon.split(';')[-1] == 'Unassigned':
                        taxon = ';'.join(taxon.split(';')[:-1])
                    if not taxon:
                        taxon = "Unknown"
                    # don't trust below domain if too far from the hit
                    #print 'distances',distances[id],max_distance
                    #print rank, taxon
                    if float(distances[id]) >= float(max_distance):
                        #taxon =~ s/;.*$//;
                        taxon = taxon.split(';')[0]
                        rank = "domain"
                    for r in reads:
                        #"\t", $r, $taxon, $rank, scalar @taxObjects, $taxReturn[1], $taxReturn[2], $taxReturn[3], $taxReturn[4], $taxReturn[5], $distances{$id}) . "\n";                   
                        
                        tagtax_fh.write(r+"\t"+taxon+"\t"+rank+"\t"+str(len(taxObjects))+"\t"+taxReturn[1]+"\t"+taxReturn[2]+"\t"+taxReturn[3]+"\t"+taxReturn[4]+"\t"+taxReturn[5]+"\t"+str(distances[id])+"\n")
                else:
                    for r in reads:
                        tagtax_fh.write(r+"\tUnknown\tNA\t0\t0\tNA\t0;0;0;0;0;0;0;0\t0;0;0;0;0;0;0;0\t100;100;100;100;100;100;100;100\t1\n")
            tagtax_fh.close()
        logger.info("Finished gast2tax")
        print "Finished gast2tax"
