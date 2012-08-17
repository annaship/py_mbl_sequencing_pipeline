import subprocess
import sys, os,stat
import time
import shutil
from pipeline.pipelinelogging import logger
import constants as C
from pipeline.db_upload import MyConnection

class FastaReader:
    def __init__(self,file_name=None):
        self.file_name = file_name
        self.h = open(self.file_name)
        self.seq = ''
        self.id = None

    def next(self): 
        def read_id():
            return self.h.readline().strip()[1:]

        def read_seq():
            ret = ''
            while True:
                line = self.h.readline()
                
                while len(line) and not len(line.strip()):
                    # found empty line(s)
                    line = self.h.readline()
                
                if not len(line):
                    # EOF
                    break
                
                if line.startswith('>'):
                    # found new defline: move back to the start
                    self.h.seek(-len(line), os.SEEK_CUR)
                    break
                    
                else:
                    ret += line.strip()
                    
            return ret
        
        self.id = read_id()
        self.seq = read_seq()
        
        if self.id:
            return True  
            
class Vamps:
    """Uploads data to the VAMPS (or vampsdev) database"""
    Name = "VAMPS"
    def __init__(self, run_object = None, idx_keys=None):

        self.runobj 	 = run_object
        
        self.basedir = self.runobj.output_dir
        #self.outdir  = os.path.join(self.runobj.output_dir,self.prefix)
        self.idx_keys = idx_keys

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
        
        # already present:
        self.global_gast_dir = os.path.join(self.basedir,C.gast_dir)
        if not os.path.exists(self.global_gast_dir):
            sys.exit("Could not find global gast dir: "+self.global_gast_dir)
        
        for key in self.idx_keys:
            out_gast_dir = os.path.join(self.global_gast_dir,key)  #directory
            gast_concat_file = os.path.join(out_gast_dir,'gast_concat')
            if not os.path.exists(gast_concat_file):
                logger.warning("Could not find gast_concat_file file: "+gast_concat_file)
            #self.gast_file = os.path.join(self.out_gast_dir,self.prefix+'.gast')
            tagtax_file = os.path.join(out_gast_dir,'tagtax_terse')
            if not os.path.exists(tagtax_file):
                logger.warning("Could not find tagtax_file file: "+tagtax_file)
            
            if self.runobj.platform == 'vamps':
                unique_file = os.path.join(out_gast_dir,'unique.fa')
                if self.runobj.fasta_file:
                    grep_cmd = ['grep','-c','>',self.runobj.fasta_file]
                else:
                    grep_cmd = ['grep','-c','>',unique_file]
            elif self.runobj.platform == 'illumina':
                #unique_file = os.path.join(self.basedir,C.gast_dir),'unique.fa')
                reads_dir = os.path.join(self.basedir,C.illumina_reads_dir)
                file_prefix = self.runobj.samples[key].file_prefix
                unique_file = os.path.join(reads_dir,file_prefix+"-PERFECT_reads.fa.unique")
                grep_cmd = ['grep','-c','>',unique_file]
            else:
                sys.exit()
                
            if not os.path.exists(unique_file):
                logger("Could not find unique_file: "+unique_file)
                
            # get dataset_count here from unique_file
            # the dataset_count should be from the non-unique file
            # but if we don't have that must use uniques
            
            try:
                dataset_count = subprocess.check_output(grep_cmd).strip()
            except:
                dataset_count = 0
            print "Dataset Count", dataset_count
            
            # to be created:
            taxes_file              = os.path.join(out_gast_dir,'vamps_data_cube_uploads.txt')
            summed_taxes_file       = os.path.join(out_gast_dir,'vamps_junk_data_cube_pipe.txt')
            distinct_taxes_file     = os.path.join(out_gast_dir,'vamps_taxonomy_pipe.txt')
            sequences_file          = os.path.join(out_gast_dir,'vamps_sequences_pipe.txt')
            export_file             = os.path.join(out_gast_dir,'vamps_export_pipe.txt')
            projects_datasets_file  = os.path.join(out_gast_dir,'vamps_projects_datasets_pipe.txt')
            project_info_file       = os.path.join(out_gast_dir,'vamps_projects_info_pipe.txt')
            
            (tax_collector,read_id_lookup) = self.taxonomy(key, tagtax_file, dataset_count, taxes_file, summed_taxes_file, distinct_taxes_file)
            
            self.sequences(key, tax_collector, read_id_lookup, unique_file, gast_concat_file, sequences_file)   
            
            self.exports()
            self.projects(key, projects_datasets_file, dataset_count)
            self.info(key, project_info_file)
            if self.runobj.load_vamps_database:
                self.load_database(key, out_gast_dir, taxes_file, summed_taxes_file, distinct_taxes_file, sequences_file, projects_datasets_file, project_info_file)
                #self.load_taxonomy
                #self.load_sequences
                #self.load_export
                #self.load_projects
                
    def taxonomy(self,key,tagtax_file,dataset_count,taxes_file,summed_taxes_file,distinct_taxes_file):
        """
        fill vamps_data_cube, vamps_junk_data_cube and vamps_taxonomy files
        """
        logger.info("Starting vamps_upload: taxonomy")
        # SUMMED create a look-up
        if self.runobj.platform == 'vamps':
            project = self.runobj.project
            dataset = self.runobj.dataset
            project_dataset = project+'--'+dataset
        elif self.runobj.platform == 'illumina':
            project = self.runobj.samples[key].project
            dataset = self.runobj.samples[key].dataset
            project_dataset = project+'--'+dataset
            
        taxa_lookup = {}
        read_id_lookup={}
        if os.path.exists(tagtax_file):
            for line in  open(tagtax_file,'r'):
                line = line.strip()
                items = line.split()
                taxa = items[1]
                if taxa[-3:] == ';NA':
                    taxa = taxa[:-3]
                read_id=items[0]
                read_id_lookup[read_id]=taxa
                
                # the count here is the frequency of the taxon in the datasets
                if taxa in taxa_lookup:                
                    taxa_lookup[taxa] += 1 
                else:
                    taxa_lookup[taxa] = 1 
                      
        #  DATA CUBE TABLE
        # taxa_lookup: {'Unknown': 146, 'Bacteria': 11888, 'Bacteria;Chloroflexi': 101}
        # dataset_count is 3 (3 taxa in this dataset)
        # frequency is 3/144
        fh1 = open(taxes_file,'w')
        
        
        fh1.write("\t".join( ["HEADER","project", "dataset", "taxonomy", "superkingdom", 
                            "phylum", "class", "orderx", "family", "genus", "species", 
                            "strain", "rank", "knt", "frequency", "dataset_count", "classifier"]) + "\n")
        tax_collector={}
        summer=0
        for tax,knt in taxa_lookup.iteritems():
            #print tax,cnt
            summer += knt
            datarow = ['',project,dataset]
            
            taxes = tax.split(';')
            
            freq = float(knt) / int(dataset_count)
            rank = C.ranks[len(taxes)-1]
            for i in range(len(C.ranks)):                
                if len(taxes) <= i:
                    taxes.append(C.ranks[i] + "_NA")

            tax_collector[tax] = {}


            datarow.append(tax)
            datarow.append("\t".join(taxes))
            datarow.append(rank)
            datarow.append(str(knt))
            datarow.append(str(freq))
            datarow.append(dataset_count)
            datarow.append("GAST")
            
            w = "\t".join(datarow)
            #print w
            fh1.write(w+"\n")
           
            tax_collector[tax]['rank'] = rank
            tax_collector[tax]['knt'] = knt
            tax_collector[tax]['freq'] = freq
        
        fh1.close()
        
        #
        # SUMMED DATA CUBE TABLE
        #
        fh2 = open(summed_taxes_file,'w')
        
        fh2.write("\t".join(["HEADER","taxonomy", "sum_tax_counts", "frequency", "dataset_count","rank", 
                            "project","dataset","project--dataset","classifier"] )+"\n")
        ranks_subarray = []
        rank_list_lookup = {}
        for i in range(0, len(C.ranks)): 
            ranks_subarray.append(C.ranks[i])
            ranks_list = ";".join(ranks_subarray) # i.e., superkingdom, phylum, class
            # open data_cube file again
            # taxes_file: data_cube_uploads
            for line in  open(taxes_file,'r'):
                line = line.strip().split("\t")
                knt = line[12]
                taxon = line[2]
                if line[0] == 'HEADER':
                    continue
                if taxon in tax_collector:
                    knt = tax_collector[taxon]['knt']
                else:
                    print 'ERROR tax not found in tax_collector: assigning zero'
                    knt = 0
                idx = len(ranks_subarray)
                l=[]
                for k in range(3,idx+3):                    
                    l.append(line[k])
                tax = ';'.join(l)
                #print 'rl tax',ranks_list,tax
                
                
                if tax in rank_list_lookup:
                    rank_list_lookup[tax] += knt
                else:
                    rank_list_lookup[tax] = knt
                    
                
          
        for tax,knt in rank_list_lookup.iteritems():
            
           
            
            #print 'tax2',tax
            rank = len( tax.split(';') ) -1
            frequency = float(knt) / int(dataset_count)
            
            
            if len(tax) - len(''.join(tax.split(';'))) >= rank:
            
                datarow = ['']
                datarow.append(tax)
                datarow.append(str(knt))
                datarow.append(str(frequency))
                datarow.append(str(dataset_count))
                datarow.append(str(rank))
                datarow.append(project)
                datarow.append(dataset)
                datarow.append(project_dataset)
                datarow.append("GAST")
            
                w = "\t".join(datarow)
                #print w
                fh2.write(w+"\n")
                

        fh2.close()
        
        
        
        
        
        #
        # DISTINCT TAXONOMY
        #
        fh3 = open(distinct_taxes_file,'w')
        fh3.write("\t".join(["HEADER","taxon_string", "rank", "num_kids"] )+"\n")
        taxon_string_lookup={}
        for line in  open(summed_taxes_file,'r'):
            if line.split()[0] == 'HEADER':
                continue
            items = line.strip().split()            
            taxon_string = items[0]
            #print taxon_string
            if taxon_string in taxon_string_lookup:
                taxon_string_lookup[taxon_string] += 1
            else:
                taxon_string_lookup[taxon_string] = 1
        
        for taxon_string,v in taxon_string_lookup.iteritems():
            datarow = ['']
            datarow.append(taxon_string)
            rank = str(len(taxon_string.split(';'))-1)
            datarow.append(rank)
            if rank==7 or taxon_string[-3:]=='_NA':
                num_kids = '0'
            else:
                num_kids = '1'
            datarow.append(num_kids)
            w = "\t".join(datarow)
            #print 'w',w
            fh3.write(w+"\n")
        fh3.close()
        
        return (tax_collector,read_id_lookup)
        
    def sequences(self,key,tax_collector, read_id_lookup, unique_file, gast_concat_file, sequences_file):
        """
        fill vamps_sequences.txt file
        """
        
        logger.info("Starting vamps_upload: sequences")
        if self.runobj.platform == 'vamps':
            project = self.runobj.project
            dataset = self.runobj.dataset
            project_dataset = project+'--'+dataset
        elif self.runobj.platform == 'illumina':
            project = self.runobj.samples[key].project
            dataset = self.runobj.samples[key].dataset
            project_dataset = project+'--'+dataset
            
        # open gast_concat table to get the distances and the ferids
        refid_collector={}
        #if os.path.exists(gast_concat_file):
        for line in  open(gast_concat_file,'r'):
            line = line.strip()
            items=line.split()
            id=items[0]
            distance=items[1]
            refhvr_ids=items[2]
            refid_collector[id]={}
            refid_collector[id]['distance']=distance
            refid_collector[id]['refhvr_ids']=refhvr_ids
                
            
        
        
        fh = open(sequences_file,'w')
        fh.write("\t".join(["HEADER","project","dataset","taxonomy","refhvr_ids", "rank",
                            "seq_count","frequency","distance","read_id","project_dataset"] )+"\n")
        
        
        
        # open uniques fa file
        if os.path.exists(unique_file) and os.path.getsize(unique_file) > 0:
            f = FastaReader(unique_file)
            
            while f.next():
                datarow = ['']
                id = f.id.split('|')[0]
                seq = f.seq
                tax = read_id_lookup[id]
                rank = tax_collector[tax]['rank']
                cnt = tax_collector[tax]['knt']
                freq = tax_collector[tax]['freq']
                if id in refid_collector:
                    distance = refid_collector[id]['distance']
                    refhvr_ids = refid_collector[id]['refhvr_ids']
                else:
                    distance = '1.0'
                    refhvr_ids = '0'
                
                datarow.append(seq)
                datarow.append(project)
                datarow.append(dataset)
                datarow.append(tax)
                datarow.append(refhvr_ids)
                datarow.append(rank)
                datarow.append(str(cnt))
                datarow.append(str(freq))
                datarow.append(distance)
                datarow.append(id)
                datarow.append(project_dataset)
                w = "\t".join(datarow)
                #print 'w',w
                fh.write(w+"\n")
                
                
            fh.close()
        logger.info("")
        
    def exports(self):
        """
        fill vamps_exports.txt file
        """
        logger.info("Starting vamps_upload: exports")
        print "TODO: upload_vamps 5- exports"
        logger.info("Finishing VAMPS exports()")
        
    def projects(self, key, projects_datasets_file, dataset_count):
        """
        fill vamps_projects_datasets.txt file
        """
        logger.info("Starting vamps_upload: projects_datasets")
        if self.runobj.platform == 'vamps':
            project = self.runobj.project
            dataset = self.runobj.dataset
            project_dataset = project+'--'+dataset
        elif self.runobj.platform == 'illumina':
            project = self.runobj.samples[key].project
            dataset = self.runobj.samples[key].dataset
            project_dataset = project+'--'+dataset
        date_trimmed = 'unknown'
        dataset_description = dataset
        dataset_count = str(dataset_count)
        has_tax = '1' # true
        fh = open(projects_datasets_file,'w')
        
        fh.write("\t".join(["HEADER","project","dataset","dataset_count","has_tax", "date_trimmed","dataset_info"] )+"\n")
        fh.write("\t"+"\t".join([project, dataset, dataset_count, has_tax, date_trimmed, dataset_description] )+"\n")
        
        fh.close()
        logger.info("Finishing VAMPS projects()")
        
        
    def info(self, key,project_info_file):
        """
        fill vamps_project_info.txt file
        """
        logger.info("Starting vamps_upload: projects_info")
        if self.runobj.platform == 'vamps':
            user = self.runobj.user
            project = self.runobj.project
            sample_source = self.runobj.env_source_id
        elif self.runobj.platform == 'illumina':
            user = self.runobj.samples[key].data_owner
            project = self.runobj.samples[key].project
            sample_source = self.runobj.samples[key].env_sample_source
        if self.runobj.site == 'vamps':
            db_host    = 'vampsdb'
            db_name    = 'vamps'
        else:
            db_host    = 'vampsdev'
            db_name    = 'vamps'
        myconn = MyConnection(host=db_host, db=db_name)
        query = "SELECT last_name,first_name,email,institution from vamps_auth where user='%s'" % (user)
        data = myconn.execute_fetch_select(query)

        
        fh = open(project_info_file,'w')
         
        title="title"
        description='description'
        contact= data[0][1]+' '+data[0][0]
        email= data[0][2]
        institution= data[0][3]
        
        fh.write("\t".join(["HEADER","project","title","description","contact", "email","institution","user","env_source_id"] )+"\n")
        fh.write("\t"+"\t".join([project, title, description, contact, email, institution, user, sample_source] )+"\n")
        # if this project already exists in the db???
        # the next step should update the table rather than add new to the db
        
        fh.close()
        logger.info("Finishing VAMPS info()")

    def load_database(self,key,out_gast_dir,taxes_file,summed_taxes_file,distinct_taxes_file,sequences_file,projects_datasets_file,project_info_file):
        """
        
        """
        logger.info("Starting load VAMPS data")
#         self.taxes_file = os.path.join(self.outdir,'vamps_data_cube_uploads.txt')
#         self.summed_taxes_file = os.path.join(self.outdir,'vamps_junk_data_cube_pipe.txt')
#         self.distinct_taxes_file = os.path.join(self.outdir,'vamps_taxonomy_pipe.txt')
#         self.sequences_file = os.path.join(self.outdir,'vamps_sequences_pipe.txt')
#         self.export_file = os.path.join(self.outdir,'vamps_export_pipe.txt')
#         self.projects_datasets_file = os.path.join(self.outdir,'vamps_projects_datasets_pipe.txt')
#         self.projects_info_file = os.path.join(self.outdir,'vamps_projects_info_pipe.txt')
        # USER: vamps_db_tables
        if self.runobj.platform == 'vamps':
            user = self.runobj.user
            project = self.runobj.project
            data_cube_table     = C.database_tables['vamps_user_uploads']['tax_dc_tbl']
            summed_cube_table   = C.database_tables['vamps_user_uploads']['tax_summed_tbl']
            taxonomy_table      = C.database_tables['vamps_user_uploads']['tax_tbl']
            sequences_table     = C.database_tables['vamps_user_uploads']['sequences_tbl']
            export_table        = C.database_tables['vamps_user_uploads']['export_tbl']
            datasets_table      = C.database_tables['vamps_user_uploads']['datasets_tbl']
            users_table         = C.database_tables['vamps_user_uploads']['users_tbl']
            
        elif self.runobj.platform == 'illumina':
            user = self.runobj[key].data_owner
            project = self.runobj.samples[key].project
            data_cube_table     = C.database_tables['vamps_mbl_origin']['tax_dc_tbl']
            summed_cube_table   = C.database_tables['vamps_mbl_origin']['tax_summed_tbl']
            taxonomy_table      = C.database_tables['vamps_mbl_origin']['tax_tbl']
            sequences_table     = C.database_tables['vamps_mbl_origin']['sequences_tbl']
            export_table        = C.database_tables['vamps_mbl_origin']['export_tbl']
            datasets_table      = C.database_tables['vamps_mbl_origin']['datasets_tbl']
            users_table         = C.database_tables['vamps_mbl_origin']['users_tbl']
            
        info_table          = C.database_tables['vamps_mbl_origin']['info_tbl']
        users_info_table    = C.database_tables['vamps_user_uploads']['info_tbl']    
        # We only have a single project and dataset here:
        # if the project is new  then we add the data to the upload_info and projects_datasets_pipe table
        # but if the project is not new:
        #   check if the existing project belongs to the user
        #   if it does then UPDATE the line in upload_info table and add line to projects_datasets_pipe table
        #       (maybe check if dataset already exists and die if yes)
        #   if the existing project doesn't belong to the owner then die with a warning to change project name
        #      (or maybe change the name by adding _user)

        
        if self.runobj.site == 'vamps':
            db_host    = 'vampsdb'
            db_name    = 'vamps'
        else:
            db_host    = 'vampsdev'
            db_name    = 'vamps'
            
        myconn = MyConnection(host=db_host, db=db_name)
        query = "SELECT project_name from %s where project_name='%s' \
                    UNION \
                 SELECT project_name from %s where project_name='%s' \
                 " % (users_info_table, project, info_table, project)
                 
        data = myconn.execute_fetch_select(query)
        if data:
            logger.info("found this project "+data[0][0]+" Exiting")
            sys.exit("Duplicate project name found; Canceling upload to database but your GASTed data are here: "+ out_gast_dir)
        else:
            # project is unknown in database - continue
 
 
            #
            #  DATA_CUBE
            #
            for line in open(taxes_file,'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab
                
                qDataCube = "insert ignore into %s (project, dataset, taxon_string,superkingdom,phylum,class,\
                                            orderx,family,genus,species,strain,rank,knt,frequency,dataset_count,classifier)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (data_cube_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],
                            line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15])
                myconn.execute_no_fetch(qDataCube)
                
            #
            # SUMMED (JUNK) DATA_CUBE
            #
            for line in open(summed_taxes_file,'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab
                #taxonomy        sum_tax_counts  frequency	dataset_count   rank    project dataset project--dataset        classifier
                qSummedCube = "insert ignore into %s (taxon_string,knt, frequency, dataset_count, rank, project, dataset, project_dataset, classifier)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (summed_cube_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7], line[8])
                myconn.execute_no_fetch(qSummedCube)    
                    
                    
            #
            #  TAXONOMY
            #
            for line in open(distinct_taxes_file,'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab    
                qTaxonomy = "insert ignore into %s (taxon_string,rank,num_kids)\
                            VALUES('%s','%s','%s')" \
                            % (taxonomy_table, line[0],line[1],line[2])
                myconn.execute_no_fetch(qTaxonomy)        
            
            #
            #  SEQUENCES
            #
            for line in open(sequences_file,'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab
                # project dataset taxonomy        refhvr_ids	rank    seq_count frequency  distance  read_id project_dataset    
                qSequences = "insert ignore into %s (sequence,project, dataset, taxonomy,refhvr_ids,rank,seq_count,frequency,distance,rep_id, project_dataset)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (sequences_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7], line[8],line[9],line[10])
                myconn.execute_no_fetch(qSequences)        
            #
            #  PROJECTS_DATASETS
            #
            for line in open(projects_datasets_file,'r'):
                line = line.strip().split("\t")
                # [1:]  # split and remove the leading 'zero'
                if line[0]=='HEADER':
                    continue
                
                qDatasets = "insert ignore into %s (project, dataset, dataset_count,has_tax,date_trimmed,dataset_info)\
                            VALUES('%s','%s','%s','%s','%s','%s')" \
                            % (datasets_table,
                            line[0],line[1],line[2],line[3],line[4],line[5])
                myconn.execute_no_fetch(qDatasets) 
            
            #
            # INFO
            #
            for line in open(project_info_file,'r'):
                line = line.strip().split("\t")
                #[1:]  # split on tab and remove the leading 'zero'
                if line[0]=='HEADER':
                    continue
                
                qInfo = "insert ignore into %s (project_name, title, description, contact, email, institution, user, env_source_id)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (users_info_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7])
                myconn.execute_no_fetch(qInfo) 
                
            #
            # USERS
            #
                
            qUser = "insert ignore into %s (project, user)\
                        VALUES('%s','%s')" \
                        % (users_table, project, user)
            myconn.execute_no_fetch(qUser) 
            
            
        
        
        logger.info("Finished load VAMPS data")
        
        
        