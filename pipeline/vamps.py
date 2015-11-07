import subprocess
import sys, os,stat
import time
import shutil
import datetime
from pipeline.pipelinelogging import logger
import constants as C
from pipeline.utils import Dirs, PipelneUtils
#from pipeline.db_upload import MyConnection
import ConMySQL

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
        
        
        if self.runobj.vamps_user_upload:
            site = self.runobj.site
            if self.runobj.mobedac:
                dir_prefix = 'mobedac_' + self.runobj.user + '_' + self.runobj.run
            else:
                dir_prefix = self.runobj.user + '_' + self.runobj.run
        else:
            site = ''
            dir_prefix = self.runobj.run
            
        dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, site = site) 
        
        
        
        
        #self.basedir = self.runobj.output_dir
        #self.outdir  = os.path.join(self.runobj.output_dir,self.prefix)
        self.idx_keys = idx_keys
        if self.runobj.vamps_user_upload:
            self.idx_keys = [self.runobj.user+self.runobj.run]
            self.iterator = self.runobj.datasets
        else:
            self.idx_keys = idx_keys
            self.iterator = self.idx_keys
        self.use_cluster = self.runobj.use_cluster
        
        
        os.environ['SGE_ROOT']='/opt/sge'
        os.environ['SGE_CELL']='grendel'
        path = os.environ['PATH']
        os.environ['PATH'] = path + ':/opt/sge/bin/lx24-amd64'
        #First step is to check for (or create via mothur)
        # a uniques fasta file and a names file 
        # one for each dataset.
        # If we are here from a vamps gast process
        # then there should be just one dataset to gast
        # but if MBL pipe then many datasets are prbably involved.
        # 1) clustergast
        # 2) gast cleanup
        # 3) gast2tax
        
        
        self.global_gast_dir = dirs.check_dir(dirs.gast_dir)
            
        
        if not os.path.exists(self.global_gast_dir):
            sys.exit("Could not find global gast dir: "+self.global_gast_dir)
        
        if self.runobj.site == 'vamps':
            db_host    = 'vampsdb'
            db_name    = 'vamps'
            db_home = '/xraid2-2/vampsweb/vamps/'
        else:
            db_host    = 'bpcweb7'
            db_name    = 'vamps'
            db_home = '/xraid2-2/vampsweb/vampsdev/'
        obj=ConMySQL.New(db_host, db_name, db_home)
        self.conn = obj.get_conn()    
        

    
    def create_vamps_files(self):
        """
        Creates the vamps files in the gast directory ready for loading into vampsdb
        """
        for key in self.iterator:
            (files, ds_count, gast_dir) = self.gather_files_per_key(key)
            if ds_count:
                 (tax_collector, read_id_lookup) = self.taxonomy(key, ds_count, files)
                 refid_collector = self.sequences(key, tax_collector, read_id_lookup, files)   
                 #self.exports(key, refid_collector, tax_collector, read_id_lookup, files)
                 self.projects(key, ds_count, files)
                 self.info(key, files)
            else:
                print "no tagtax file found or no dataset_count -- continuing to next dataset..."
            
    def load_vamps_db(self):
        """
        Loads files into vampsdb tables.
        """

        ret_value=''
        for key in self.iterator:
            (files, ds_count, gast_dir) = self.gather_files_per_key(key)
            if ds_count:
                if self.runobj.load_vamps_database:
                     ret_value = self.load_database(key, gast_dir, files)
                     if ret_value[:5] == 'ERROR':
                        return ret_value
            else:
                print "no tagtax file found or no dataset_count -- continuing to next dataset..."
        return 'SUCCESS'

    def gather_files_per_key(self, key):
    
        file_collector={}
        out_gast_dir = os.path.join(self.global_gast_dir,key)  #directory
        file_collector['gast_concat_file'] = os.path.join(out_gast_dir,'gast_concat')
        file_collector['tagtax_file'] = os.path.join(out_gast_dir,'tagtax_terse')
        if not os.path.exists(file_collector['gast_concat_file']):
                logger.warning("Could not find gast_concat_file file: "+file_collector['gast_concat_file'])
                 
        if not os.path.exists(file_collector['tagtax_file']):
            logger.warning("Could not find tagtax_file file: "+file_collector['tagtax_file'])
        #print key,self.runobj.platform
        
        if self.runobj.vamps_user_upload:
            
            file_collector['unique_file'] = os.path.join(out_gast_dir,'unique.fa')
            if self.runobj.project[:3] == 'MBE':   # Mobedac: use uniques in export
                file_collector['original_fa_file'] = os.path.join(out_gast_dir,'unique.fa')
            else:
                file_collector['original_fa_file'] = os.path.join(out_gast_dir,'fasta.fa')

            if self.runobj.fasta_file:
                grep_cmd = ['grep','-c','>',self.runobj.fasta_file]
            else:
                grep_cmd = ['grep','-c','>',file_collector['unique_file']]
        else:
            if self.runobj.platform == 'illumina':
                
                #unique_file = os.path.join(self.basedir,C.gast_dir),'unique.fa')
                reads_dir = dirs.check_dir(dirs.reads_overlap_dir)
                file_prefix = self.runobj.samples[key].file_prefix
                file_collector['unique_file'] = os.path.join(reads_dir,file_prefix+"-PERFECT_reads.fa.unique")
                # ANNA What is the correct file here:
                file_collector['original_fa_file'] = os.path.join(reads_dir,file_prefix+"-PERFECT_reads.fa.unique")
                grep_cmd = ['grep','-c','>',file_collector['unique_file']]
            elif self.runobj.platform == '454':
                pass
            else:
                sys.exit("no usable platform found")
            
        if not os.path.exists(file_collector['unique_file']):
            logger.error("Could not find unique_file: "+file_collector['unique_file'])
            
        # get dataset_count here from unique_file
        # the dataset_count should be from the non-unique file
        # but if we don't have that must use uniques
        
        try:
            dataset_count = subprocess.check_output(grep_cmd).strip()
        except:
            dataset_count = 0
        print key,": Sequence Count", dataset_count
        
        # output files to be created:            
        file_collector['taxes_file']                = os.path.join(out_gast_dir,'vamps_data_cube_uploads.txt')
        file_collector['summed_taxes_file']         = os.path.join(out_gast_dir,'vamps_junk_data_cube_pipe.txt')
        file_collector['distinct_taxes_file']       = os.path.join(out_gast_dir,'vamps_taxonomy_pipe.txt')
        file_collector['sequences_file']            = os.path.join(out_gast_dir,'vamps_sequences_pipe.txt')
        file_collector['export_file']               = os.path.join(out_gast_dir,'vamps_export_pipe.txt')
        file_collector['projects_datasets_file']    = os.path.join(out_gast_dir,'vamps_projects_datasets_pipe.txt')
        file_collector['project_info_file']         = os.path.join(out_gast_dir,'vamps_projects_info_pipe.txt')
    
    
        return (file_collector, dataset_count, out_gast_dir)
    

    
    
    def taxonomy(self,key, dataset_count, file_collector):
        """
        fill vamps_data_cube, vamps_junk_data_cube and vamps_taxonomy files
        """
        logger.info("Starting vamps_upload: taxonomy")
        print "Starting vamps_upload: taxonomy"
        # SUMMED create a look-up
        if self.runobj.vamps_user_upload:
            project = self.runobj.project
            dataset = key
        else:
            if self.runobj.platform == 'illumina':
                project = self.runobj.samples[key].project
                dataset = self.runobj.samples[key].dataset
            elif self.runobj.platform == '454':
                pass
            else:
                pass
            
        project = project[0].capitalize() + project[1:]
        project_dataset = project+'--'+dataset
        taxa_lookup = {}
        read_id_lookup={}
        if os.path.exists(file_collector['tagtax_file']):
            for line in  open(file_collector['tagtax_file'],'r'):
                line = line.strip()
                items = line.split("\t")
                taxa = items[1]
                if taxa[-3:] == ';NA':
                    taxa = taxa[:-3]
                read_id=items[0]  # this is the whole id from mothur -> won't match the clean id 
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
        fh1 = open(file_collector['taxes_file'],'w')
        
        
        fh1.write("\t".join( ["HEADER","project", "dataset", "taxonomy", "superkingdom", 
                            "phylum", "class", "orderx", "family", "genus", "species", 
                            "strain", "rank", "knt", "frequency", "dataset_count", "classifier"]) + "\n")
        tax_collector={}
        summer=0
        for tax,knt in taxa_lookup.iteritems():
            print tax,knt
            summer += knt
            datarow = ['',project,dataset]
            
            taxa = tax.split(';')
            #if taxa[0] in C.domains:
            freq = float(knt) / int(dataset_count)
            rank = C.ranks[len(taxa)-1]
            for i in range(len(C.ranks)):                
                if len(taxa) <= i:
                    if C.ranks[i] == 'orderx':
                    	taxa.append("order_NA")
                    else:
                    	taxa.append(C.ranks[i] + "_NA")

            tax_collector[tax] = {}


            datarow.append(tax)
            datarow.append("\t".join(taxa))
            datarow.append(rank)
            datarow.append(str(knt))
            datarow.append(str(freq))
            datarow.append(dataset_count)
            datarow.append(self.runobj.classifier)
            
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
        fh2 = open(file_collector['summed_taxes_file'],'w')
        
        fh2.write("\t".join(["HEADER","taxonomy", "sum_tax_counts", "frequency", "dataset_count","rank", 
                            "project","dataset","project--dataset","classifier"] )+"\n")
        ranks_subarray = []
        rank_list_lookup = {}
        for i in range(0, len(C.ranks)): 
            ranks_subarray.append(C.ranks[i])
            ranks_list = ";".join(ranks_subarray) # i.e., superkingdom, phylum, class
            # open data_cube file again
            # taxes_file: data_cube_uploads
            for line in  open(file_collector['taxes_file'],'r'):
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
            taxa = tax.split(';')
            #if taxa[0] in C.domains:
            rank = len( taxa ) -1
            
            frequency = float(knt) / int(dataset_count)
            
            if len(tax) - len(''.join(taxa)) >= rank:
            
                datarow = ['']
                datarow.append(tax)
                datarow.append(str(knt))
                datarow.append(str(frequency))
                datarow.append(str(dataset_count))
                datarow.append(str(rank))
                datarow.append(project)
                datarow.append(dataset)
                datarow.append(project_dataset)
                datarow.append(self.runobj.classifier)
            
                w = "\t".join(datarow)
                #print w
                fh2.write(w+"\n")
                

        fh2.close()
        
        
                
        #
        # DISTINCT TAXONOMY
        #
        fh3 = open(file_collector['distinct_taxes_file'],'w')
        fh3.write("\t".join(["HEADER","taxon_string", "rank", "num_kids"] )+"\n")
        taxon_string_lookup={}
        for line in  open(file_collector['summed_taxes_file'],'r'):
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
            taxa = taxon_string.split(';')
            if taxa[0] in C.domains:
                rank = str(len(taxa)-1)
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
        
    def sequences(self,key,tax_collector, read_id_lookup, file_collector):
        """
        fill vamps_sequences.txt file
        
        """
        
        logger.info("Starting vamps_upload: sequences")
        print "Starting vamps_upload: sequences"
        if self.runobj.vamps_user_upload:
            project = self.runobj.project
            dataset = key
        else:
            if self.runobj.platform == 'illumina':
                project = self.runobj.samples[key].project
                dataset = self.runobj.samples[key].dataset
            elif self.runobj.platform == '454':
                pass
            else:
                pass
        
        project = project[0].capitalize() + project[1:]
        project_dataset = project+'--'+dataset
        # open gast_concat table to get the distances and the ferids
        refid_collector={}
        #if os.path.exists(gast_concat_file):
        for line in  open(file_collector['gast_concat_file'],'r'):
            line = line.strip()
            items=line.split()
            id=items[0]
            distance=items[1]
            refhvr_ids=items[2]
            refid_collector[id]={}
            refid_collector[id]['distance']=distance
            refid_collector[id]['refhvr_ids']=refhvr_ids
                
            
        
        
        fh = open(file_collector['sequences_file'],'w')
        fh.write("\t".join(["HEADER","project","dataset","taxonomy","refhvr_ids", "rank",
                            "seq_count","frequency","distance","read_id","project_dataset"] )+"\n")
        
        
        
        # open uniques fa file
        if os.path.exists(file_collector['unique_file']) and os.path.getsize(file_collector['unique_file']) > 0:
            f = FastaReader(file_collector['unique_file'])
            
            while f.next():
                datarow = ['']
                defline_items = f.id.split('|')
                id = defline_items[0].split()[0]                
                cnt = defline_items[-1].split(':')[1]
                seq = f.seq.upper()
                if id in read_id_lookup:
                    print 'FOUND TAX for sequences file'
                    tax = read_id_lookup[id]
                else: 
                    print 'NO TAX for sequences file'
                    tax = ''
                    
                if tax in tax_collector:
                    rank = tax_collector[tax]['rank']
                    #cnt = tax_collector[tax]['knt']
                    freq = tax_collector[tax]['freq']
                else:
                    rank = 'NA'
                    cnt  = 0
                    freq = 0
                    
                if id in refid_collector:
                    distance = refid_collector[id]['distance']
                    refhvr_ids = refid_collector[id]['refhvr_ids']
                else:
                    distance = '1.0'
                    refhvr_ids = '0'
                if not cnt:
                		cnt = 1
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
        return refid_collector
        
    def exports(self, key, refid_collector, tax_collector, read_id_lookup, file_collector):
        """
        fill vamps_exports.txt file
    
        """
        
        logger.info("Skipping VAMPS exports()")
        
    def projects(self, key, dataset_count, file_collector):
        """
        fill vamps_projects_datasets.txt file
        """
        logger.info("Starting vamps_upload: projects_datasets")
        if self.runobj.vamps_user_upload:
            project = self.runobj.project
            dataset = key
        else:
            if self.runobj.platform == 'illumina':
                project = self.runobj.samples[key].project
                dataset = self.runobj.samples[key].dataset
            elif self.runobj.platform == '454':
                pass
            else:
                pass
        project = project[0].capitalize() + project[1:]
        project_dataset = project+'--'+dataset
        date_trimmed = 'unknown'
        dataset_description = dataset
        dataset_count = str(dataset_count)
        has_tax = '1' # true
        fh = open(file_collector['projects_datasets_file'],'w')
        
        fh.write("\t".join(["HEADER","project","dataset","dataset_count","has_tax", "date_trimmed","dataset_info"] )+"\n")
        fh.write("\t"+"\t".join([project, dataset, dataset_count, has_tax, date_trimmed, dataset_description] )+"\n")
        
        fh.close()
        logger.info("Finishing VAMPS projects()")
        
        
    def info(self, key, file_collector):
        """
        fill vamps_project_info.txt file
        """
        logger.info("Starting vamps_upload: projects_info")
        print "Starting vamps_upload: projects_info"
        if self.runobj.vamps_user_upload:
            user = self.runobj.user
            project = self.runobj.project
            sample_source_id = self.runobj.env_source_id
        else:
            if self.runobj.platform == 'illumina':
                user = self.runobj.samples[key].data_owner
                project = self.runobj.samples[key].project
                sample_source_id = self.runobj.samples[key].env_sample_source_id
            elif self.runobj.platform == '454':
                pass
            else:
                pass
        project = project[0].capitalize() + project[1:]
        cursor = self.conn.cursor()
        
        query = "SELECT last_name,first_name,email,institution from vamps_auth where user='%s'" % (user)
        #data = myconn.execute_fetch_select(query)
        cursor.execute(query)
        data = cursor.fetchone()
        
        fh = open(file_collector['project_info_file'],'w')
        title="title"
        description='description'
        contact= data[1]+' '+data[0]
        email= data[2]
        institution= data[3]
        
        fh.write("\t".join(["HEADER","project","title","description","contact", "email","institution","user","env_source_id"] )+"\n")
        fh.write("\t"+"\t".join([project, title, description, contact, email, institution, user, sample_source_id] )+"\n")
        # if this project already exists in the db???
        # the next step should update the table rather than add new to the db
        
        fh.close()
        self.conn.commit()
        cursor.close()
        
        logger.info("Finishing VAMPS info()")


    def load_database(self, key, out_gast_dir, file_collector):
        """
        
        """
        logger.info("Starting load VAMPS data")
        print "Starting load VAMPS data"

        # USER: vamps_db_tables
        execute_error = False
        if self.runobj.vamps_user_upload:
            user = self.runobj.user
            project = self.runobj.project
            data_cube_table     = C.database_tables['vamps_user_uploads']['tax_dc_tbl']
            summed_cube_table   = C.database_tables['vamps_user_uploads']['tax_summed_tbl']
            taxonomy_table      = C.database_tables['vamps_user_uploads']['tax_tbl']
            sequences_table     = C.database_tables['vamps_user_uploads']['sequences_tbl']
            export_table        = C.database_tables['vamps_user_uploads']['export_tbl']
            datasets_table      = C.database_tables['vamps_user_uploads']['datasets_tbl']
            users_table         = C.database_tables['vamps_user_uploads']['users_tbl']
        else:    
            if self.runobj.platform == 'illumina':
                user = self.runobj[key].data_owner
                project = self.runobj.samples[key].project
                data_cube_table     = C.database_tables['vamps_mbl_origin']['tax_dc_tbl']
                summed_cube_table   = C.database_tables['vamps_mbl_origin']['tax_summed_tbl']
                taxonomy_table      = C.database_tables['vamps_mbl_origin']['tax_tbl']
                sequences_table     = C.database_tables['vamps_mbl_origin']['sequences_tbl']
                export_table        = C.database_tables['vamps_mbl_origin']['export_tbl']
                datasets_table      = C.database_tables['vamps_mbl_origin']['datasets_tbl']
                users_table         = C.database_tables['vamps_mbl_origin']['users_tbl']
            elif self.runobj.platform == '454':
                pass
            else:
                pass
        info_table          = C.database_tables['vamps_mbl_origin']['info_tbl']
        users_info_table    = C.database_tables['vamps_user_uploads']['info_tbl']    
        
        
        cursor = self.conn.cursor()
        
        
 
        #
        #  DATA_CUBE
        #
        print "loading "+key+": data_cube"
        if os.path.exists(file_collector['taxes_file']):
            for line in open(file_collector['taxes_file'],'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                qDataCube = "insert ignore into %s (project, dataset, taxon_string,superkingdom,phylum,class, orderx,family,genus,species,strain,\
                            rank,knt,frequency,dataset_count,classifier)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (data_cube_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],
                            line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15])
                #myconn.execute_no_fetch(qDataCube)
                #print qDataCube
                try:
                    rows_affected = cursor.execute(qDataCube)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit()  
        else:
            print "taxes file not found for dataset "+key
            execute_error=True
        if execute_error:
            return 'ERROR: loading data_cube table'
            
        #
        # SUMMED (JUNK) DATA_CUBE
        #
        print "loading "+key+": junk_data_cube"
        if os.path.exists(file_collector['summed_taxes_file']):
            for line in open(file_collector['summed_taxes_file'],'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab
                #taxonomy        sum_tax_counts  frequency	dataset_count   rank    project dataset project--dataset        classifier
                qSummedCube = "insert ignore into %s (taxon_string,knt, frequency, dataset_count, rank, project, dataset, project_dataset, classifier)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (summed_cube_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7], line[8])
                #myconn.execute_no_fetch(qSummedCube) 
                #print qSummedCube 
                try:
                    cursor.execute(qSummedCube)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit()  
        else:
            print "summed taxes file not found for dataset "+key
            execute_error=True
        if execute_error:
            return 'ERROR: loading summed(junk)_data_cube table'    

        #
        #  TAXONOMY
        #
        print "loading "+key+": taxonomy"
        if os.path.exists(file_collector['distinct_taxes_file']):
            for line in open(file_collector['distinct_taxes_file'],'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab    
                qTaxonomy = "insert ignore into %s (taxon_string,rank,num_kids)\
                            VALUES('%s','%s','%s')" \
                            % (taxonomy_table, line[0],line[1],line[2])
                #myconn.execute_no_fetch(qTaxonomy)
                try:
                    cursor.execute(qTaxonomy)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit()  
        else:
            print "distinct taxes file not found for dataset "+key
            execute_error=True
        if execute_error:
            return 'ERROR: loading taxonomy table'   

        #
        #  SEQUENCES
        #
        print "loading "+key+": sequences"
        if os.path.exists(file_collector['sequences_file']):
            for line in open(file_collector['sequences_file'],'r'):
                line = line.strip().split("\t")
                if line[0]=='HEADER':
                    continue
                #line = line[1:] # remove leading empty tab
                # project dataset taxonomy        refhvr_ids	rank    seq_count frequency  distance  read_id project_dataset    
                qSequences = "insert ignore into %s (sequence,project, dataset, taxonomy,refhvr_ids,rank,seq_count,frequency,distance,rep_id, project_dataset)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (sequences_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7], line[8],line[9],line[10])
                #myconn.execute_no_fetch(qSequences) 
                try:
                    cursor.execute(qSequences)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit()  
        else:
            print "sequences file not found for dataset "+key  
            execute_error=True
        if execute_error:
            return 'ERROR: loading sequences table'     

        #
        #  EXPORT
        #
#         print "loading "+key+": export"
#         if os.path.exists(file_collector['export_file']):
#             for line in open(file_collector['export_file'],'r'):
#                 line = line.strip().split("\t")
#                 if line[0]=='HEADER':
#                     continue
#                 #line = line[1:] # remove leading empty tab
#                 # t.read_id, t.project, t.dataset, g.refhvr_ids, x.distance, x.taxonomy, t.sequence, x.rank," " t.entry_date
#                 # project dataset taxonomy        refhvr_ids	rank    seq_count frequency  distance  read_id project_dataset    
#                 qExport = "insert ignore into %s (read_id, project, dataset, refhvr_ids, distance, taxonomy, sequence, rank, date_trimmed)\
#                             VALUES('%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
#                             % (export_table,
#                             line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7], line[8])
# 
#                 #myconn.execute_no_fetch(qExport) 
# 
#                 try:
#                     cursor.execute(qExport)
#                 except:
#                     self.conn.rollback()  
#                     execute_error=True
#                 else:
#                     self.conn.commit() 
#         else:
#             print "export file not found for dataset "+key
#             execute_error=True
#         if execute_error:
#             return 'ERROR: loading export table'    

        #
        #  PROJECTS_DATASETS
        #
        print "loading "+key+": projects_datasets"
        if os.path.exists(file_collector['projects_datasets_file']):
            for line in open(file_collector['projects_datasets_file'],'r'):
                line = line.strip().split("\t")
                # [1:]  # split and remove the leading 'zero'
                if line[0]=='HEADER':
                    continue
                
                qDatasets = "insert ignore into %s (project, dataset, dataset_count,has_tax,date_trimmed,dataset_info)\
                            VALUES('%s','%s','%s','%s','%s','%s')" \
                            % (datasets_table,
                            line[0],line[1],line[2],line[3],line[4],line[5])
                #myconn.execute_no_fetch(qDatasets) 
                try:
                    cursor.execute(qDatasets)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit() 
                
                qDatasets = "update %s set has_tax='1' where project='%s'" \
                            % (datasets_table, line[0])
                #myconn.execute_no_fetch(qDatasets)
                try:
                    cursor.execute(qDatasets)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit() 
        else:
            print "project_datasets file not found for dataset "+key 
            execute_error=True
        if execute_error:
            return 'ERROR: loading datasets table(s)'     

        #
        # INFO
        #
        print "loading "+key+": info"
        if os.path.exists(file_collector['project_info_file']):
            for line in open(file_collector['project_info_file'],'r'):
                line = line.strip().split("\t")
                #[1:]  # split on tab and remove the leading 'zero'
                if line[0]=='HEADER':
                    continue
                
                qInfo = "insert ignore into %s (project_name, title, description, contact, email, institution, user, env_source_id)\
                            VALUES('%s','%s','%s','%s','%s','%s','%s','%s')" \
                            % (users_info_table,
                            line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7])
                #myconn.execute_no_fetch(qInfo) 
                try:
                    cursor.execute(qInfo)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit() 
                
                qInfo = "update %s set has_tax='1' where project_name='%s'" \
                            % (users_info_table, line[0])
                #myconn.execute_no_fetch(qInfo) 
                try:
                    cursor.execute(qInfo)
                except:
                    self.conn.rollback()  
                    execute_error=True
                else:
                    self.conn.commit() 
        else:
            print "upload_info file not found for dataset "+key 
            execute_error=True
        if execute_error:
            return 'ERROR: loading info table'          

        #
        # USERS
        #
        print "loading users:"+key
        
        qUser = "insert ignore into %s (project, user)\
                    VALUES('%s','%s')" \
                    % (users_table, project, user)
        #myconn.execute_no_fetch(qUser) 
        try:
            cursor.execute(qUser)
        except:
            self.conn.rollback()  
            execute_error=True
        else:
            self.conn.commit()      
        if execute_error:
            return 'ERROR: loading info table'     
        logger.info("Finished load VAMPS data")
        
        #self.conn.commit()
        cursor.close()
        return 'SUCCESS'
