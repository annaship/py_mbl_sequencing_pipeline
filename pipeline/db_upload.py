import sys

import os
from pipeline.get_ini import readCSV
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils

try:
    import MySQLdb
except:
    sys.exit("""
    MySQLdb ERROR
      To load the correct module, try running these commands before running the pipeline:
      
source /xraid/bioware/Modules/etc/profile.modules
module load bioware

    """)
class MyConnection:
    """
    Connection to env454
    Takes parameters from ~/.my.cnf, default host = "vampsdev", db="test"
    if different use my_conn = MyConnection(host, db)
    """
    def __init__(self, host="vampsdev", db="test"):
        self.conn   = None
        self.cursor = None
        self.rows   = 0
        self.new_id = None
        self.lastrowid = None
                
        try:
            print "=" * 40
            print "host = " + str(host) + ", db = "  + str(db)
            print "=" * 40

            self.conn   = MySQLdb.connect(host=host, db=db, read_default_file="~/.my.cnf")
            self.cursor = self.conn.cursor()
                   
        except MySQLdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            raise
        except:                       # catch everything
            print "Unexpected:"         # handle unexpected exceptions
            print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
            raise                       # re-throw caught exception   

    def execute_fetch_select(self, sql):
        if self.cursor:
            self.cursor.execute(sql)
            res = self.cursor.fetchall ()
            return res

    def execute_no_fetch(self, sql):
        if self.cursor:
            self.cursor.execute(sql)
            self.conn.commit()
#            if (self.conn.affected_rows()):
#            print dir(self.cursor)
            return self.cursor.lastrowid
#        logger.debug("rows = "  + str(self.rows))
 

class dbUpload:
    """db upload methods"""
    Name = "dbUpload"
    """
    TODO: add tests and test case
    TODO: change hardcoded values to args: 
        self.sequence_table_name = "sequence_ill", 
        self.sequence_field_name = "sequence_comp"  
    TODO: generalize all bulk uploads and all inserts? to not copy and paste
    TODO: add refssu_id
    TODO: change csv validaton for new fields
    Order:
        # put_run_info
        # insert_seq()
        # insert_pdr_info()
        # gast
        # insert_taxonomy()
        # insert_sequence_uniq_info_ill()

    """
    def __init__(self, runobj = None):
        self.runobj      = runobj
        self.rundate     = self.runobj.run
        self.use_cluster = 1       
#        if self.runobj.vamps_user_upload:
#            site       = self.runobj.site
#            dir_prefix = self.runobj.user + '_' + self.runobj.run
#        else:
#            site = ''
#            dir_prefix = self.runobj.run         
#        dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, site = site)

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
 
        
        self.analysis_dir = dirs.check_dir(dirs.analysis_dir)
        self.fasta_dir    = dirs.check_dir(dirs.reads_overlap_dir)
        self.gast_dir     = dirs.check_dir(dirs.gast_dir)

        host_name     = runobj.database_host
        database_name = runobj.database_name
        
        self.filenames   = []
#        self.my_conn     = MyConnection(host = 'newbpcdb2', db="env454")
        self.my_conn     = MyConnection()    
        self.sequence_table_name = "sequence_ill" 
        self.sequence_field_name = "sequence_comp" 
        self.my_csv      = None

        self.unique_file_counts = dirs.unique_file_counts
        dirs.delete_file(self.unique_file_counts)
        self.seq_id_dict = {}
        self.tax_id_dict = {}
        self.run_id      = None
        
#        self.refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
   
    def get_fasta_file_names(self):
        fa_files = []
        pipelne_utils   = PipelneUtils()
        files = pipelne_utils.get_all_files(self.fasta_dir)
        for full_name in files.keys():    
            if (files[full_name][1] == ".unique") and ((files[full_name][0].split(".")[-1].strip() == "fa") or (files[full_name][0].split("_")[-1] == "MERGED")):
                fa_files.append(full_name)
        return fa_files
        
    def get_run_info_ill_id(self, filename_base):
        
        my_sql = """SELECT run_info_ill_id FROM run_info_ill 
                    JOIN run using(run_id)
                    WHERE file_prefix = '%s'
                    and run = '%s'
        """ % (filename_base, self.rundate)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])
        
    def insert_seq(self, sequences):
        query_tmpl = "INSERT IGNORE INTO %s (%s) VALUES (COMPRESS(%s))"
        val_tmpl   = "'%s'"
        my_sql     = query_tmpl % (self.sequence_table_name, self.sequence_field_name, ')), (COMPRESS('.join([val_tmpl % key for key in sequences]))
        seq_id     = self.my_conn.execute_no_fetch(my_sql)
        print "sequences in file: %s" % (len(sequences))
        return seq_id
        
    def get_seq_id_dict(self, sequences):
        id_name    = self.sequence_table_name + "_id" 
        query_tmpl = """SELECT %s, uncompress(%s) FROM %s WHERE %s in (COMPRESS(%s))"""
        val_tmpl   = "'%s'"
        my_sql     = query_tmpl % (id_name, self.sequence_field_name, self.sequence_table_name, self.sequence_field_name, '), COMPRESS('.join([val_tmpl % key for key in sequences]))
        res        = self.my_conn.execute_fetch_select(my_sql)
        self.seq_id_dict = dict((y, int(x)) for x, y in res)

    def get_id(self, table_name, value):
        id_name = table_name + '_id'
        my_sql  = """SELECT %s FROM %s WHERE %s = '%s'""" % (id_name, table_name, table_name, value)
        res     = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])         
            
    def get_sequence_id(self, seq):
        my_sql = """SELECT sequence_ill_id FROM sequence_ill WHERE COMPRESS('%s') = sequence_comp""" % (seq)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])     
    
    def insert_pdr_info(self, fasta, run_info_ill_id):
        # ------- insert sequence info per run/project/dataset --------
        sequence_ill_id = self.seq_id_dict[fasta.seq]

        seq_count       = int(fasta.id.split('|')[-1].split(':')[-1])
#        print run_info_ill_id, sequence_ill_id, seq_count
        my_sql          = """INSERT IGNORE INTO sequence_pdr_info_ill (run_info_ill_id, sequence_ill_id, seq_count) 
                             VALUES (%s, %s, %s)""" % (run_info_ill_id, sequence_ill_id, seq_count)
        res_id = self.my_conn.execute_no_fetch(my_sql)
        return res_id
 
    def get_gasta_result(self, filename):
        gast_file_name = os.path.join(self.gast_dir, filename + '.gast')
        try:
            with open(gast_file_name) as fd:
                gast_dict = dict([(l.split("\t")[0], l.split("\t")[1:]) for l in fd])    
            return gast_dict
        except IOError, e:
#            print dir(e)
#['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__getitem__', '__getslice__', '__hash__', '__init__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', '__unicode__', 'args', 'errno', 'filename', 'message', 'strerror']
#            print "errno = %s" % e.errno
            logger.debug("errno = %s" % e.errno)
            if e.errno == 2:
                # suppress "No such file or directory" error
                pass            
        except OSError, e:
                # reraise the exception, as it's an unexpected error
                raise
        
    """TODO: why run_info_ill
    1    1529    2164    8    6951    2411    83            19    JV    JV    GCCTAA    0    230    6_FP1BermC_6_14_10_CGCTC    101    23
    """
    def insert_taxonomy(self, fasta, gast_dict):
        if gast_dict:
            (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast_dict[fasta.id]
            "if we already had this taxonomy in this run, just skip it"
            if taxonomy in self.tax_id_dict:
                next
            else:
                tax_id = self.get_id("taxonomy", taxonomy)
                if tax_id:
                    self.tax_id_dict[taxonomy] = tax_id
                else:
                    my_sql = """INSERT IGNORE INTO taxonomy (taxonomy) VALUES ('%s')""" % (taxonomy.rstrip())
                    tax_id = self.my_conn.execute_no_fetch(my_sql)
                    self.tax_id_dict[taxonomy] = tax_id
                return tax_id
            
    def insert_sequence_uniq_info_ill(self, fasta, gast_dict):
        if gast_dict:
            (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast_dict[fasta.id]
            sequence_ill_id = self.seq_id_dict[fasta.seq]
            if taxonomy in self.tax_id_dict:
                try:
                    taxonomy_id = self.tax_id_dict[taxonomy] 
                except Exception, e:
                    logger.debug("Error = %s" % e)
                    raise
            my_sql = """INSERT IGNORE INTO sequence_uniq_info_ill (sequence_ill_id, taxonomy_id, gast_distance, refssu_count, rank_id, refhvr_ids) VALUES
                   (
                    %s,
                    %s,
                    '%s',
                    '%s',
                    (SELECT rank_id FROM rank WHERE rank = '%s'),
                    '%s'                
                   )
                   """ % (sequence_ill_id, taxonomy_id, distance, refssu_count, rank, refhvr_ids.rstrip())
            res_id = self.my_conn.execute_no_fetch(my_sql)
            return res_id
    
    def put_run_info(self, content = None):

        values = list(set([run_key.split('_')[1] for run_key in self.runobj.run_keys]))
        self.insert_bulk_data('run_key', values)
        values = list(set([self.runobj.samples[key].dna_region for key in self.runobj.samples]))
        self.insert_bulk_data('dna_region', values)
        self.insert_rundate()
        self.insert_test_contact()
       
#        --------- indiv_inserts --------- 

#        for k, v in content.items():
#            print k, v['dataset'], v
        for key in self.runobj.samples:
            value = self.runobj.samples[key]
#            print dict(run.samples[key])
#            print key, value.dataset
            self.get_contact_v_info()
#            self.insert_contact()
            contact_id = self.get_contact_id(value.data_owner)
            self.insert_project(value, contact_id)
            self.insert_dataset(value) 

            self.insert_run_info(value)
#            self.insert_primer()
#        return content[1].keys()

    def insert_bulk_data(self, key, values):
        query_tmpl = "INSERT IGNORE INTO %s (%s) VALUES (%s)"
        val_tmpl   = "'%s'"
        my_sql     = query_tmpl % (key, key, '), ('.join([val_tmpl % key for key in values]))
        self.my_conn.execute_no_fetch(my_sql)
    
    def get_contact_v_info(self):
        """
        TODO: get info from Hilary? from vamps?
        """
        pass
    def insert_test_contact(self):
        my_sql = '''INSERT IGNORE INTO contact (contact, email, institution, vamps_name, first_name, last_name)
                VALUES ("guest user", "guest@guest.com", "guest institution", "guest", "guest", "user")'''
        self.my_conn.execute_no_fetch(my_sql)        
        
    def get_contact_id(self, data_owner):
        my_sql = """SELECT contact_id FROM contact WHERE vamps_name = '%s'""" % (data_owner)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])        

    def insert_rundate(self):
        my_sql = """INSERT IGNORE INTO run (run, run_prefix) VALUES
            ('%s', 'illumina')""" % (self.rundate)
        self.run_id = self.my_conn.execute_no_fetch(my_sql)
        
    def insert_project(self, content_row, contact_id):
        my_sql = """INSERT IGNORE INTO project (project, title, project_description, rev_project_name, funding, env_sample_source_id, contact_id) VALUES
        ('%s', '%s', '%s', reverse('%s'), '%s', '%s', %s)
        """ % (content_row.project, content_row.project_title, content_row.project_description, content_row.project, content_row.funding, content_row.env_sample_source_id, contact_id)
#        print my_sql
        self.my_conn.execute_no_fetch(my_sql)

    def insert_dataset(self, content_row):
        """
        TODO: get dataset_description
        """        
        my_sql = """INSERT IGNORE INTO dataset (dataset, dataset_description) VALUES
        ('%s', '')
        """ % (content_row.dataset)
        self.my_conn.execute_no_fetch(my_sql)
    
    def insert_run_info(self, content_row):
        run_key_id      = self.get_id('run_key',      content_row.run_key)
        if not (self.run_id):
            self.run_id = self.get_id('run',          self.rundate)
        dataset_id      = self.get_id('dataset',      content_row.dataset)
        project_id      = self.get_id('project',      content_row.project)
        dna_region_id   = self.get_id('dna_region',   content_row.dna_region)
        primer_suite_id = self.get_id('primer_suite', content_row.primer_suite)
        file_prefix     = content_row.barcode_index + "_" + content_row.run_key + "_" + content_row.lane
        #overlap = content_row.overlap
        #if (content_row.overlap == 'complete'):
        #    overlap = 0
        
        my_sql = """INSERT IGNORE INTO run_info_ill (run_key_id, run_id, lane, dataset_id, project_id, tubelabel, barcode, 
                                                    adaptor, dna_region_id, amp_operator, seq_operator, barcode_index, overlap, insert_size, 
                                                    file_prefix, read_length, primer_suite_id) 
                                            VALUES (%s, %s, %s, %s, %s, '%s', '%s',  
                                                    '%s', %s, '%s', '%s', '%s', '%s', %s, 
                                                    '%s', %s, %s)
        """ % (run_key_id, self.run_id, content_row.lane, dataset_id, project_id, content_row.tubelabel, content_row.barcode, 
               content_row.adaptor, dna_region_id, content_row.amp_operator, content_row.seq_operator, content_row.barcode_index, content_row.overlap, content_row.insert_size,
                                                    file_prefix, content_row.read_length, primer_suite_id)
        
        self.my_conn.execute_no_fetch(my_sql)

    def insert_primer(self):
        pass
        
    def del_sequence_pdr_info_by_project_dataset(self, projects = "", datasets = "", primer_suite = ""):
        my_sql1 = """DELETE FROM sequence_pdr_info_ill
                    USING sequence_pdr_info_ill JOIN run_info_ill USING (run_info_ill_id) 
                    JOIN run USING(run_id) 
                    JOIN project using(project_id)
                    JOIN dataset using(dataset_id)
                    JOIN primer_suite using(primer_suite_id)
                    WHERE primer_suite = "%s"
                    AND run = "%s"
                """ % (primer_suite, self.rundate)
        my_sql2 = " AND project in (" + projects + ")"
        my_sql3 = " AND dataset in (" + datasets + ")"
        if (projects == "") and (datasets == ""):
            my_sql = my_sql1
        elif (projects != "") and (datasets == ""):
            my_sql = my_sql1 + my_sql2
        elif (projects == "") and (datasets != ""):
            my_sql = my_sql1 + my_sql3
        elif (projects != "") and (datasets != ""):
            my_sql = my_sql1 + my_sql2 + my_sql3
        self.my_conn.execute_no_fetch(my_sql)

#    def del_sequence_pdr_info(self):
#        my_sql = """DELETE FROM sequence_pdr_info_ill
#                    USING sequence_pdr_info_ill JOIN run_info_ill USING (run_info_ill_id) JOIN run USING(run_id) WHERE run = "%s"
#                """ % self.rundate
#        self.my_conn.execute_no_fetch(my_sql)
        
#    def del_run_info(self):
#        my_sql = """DELETE FROM run_info_ill
#                    USING run_info_ill JOIN run USING(run_id) WHERE run = "%s"
#                """ % self.rundate
#        self.my_conn.execute_no_fetch(my_sql)

    def del_run_info_by_project_dataset(self, projects = "", datasets = "", primer_suite = ""):
        my_sql1 = """DELETE FROM run_info_ill
                    USING run_info_ill 
                    JOIN run USING(run_id) 
                    JOIN project using(project_id)
                    JOIN primer_suite using(primer_suite_id)
                    WHERE primer_suite = "%s"
                    AND run = "%s"
                """ % (primer_suite, self.rundate)
        my_sql2 = " AND project in (" + projects + ")"
        my_sql3 = " AND dataset in (" + datasets + ")"
        if (projects == "") and (datasets == ""):
            my_sql = my_sql1
        elif (projects != "") and (datasets == ""):
            my_sql = my_sql1 + my_sql2
        elif (projects == "") and (datasets != ""):
            my_sql = my_sql1 + my_sql3
        elif (projects != "") and (datasets != ""):
            my_sql = my_sql1 + my_sql2 + my_sql3
        self.my_conn.execute_no_fetch(my_sql)


    def del_sequence_uniq_info(self):
        my_sql = """DELETE FROM sequence_uniq_info_ill 
                    USING sequence_uniq_info_ill 
                    LEFT JOIN sequence_pdr_info_ill USING(sequence_ill_id) 
                    WHERE sequence_pdr_info_ill_id is NULL"""
        self.my_conn.execute_no_fetch(my_sql)

    def del_sequences(self):
        my_sql = """DELETE FROM sequence_ill 
                    USING sequence_ill 
                    LEFT JOIN sequence_pdr_info_ill USING(sequence_ill_id) 
                    WHERE sequence_pdr_info_ill_id IS NULL
                """
        self.my_conn.execute_no_fetch(my_sql)


    def count_sequence_pdr_info_ill(self):
        primer_suite = self.get_primer_suite_name()
        lane         = self.get_lane().pop()   

        my_sql = """SELECT count(sequence_pdr_info_ill_id) 
                    FROM sequence_pdr_info_ill 
                      JOIN run_info_ill USING(run_info_ill_id) 
                      JOIN run USING(run_id) 
                      JOIN primer_suite using(primer_suite_id) 
                    WHERE run = '%s' 
                      AND lane = %s
                      AND primer_suite = '%s'
                      """ % (self.rundate, lane, primer_suite)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])   
    
    def get_primer_suite_name(self):
        primer_suite = [v.primer_suite for v in self.runobj.samples.itervalues()]
        return list(set(primer_suite))[0]
        
    def get_project_names(self):
        projects = [v.project for v in self.runobj.samples.itervalues()]
        return '", "'.join(set(projects))

    def get_dataset_names(self):
        datasets = [v.dataset for v in self.runobj.samples.itervalues()]
        return '", "'.join(set(datasets))

    def get_lane(self):
        lane = [v.lane for v in self.runobj.samples.itervalues()]
        return set(lane)

    def count_seq_from_file(self):
        try:
            with open(self.unique_file_counts) as fd:
                file_seq_orig = dict(line.strip().split(None, 1) for line in fd)
            file_seq_orig_count = sum([int(x) for x in file_seq_orig.values()])
            return file_seq_orig_count
        except IOError as e:
            print "Can't open file %s, error = %s" % (self.unique_file_counts, e)         
        except Exception:
            print "Unexpected error from 'count_seq_from_file':", sys.exc_info()[0]
            raise

    def check_seq_upload(self):
        file_seq_db_count   = self.count_sequence_pdr_info_ill()
#        print "file_seq_db_count = %s" % file_seq_db_count
        file_seq_orig_count = self.count_seq_from_file()
        if (file_seq_orig_count == file_seq_db_count):
            print "All sequences from files made it to the db: %s == %s" % (file_seq_orig_count, file_seq_db_count)
        else:
            print "Oops, amount of sequences from files not equal to the one in the db.\nIn file: %s != in db: %s\n==============" % (file_seq_orig_count, file_seq_db_count)
            
    def put_seq_statistics_in_file(self, filename, seq_in_file):
        pipelne_utils   = PipelneUtils()
#        if os.path.exists(file_full):
#            os.remove(file_full)
        pipelne_utils.write_seq_frequencies_in_file(self.unique_file_counts, filename, seq_in_file)       
        
