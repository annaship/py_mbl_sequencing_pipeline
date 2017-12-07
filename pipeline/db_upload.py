import sys
import os
import constants as C
from subprocess import Popen, PIPE
from shlex import split

from pipeline.get_ini import readCSV
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils
import IlluminaUtils.lib.fastalib as fastalib

try:
    import MySQLdb
except MySQLdb.Error, e:
    message = """
    MySQLdb ERROR
      To load the correct module, try running these commands before running the pipeline:
       
source /xraid/bioware/Modules/etc/profile.modules
module load bioware
    """
    PipelneUtils.print_both(message)
    PipelneUtils.print_both("Error %d: %s" % (e.args[0], e.args[1]))
    raise
except:                       # catch everything
    PipelneUtils.print_both("Unexpected:")
#     print "Unexpected:"         # handle unexpected exceptions
    PipelneUtils.print_both(sys.exc_info()[0])
#     print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
    raise          

#     sys.exit("""
#     MySQLdb ERROR
#       To load the correct module, try running these commands before running the pipeline:
#       
# source /xraid/bioware/Modules/etc/profile.modules
# module load bioware
# 
#     """)
class MyConnection:
    """
    Connection to env454
    Takes parameters from ~/.my.cnf, default host = "vampsdev", db="test"
    if different use my_conn = MyConnection(host, db)
    """
    def __init__(self, host="bpcweb7", db="test"):
# , read_default_file=os.path.expanduser("~/.my.cnf"), port = 3306
        
        self.utils  = PipelneUtils()        
        self.conn   = None
        self.cursor = None
        self.rows   = 0
        self.new_id = None
        self.lastrowid = None
        
        try:           
            self.utils.print_both("=" * 40)
            self.utils.print_both("host = " + str(host) + ", db = "  + str(db))
            self.utils.print_both("=" * 40)
            read_default_file = os.path.expanduser("~/.my.cnf")
            port_env = 3306
            
            if self.utils.is_local():
                host = "127.0.0.1"
#                 if db == "env454":
#                     port_env = 3308
#                     read_default_file = os.path.expanduser("~/.my.cnf_server")
#                 else:
#                     db = "test_env454"
                read_default_file = "~/.my.cnf_local"
            self.conn   = MySQLdb.connect(host = host, db = db, read_default_file = read_default_file, port = port_env)
            self.cursor = self.conn.cursor()
            # self.escape = self.conn.escape()
                   
        except MySQLdb.Error, e:
            self.utils.print_both("Error %d: %s" % (e.args[0], e.args[1]))
            raise
        except:                       # catch everything
            self.utils.print_both("Unexpected:")
            self.utils.print_both(sys.exc_info()[0])
#             print "Unexpected:"         # handle unexpected exceptions
#             print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
            raise                       # re-throw caught exception   

    def execute_fetch_select(self, sql):
        if self.cursor:
          try:
            # sql = self.conn.escape(sql)
            self.cursor.execute(sql)
            res = self.cursor.fetchall ()
          except:
            self.utils.print_both(("ERROR: query = %s") % sql)
            raise
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
        self.utils       = PipelneUtils()
        self.runobj      = runobj
        self.rundate     = self.runobj.run
        self.use_cluster = 1       
        self.unique_fasta_files = []
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
        
        self.dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, lane_name = lane_name, site = site) 
 
        
        self.analysis_dir = self.dirs.check_dir(self.dirs.analysis_dir)
        self.fasta_dir    = self.dirs.check_dir(self.dirs.reads_overlap_dir)
        self.gast_dir     = self.dirs.check_dir(self.dirs.gast_dir)

        host_name     = runobj.database_host
        database_name = runobj.database_name
        
        self.filenames   = []
        # logger.error("self.utils.is_local() LLL1 db upload")
        # logger.error(self.utils.is_local())
        
        if self.utils.is_local():
            self.my_conn = MyConnection(host = 'localhost', db="test_env454")
        else:
            self.my_conn = MyConnection(host='bpcdb1', db="env454")
#             self.my_conn = MyConnection(host='bpcdb1.jbpc-np.mbl.edu', db="env454")
        self.sequence_table_name = "sequence_ill" 
        self.sequence_field_name = "sequence_comp" 
        self.my_csv              = None

        self.unique_file_counts = self.dirs.unique_file_counts
        self.dirs.delete_file(self.unique_file_counts)
        self.seq_id_dict = {}
        self.tax_id_dict = {}
        self.run_id      = None
#        self.nonchimeras_suffix = ".nonchimeric.fa"
        self.nonchimeric_suffix = "." + C.nonchimeric_suffix #".nonchimeric.fa"
        self.fa_unique_suffix   = ".fa." + C.unique_suffix #.fa.unique
        self.v6_unique_suffix   = "MERGED_V6_PRIMERS_REMOVED." + C.unique_suffix
        self.suff_list = [self.nonchimeric_suffix, self.fa_unique_suffix, self.v6_unique_suffix]

#         self.merge_unique_suffix = "." + C.filtered_suffix + "." + C.unique_suffix #.MERGED-MAX-MISMATCH-3.unique
        self.suffix_used        = ""
        
#        self.refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
   
   
    def get_fasta_file_names(self):
        files_names = self.dirs.get_all_files(self.fasta_dir)
        self.unique_fasta_files = [f for f in files_names.keys() if f.endswith(tuple(self.suff_list))]
# needs return because how it's called from pipelineprocesor
        return self.unique_fasta_files
        

    def get_run_info_ill_id(self, filename_base):
        
        my_sql = """SELECT run_info_ill_id FROM run_info_ill 
                    JOIN run using(run_id)
                    WHERE file_prefix = '%s'
                    and run = '%s'
        """ % (filename_base, self.rundate)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])
        
    def make_seq_upper(self, filename):
        read_fasta = fastalib.ReadFasta(filename)
        sequences  = [seq.upper() for seq in read_fasta.sequences] #here we make uppercase for VAMPS compartibility    
        read_fasta.close()
        return sequences 
        
    def insert_seq(self, sequences):
      query_tmpl = "INSERT INTO %s (%s) VALUES (COMPRESS(%s))"
      val_tmpl   = "'%s'"
      my_sql     = query_tmpl % (self.sequence_table_name, self.sequence_field_name, ')), (COMPRESS('.join([val_tmpl % key for key in sequences]))
      my_sql     = my_sql + " ON DUPLICATE KEY UPDATE %s = VALUES(%s)" % (self.sequence_field_name, self.sequence_field_name)
#       print "MMM my_sql = %s" % my_sql
      seq_id     = self.my_conn.execute_no_fetch(my_sql)
      self.utils.print_both("sequences in file: %s\n" % (len(sequences)))
      return seq_id
    #     try:
    #         query_tmpl = "INSERT IGNORE INTO %s (%s) VALUES (COMPRESS(%s))"
    #         val_tmpl   = "'%s'"
    #         my_sql     = query_tmpl % (self.sequence_table_name, self.sequence_field_name, ')), (COMPRESS('.join([val_tmpl % key for key in sequences]))
    #         seq_id     = self.my_conn.execute_no_fetch(my_sql)
    # #         print "sequences in file: %s" % (len(sequences))
    #         self.utils.print_both("sequences in file: %s\n" % (len(sequences)))
    #         return seq_id
    #     except self.my_conn.conn.cursor._mysql_exceptions.Error as err:
    #         if err.errno == 1582:
    #             self.utils.print_both(("ERROR: _mysql_exceptions.OperationalError: (1582, \"Incorrect parameter count in the call to native function 'COMPRESS'\"), there is an empty fasta in %s") % self.fasta_dir)
    #         else:
    #             raise
    #     except:
    #         if len(sequences) == 0:
    #             self.utils.print_both(("ERROR: There are no sequences, please check if there are correct fasta files in the directory %s") % self.fasta_dir)
    #         raise
        
    def get_seq_id_dict(self, sequences):
        id_name    = self.sequence_table_name + "_id" 
        query_tmpl = """SELECT %s, uncompress(%s) FROM %s WHERE %s in (COMPRESS(%s))"""
        val_tmpl   = "'%s'"
        try:
            my_sql     = query_tmpl % (id_name, self.sequence_field_name, self.sequence_table_name, self.sequence_field_name, '), COMPRESS('.join([val_tmpl % key for key in sequences]))
            res        = self.my_conn.execute_fetch_select(my_sql)
            one_seq_id_dict = dict((y, int(x)) for x, y in res)
            self.seq_id_dict.update(one_seq_id_dict)
        except:
            if len(sequences) == 0:
                self.utils.print_both(("ERROR: There are no sequences, please check if there are correct fasta files in the directory %s") % self.fasta_dir)
            raise


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
        res_id = ""
        if (not run_info_ill_id):
            self.utils.print_both("ERROR: There is no run info yet, please check if it's uploaded to env454")
            
        # ------- insert sequence info per run/project/dataset --------
        seq_upper = fasta.seq.upper()
        sequence_ill_id = self.seq_id_dict[seq_upper]

        seq_count       = int(fasta.id.split('|')[-1].split(':')[-1])
#        print run_info_ill_id, sequence_ill_id, seq_count
        my_sql          = "INSERT INTO sequence_pdr_info_ill (run_info_ill_id, sequence_ill_id, seq_count) VALUES (%s, %s, %s)" % (run_info_ill_id, sequence_ill_id, seq_count)
        my_sql          = my_sql + " ON DUPLICATE KEY UPDATE run_info_ill_id = VALUES(run_info_ill_id), sequence_ill_id = VALUES(sequence_ill_id), seq_count = VALUES(seq_count);"
#         print "MMM1 my_sql = %s" % my_sql
#         try:
#             res_id = self.my_conn.execute_no_fetch(my_sql)
#             return res_id
#         except:
#             self.utils.print_both("Offensive query: %s" % my_sql)
#             raise
        return my_sql
        
    def make_gast_files_dict(self):
        return self.dirs.get_all_files(self.gast_dir, "gast")
        
        
    def gast_filename(self, filename):
#         todo: if filename in make_gast_files_dict, use it full path
        gast_file_names = self.make_gast_files_dict()
        gast_file_name_path = ""
        for gast_file_name_path, tpls in gast_file_names.iteritems():
            if any(t.endswith(filename) for t in tpls):
                return gast_file_name_path 
    
    def get_gast_result(self, filename):
        gast_file_name = self.gast_filename(filename)
        self.utils.print_both("current gast_file_name = %s." % gast_file_name)
        
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
#         except OSError, e:
        except TypeError, e:
            self.utils.print_both("Check if there is a gast file under %s for %s." % (self.gast_dir, filename))
            pass            
        except:
            # reraise the exception, as it's an unexpected error
            raise

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
#                     tax_id = self.my_conn.execute_no_fetch(my_sql)
#                     self.tax_id_dict[taxonomy] = tax_id
#                 return tax_id
                    return my_sql

        else:
            self.utils.print_both("ERROR: can't read gast files! No taxonomy information will be processed. Please check if gast results are in analysis/gast")
#             logger.debug("ERROR: can't read gast files! No taxonomy information will be processed.")            

    def insert_sequence_uniq_info_ill(self, fasta, gast_dict):
#     def make_insert_sequence_uniq_info_ill(self, fasta, gast_dict):
        cnt = 0
        if gast_dict:
            (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast_dict[fasta.id]
            seq_upper = fasta.seq.upper()
            sequence_ill_id = self.seq_id_dict[seq_upper]
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
                   ON DUPLICATE KEY UPDATE
                       updated = (CASE WHEN taxonomy_id <> %s THEN NOW() ELSE updated END),
                       taxonomy_id = %s,
                       gast_distance = '%s',
                       refssu_count = '%s',
                       rank_id = (SELECT rank_id FROM rank WHERE rank = '%s'),
                       refhvr_ids = '%s';
                   """ % (sequence_ill_id, taxonomy_id, distance, refssu_count, rank, refhvr_ids.rstrip(), taxonomy_id, taxonomy_id, distance, refssu_count, rank, refhvr_ids.rstrip())
                     
#             my_sql = """INSERT INTO sequence_uniq_info_ill (sequence_ill_id, taxonomy_id, gast_distance, refssu_count, rank_id, refhvr_ids) VALUES
#                    (
#                     %s,
#                     %s,
#                     '%s',
#                     '%s',
#                     (SELECT rank_id FROM rank WHERE rank = '%s'),
#                     '%s'                
#                    )
#                    """ % (sequence_ill_id, taxonomy_id, distance, refssu_count, rank, refhvr_ids.rstrip())
#             my_sql          = my_sql + " ON DUPLICATE KEY UPDATE sequence_ill_id = VALUES(sequence_ill_id)"
#             my_sql          = my_sql + """ ON DUPLICATE KEY UPDATE sequence_ill_id = VALUES(sequence_ill_id), 
#             taxonomy_id = VALUES(taxonomy_id), 
#             gast_distance = VALUES(gast_distance), 
#             refssu_count = VALUES(refssu_count), 
#             rank_id = VALUES(rank_id), 
#             refhvr_ids = VALUES(refhvr_ids),
#             updated = (CASE WHEN VALUES(taxonomy_id) <> %s THEN NOW() ELSE updated END)
#             """
            res_id = self.my_conn.execute_no_fetch(my_sql)
            return res_id

    # def update_sequence_uniq_info_ill(self, fasta, gast_dict):
    #     if gast_dict:
    #         (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast_dict[fasta.id]
    #         seq_upper = fasta.seq.upper()
    #         sequence_ill_id = self.seq_id_dict[seq_upper]
    #         if taxonomy in self.tax_id_dict:
    #             try:
    #                 taxonomy_id = self.tax_id_dict[taxonomy] 
    #             except Exception, e:
    #                 logger.debug("Error = %s" % e)
    #                 raise
    #                   
    #         my_sql = """INSERT IGNORE INTO sequence_uniq_info_ill (sequence_ill_id, taxonomy_id, gast_distance, refssu_count, rank_id, refhvr_ids) VALUES
    #                (
    #                 %s,
    #                 %s,
    #                 '%s',
    #                 '%s',
    #                 (SELECT rank_id FROM rank WHERE rank = '%s'),
    #                 '%s'                
    #                )
    #                ON DUPLICATE KEY UPDATE
    #                    updated = (CASE WHEN taxonomy_id <> %s THEN NOW() ELSE updated END),
    #                    taxonomy_id = %s,
    #                    gast_distance = '%s',
    #                    refssu_count = '%s',
    #                    rank_id = (SELECT rank_id FROM rank WHERE rank = '%s'),
    #                    refhvr_ids = '%s'
    #                """ % (sequence_ill_id, taxonomy_id, distance, refssu_count, rank, refhvr_ids.rstrip(), taxonomy_id, taxonomy_id, distance, refssu_count, rank, refhvr_ids.rstrip())
    #                 
    #         res_id = self.my_conn.execute_no_fetch(my_sql)
    #         return res_id
    

    def put_run_info(self, content = None):

        run_keys = list(set([run_key.split('_')[1] for run_key in self.runobj.run_keys]))
        self.insert_bulk_data('run_key', run_keys)
        dna_regions = list(set([self.runobj.samples[key].dna_region for key in self.runobj.samples]))
        self.insert_bulk_data('dna_region', dna_regions)
        self.insert_rundate()

        for key in self.runobj.samples:
            value = self.runobj.samples[key]
            self.get_contact_v_info()
            contact_id = self.get_contact_id(value.data_owner)
            self.insert_project(value, contact_id)
            self.insert_dataset(value) 

            self.insert_run_info(value)

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
        my_sql = """INSERT IGNORE INTO run (run, run_prefix, platform) VALUES
            ('%s', 'illumin', '%s')""" % (self.rundate, self.runobj.platform)
        self.run_id = self.my_conn.execute_no_fetch(my_sql)
        
    def insert_project(self, content_row, contact_id):
        if (not contact_id):
            self.utils.print_both("ERROR: There is no such contact info on env454, please check if the user has an account on VAMPS")        
        my_sql = """INSERT IGNORE INTO project (project, title, project_description, rev_project_name, funding, env_sample_source_id, contact_id) VALUES
        ('%s', '%s', '%s', reverse('%s'), '%s', '%s', %s)
        """ % (content_row.project, content_row.project_title, content_row.project_description, content_row.project, content_row.funding, content_row.env_sample_source_id, contact_id)
        self.utils.print_both(my_sql)
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
        illumina_index_id = self.get_id('illumina_index', content_row.barcode_index)
        file_prefix     = content_row.barcode_index + "_" + content_row.run_key + "_" + content_row.lane
        #overlap = content_row.overlap
        #if (content_row.overlap == 'complete'):
        #    overlap = 0
        
        my_sql = """INSERT IGNORE INTO run_info_ill (run_key_id, run_id, lane, dataset_id, project_id, tubelabel, barcode, 
                                                    adaptor, dna_region_id, amp_operator, seq_operator, overlap, insert_size, 
                                                    file_prefix, read_length, primer_suite_id, platform, illumina_index_id) 
                                            VALUES (%s, %s, %s, %s, %s, '%s', '%s',  
                                                    '%s', %s, '%s', '%s', '%s', %s, 
                                                    '%s', %s, %s, '%s', %s)
        """ % (run_key_id, self.run_id, content_row.lane, dataset_id, project_id, content_row.tubelabel, content_row.barcode, 
               content_row.adaptor, dna_region_id, content_row.amp_operator, content_row.seq_operator, content_row.overlap, content_row.insert_size,
                                                    file_prefix, content_row.read_length, primer_suite_id, self.runobj.platform, illumina_index_id)
        
        self.utils.print_both("insert run_info sql = %s" % my_sql)
        
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
        results = {}
        primer_suites = self.get_primer_suite_name()
        lane          = self.get_lane().pop()
        for primer_suite in primer_suites:
            primer_suite_lane = primer_suite + ", lane " + lane
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
            try:
                if (int(res[0][0]) > 0):
                    results[primer_suite_lane] = int(res[0][0])
#                     results.append(int(res[0][0]))
            except Exception:
                self.utils.print_both("Unexpected error from 'count_sequence_pdr_info_ill':", sys.exc_info()[0])
                raise                
        return results
#             int(res[0][0])   
    
    def get_primer_suite_name(self):
        primer_suites = [v.primer_suite for v in self.runobj.samples.itervalues()]
        return list(set(primer_suites))
        
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
            self.utils.print_both("Can't open file %s, error = %s" % (self.unique_file_counts, e))         
        except Exception:
            self.utils.print_both("Unexpected error from 'count_seq_from_file':", sys.exc_info()[0])
            raise
        
    def count_seq_from_files_grep(self):
#         grep '>' *-PERFECT_reads.fa.unique
#       or
#         cd /xraid2-2/g454/run_new_pipeline/illumina/20130607/lane_5_A/analysis/reads_overlap/; grep '>' *_MERGED-MAX-MISMATCH-3.unique.nonchimeric.fa | wc -l; date
        try:
            self.suffix_used = list(set([ext for f in self.unique_fasta_files for ext in self.suff_list if f.endswith(ext)]))[0] 
        except:
            print "self.unique_fasta_files = %s, self.suff_list = %s" % (self.unique_fasta_files, self.suff_list)
            self.suffix_used = ""
        print self.suffix_used
        suffix = self.fasta_dir + "/*" + self.suffix_used 
        program_name = "grep"
        call_params  = " '>' " + suffix
        command_line = program_name + call_params
        p1 = Popen(command_line, stdout=PIPE, shell=True)
        p2 = Popen(split("wc -l"), stdin=p1.stdout, stdout=PIPE)
#         output = p2.stdout.read().split(" ")[0].strip()
        output, err = p2.communicate()
#         print output
        return int(output.strip())           


    def check_seq_upload(self):
        file_seq_db_counts   = self.count_sequence_pdr_info_ill()
#        print "file_seq_db_count = %s" % file_seq_db_count
#         file_seq_orig_count = self.count_seq_from_file()
        file_seq_orig_count = self.count_seq_from_files_grep()
        
        for pr_suite, file_seq_db_count in file_seq_db_counts.iteritems():
            if (file_seq_orig_count == file_seq_db_count):
                self.utils.print_both("All sequences from files made it to the db for %s: %s == %s\n" % (pr_suite, file_seq_orig_count, file_seq_db_count))
            else:
                self.utils.print_both("Warning: Amount of sequences from files not equal to the one in the db for %s: %s != %s\n" % (pr_suite, file_seq_orig_count, file_seq_db_count))
#                 
#                 ("Oops, amount of sequences from files not equal to the one in the db for %.\nIn file: %s != in db: %s\n==============" % (pr_suite, file_seq_orig_count, file_seq_db_count))
            
    def put_seq_statistics_in_file(self, filename, seq_in_file):
#        if os.path.exists(file_full):
#            os.remove(file_full)
        self.utils.write_seq_frequencies_in_file(self.unique_file_counts, filename, seq_in_file)       
        
