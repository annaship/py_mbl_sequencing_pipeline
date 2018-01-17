import sys
import os
import constants as C
from subprocess import Popen, PIPE
from shlex import split
from time import sleep, time, gmtime, strftime

from pipeline.get_ini import readCSV
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils
import IlluminaUtils.lib.fastalib as fastalib
from collections import defaultdict
from itertools import izip_longest

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
#     PipelneUtils.print_both("Unexpected:")
    print "Unexpected:"         # handle unexpected exceptions
#     PipelneUtils.print_both(sys.exc_info()[0])
    print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
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
        self.cursorD = None
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
            self.cursorD = self.conn.cursor (MySQLdb.cursors.DictCursor)
            # self.escape = self.conn.escape()

        except (AttributeError, MySQLdb.OperationalError):
            self.conn = MySQLdb.connect(host=host, db=db, read_default_file=read_default_file, port=port_env)
            self.cursor = self.conn.cursor()


        except MySQLdb.Error, e:
            self.utils.print_both("Error %d: %s" % (e.args[0], e.args[1]))
            raise
        except:                       # catch everything
            self.utils.print_both("Unexpected:")
            self.utils.print_both(sys.exc_info()[0])
#             print "Unexpected:"         # handle unexpected exceptions
#             print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
            raise                       # re-throw caught exception


    def connect(self, host, db, read_default_file, port_env):
        return MySQLdb.connect(host = host, db = db, read_default_file = read_default_file, port = port_env)


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
#             if (self.conn.affected_rows()):
#                 logger.debug("affected_rows = "  + str(self.conn.affected_rows()))
#             logger.debug(self.cursor._info)
#             _info    str: Records: 238  Duplicates: 66  Warnings: 0
            return self.cursor._info


    def execute_insert(self, table_name, field_name, val_list, ignore = "IGNORE"):
        try:
            sql = "INSERT %s INTO %s (%s) VALUES (%s) " % (ignore, table_name, field_name, val_list)
            sql = sql + " ON DUPLICATE KEY UPDATE %s = VALUES(%s);" % (field_name, field_name)

#             print 'sql',sql
            #if table_name == 'dataset' or table_name == 'project':
            #    print 'sql',sql
            if self.cursor:
                self.cursor.execute(sql)
                self.conn.commit()
                return (self.cursor.rowcount, self.cursor.lastrowid)
        except:
            self.utils.print_both(("ERROR: query = %s") % sql)
            raise

    def get_all_name_id(self, table_name, id_name = "", field_name = "", where_part = ""):
        if (field_name == ""):
            field_name = table_name
        if (id_name == ""):
            id_name = table_name + '_id'
        my_sql  = """SELECT %s, %s FROM %s %s""" % (field_name, id_name, table_name, where_part)
#         self.utils.print_both(("my_sql from get_all_name_id = %s") % my_sql)
        res     = self.execute_fetch_select(my_sql)

        if res:
            return res

    def run_groups(self, group_vals, query_tmpl, join_xpr = ', '):
        for group in group_vals:
            val_part = join_xpr.join([key for key in group if key is not None])
            my_sql = query_tmpl % (val_part)
            insert_info = self.execute_no_fetch(my_sql)
            logger.debug("insert info = %s" % insert_info)

    def make_sql_for_groups(self, table_name, fields):
        field_list = fields.split(",")
        my_sql_1 = "INSERT IGNORE INTO %s (%s) VALUES " % (table_name, fields)
        my_sql_2 =  " ON DUPLICATE KEY UPDATE "
        for field_name in field_list[:-1]:
            my_sql_2 = my_sql_2 + " %s = VALUES(%s), " % (field_name.strip(), field_name.strip())
        my_sql_2 = my_sql_2 + "  %s = VALUES(%s);" % (field_list[-1].strip(), field_list[-1].strip())
        return my_sql_1 + " %s " + my_sql_2


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
    def __init__(self, runobj = None, db_server = None):
        if db_server is None:
            db_server = "env454"

        self.db_server   = db_server
        self.utils       = PipelneUtils()
        self.runobj      = runobj
        self.rundate     = self.runobj.run
        self.use_cluster = 1
        self.unique_fasta_files = []

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

        self.filenames   = []
        # logger.error("self.utils.is_local() LLL1 db upload")
        # logger.error(self.utils.is_local())
        self.table_names_dict = {"vamps2": {"sequence_field_name": "sequence_comp", "sequence_table_name": "sequence",
                                             "sequence_pdr_info_table_name": "sequence_pdr_info", "contact": "user", "username": "username"},
                                 "env454": {"sequence_field_name": "sequence_comp", "sequence_table_name": "sequence_ill",
                                             "sequence_pdr_info_table_name": "sequence_pdr_info_ill", "contact": "contact", "username": "vamps_name"}}

        self.db_cnf = {
            "vamps2": {"local":      {"host": "localhost", "db": "vamps2"},
                       "production": {"host": "vampsdev",   "db": "vamps2"}
                       },
            "env454": {"local":      {"host": "localhost", "db": "test_env454"},
                       "production": {"host": "bpcdb1",    "db": "env454"}
                       }
                    }

        if self.utils.is_local():
            is_local = "local"
        else:
            is_local = "production"

        self.table_names = self.table_names_dict[self.db_server]

        host = self.db_cnf[self.db_server][is_local]["host"]
        db   = self.db_cnf[self.db_server][is_local]["db"]
        self.my_conn = MyConnection(host, db)

        self.taxonomy = Taxonomy(self.my_conn)
        self.seq      = Seq(self.taxonomy, self.table_names)

        self.gast_dict = {}
        self.silva_taxonomy_info_per_seq_list = []

        self.unique_file_counts = self.dirs.unique_file_counts
        self.dirs.delete_file(self.unique_file_counts)
        self.taxonomies = set()
        self.run_id      = None
        self.nonchimeric_suffix = "." + C.nonchimeric_suffix #".nonchimeric.fa"
        self.fa_unique_suffix   = ".fa." + C.unique_suffix #.fa.unique
        self.v6_unique_suffix   = "MERGED_V6_PRIMERS_REMOVED." + C.unique_suffix
        self.suff_list = [self.nonchimeric_suffix, self.fa_unique_suffix, self.v6_unique_suffix]
        self.suffix_used        = ""
        self.all_dataset_ids = self.my_conn.get_all_name_id("dataset")
        if db_server == "vamps2":
            self.put_run_info()
        self.all_dataset_run_info_dict = self.get_dataset_per_run_info_id()


    def get_fasta_file_names(self):
        files_names = self.dirs.get_all_files(self.fasta_dir)
        self.unique_fasta_files = [f for f in files_names.keys() if f.endswith(tuple(self.suff_list))]
# needs return because how it's called from pipelineprocesor
        return self.unique_fasta_files

    def send_mail(self):
        program_name1 = "echo"
        call_params1 = " 'Projects ? were uploaded to VAMPS2'"
        command_line1 = program_name1 + call_params1

        program_name2 = "mail"
        call_params2 = ' -s "Vamps2 upload" ashipunova@mbl.edu'
        command_line2 = program_name2 + call_params2

        p1 = Popen(command_line1, stdout=PIPE, shell=True)
        p2 = Popen(command_line2, stdin=p1.stdout, stdout=PIPE)
        #         output = p2.stdout.read().split(" ")[0].strip()
        output, err = p2.communicate()
        print output

    def get_run_info_ill_id(self, filename_base):

        my_sql = """SELECT run_info_ill_id FROM run_info_ill
                    JOIN run using(run_id)
                    WHERE file_prefix = '%s'
                    and run = '%s';
        """ % (filename_base, self.rundate)
        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])

    def get_dataset_per_run_info_id(self):
        all_dataset_run_info_sql = "SELECT run_info_ill_id, dataset_id FROM run_info_ill"
        res = self.my_conn.execute_fetch_select(all_dataset_run_info_sql)
        return dict([(r, d) for r, d in res])

    def get_id(self, table_name, value):
        id_name = table_name + '_id'
        my_sql  = """SELECT %s FROM %s WHERE %s = '%s';""" % (id_name, table_name, table_name, value)
        res     = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])


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
                gast_content = fd.readlines()
            self.gast_dict = dict([(l.split("\t")[0], l.split("\t")[1:]) for l in gast_content[1:]])
#             gast_dict.remove([k for k in gast_dict if k[0] == 'taxonomy'][0])
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

    def put_run_info(self, content = None):

        run_keys = list(set([run_key.split('_')[1] for run_key in self.runobj.run_keys]))
        self.insert_bulk_data('run_key', run_keys)
        dna_regions = list(set([self.runobj.samples[key].dna_region for key in self.runobj.samples]))
        self.insert_bulk_data('dna_region', dna_regions)
        self.insert_rundate()

        for key in self.runobj.samples:
            value = self.runobj.samples[key]
#             self.get_contact_v_info()
            contact_id = self.get_contact_id(value.data_owner)
            if (not contact_id):
                logger.error("""ERROR: There is no such contact info on %s,
                    please check if the user %s has an account on VAMPS""" % (self.db_server, value.data_owner))
                sys.exit("""ERROR: There is no such contact info on %s,
                    please check if the user %s has an account on VAMPS""" % (self.db_server, value.data_owner))
            self.insert_project(value, contact_id)
            self.insert_dataset(value)

            self.insert_run_info(value)

    def insert_bulk_data(self, key, values):
        query_tmpl = "INSERT IGNORE INTO %s (%s) VALUES (%s)"
        val_tmpl   = "'%s'"
        my_sql     = query_tmpl % (key, key, '), ('.join([val_tmpl % v for v in values]))
        my_sql     = my_sql + " ON DUPLICATE KEY UPDATE %s = VALUES(%s);" % (key, key)

        cursor_info = self.my_conn.execute_no_fetch(my_sql)

    def get_contact_v_info(self):
        """
        TODO: get info from Hilary? from vamps?
        """
        pass
    def insert_test_contact(self):
        my_sql = '''INSERT IGNORE INTO contact (contact, email, institution, vamps_name, first_name, last_name)
                VALUES ("guest user", "guest@guest.com", "guest institution", "guest", "guest", "user");'''
        return self.my_conn.execute_no_fetch(my_sql)

    def get_contact_id(self, data_owner):
        my_sql = """SELECT %s_id FROM %s WHERE %s = '%s';""" % (self.table_names["contact"], self.table_names["contact"], self.table_names["username"], data_owner)

        res    = self.my_conn.execute_fetch_select(my_sql)
        if res:
            return int(res[0][0])

    def insert_rundate(self):
        my_sql = """INSERT IGNORE INTO run (run, run_prefix, platform) VALUES
            ('%s', 'illumin', '%s');""" % (self.rundate, self.runobj.platform)
        return self.my_conn.execute_no_fetch(my_sql)

    def insert_project(self, content_row, contact_id):
        if (not contact_id):
            self.utils.print_both("ERROR: There is no such contact info on env454, please check if the user has an account on VAMPS")

        if self.db_server == "vamps2":
            fields = "project, title, project_description, rev_project_name, funding, owner_user_id, created_at"
            vals = """('%s', '%s', '%s', reverse('%s'), '%s', '%s', NOW())
            """ % (content_row.project, content_row.project_title, content_row.project_description, content_row.project, content_row.funding, contact_id)
            group_vals = self.utils.grouper([vals], 1)
            query_tmpl = self.my_conn.make_sql_for_groups("project", fields)
            self.my_conn.run_groups(group_vals, query_tmpl)

        elif self.db_server == "env454":
            my_sql = """INSERT IGNORE INTO project (project, title, project_description, rev_project_name, funding, env_sample_source_id, contact_id) VALUES
                ('%s', '%s', '%s', reverse('%s'), '%s', '%s', %s);
                """ % (content_row.project, content_row.project_title, content_row.project_description, content_row.project, content_row.funding, content_row.env_sample_source_id, contact_id)
#         TODO: change! what if we have more self.db_server?
            self.utils.print_both(my_sql)
            cursor_info = self.my_conn.execute_no_fetch(my_sql)

    def insert_dataset(self, content_row):
        if self.db_server == "vamps2":
            project_id = self.get_id('project', content_row.project)
            my_sql = """INSERT IGNORE INTO dataset (dataset, dataset_description, project_id, created_at) VALUES
                ('%s', '%s', %s, NOW());
                """ % (content_row.dataset, content_row.dataset_description, project_id)
        elif self.db_server == "env454":
            my_sql = """INSERT IGNORE INTO dataset (dataset, dataset_description) VALUES
                ('%s', '%s');
                """ % (content_row.dataset, content_row.dataset_description)
        return self.my_conn.execute_no_fetch(my_sql)

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
        
        # no project_id
        if (self.db_server == "vamps2"):
            my_sql = """INSERT IGNORE INTO run_info_ill (run_key_id, run_id, lane, dataset_id, tubelabel, barcode,
                                                    adaptor, dna_region_id, amp_operator, seq_operator, overlap, insert_size,
                                                    file_prefix, read_length, primer_suite_id, platform, illumina_index_id)
                                            VALUES (%s, %s, %s, %s, '%s', '%s',
                                                    '%s', %s, '%s', '%s', '%s', %s,
                                                    '%s', %s, %s, '%s', %s);
        """ % (run_key_id, self.run_id, content_row.lane, dataset_id, content_row.tubelabel, content_row.barcode,
               content_row.adaptor, dna_region_id, content_row.amp_operator, content_row.seq_operator, content_row.overlap, content_row.insert_size,
                                                    file_prefix, content_row.read_length, primer_suite_id, self.runobj.platform, illumina_index_id)

        elif (self.db_server == "env454"):
            my_sql = """INSERT IGNORE INTO run_info_ill (run_key_id, run_id, lane, dataset_id, project_id, tubelabel, barcode,
                                                    adaptor, dna_region_id, amp_operator, seq_operator, overlap, insert_size,
                                                    file_prefix, read_length, primer_suite_id, platform, illumina_index_id)
                                            VALUES (%s, %s, %s, %s, %s, '%s', '%s',
                                                    '%s', %s, '%s', '%s', '%s', %s,
                                                    '%s', %s, %s, '%s', %s);
        """ % (run_key_id, self.run_id, content_row.lane, dataset_id, project_id, content_row.tubelabel, content_row.barcode,
               content_row.adaptor, dna_region_id, content_row.amp_operator, content_row.seq_operator, content_row.overlap, content_row.insert_size,
                                                    file_prefix, content_row.read_length, primer_suite_id, self.runobj.platform, illumina_index_id)


#         my_sql = """INSERT IGNORE INTO run_info_ill (run_key_id, run_id, lane, dataset_id, project_id, tubelabel, barcode,
#                                                     adaptor, dna_region_id, amp_operator, seq_operator, overlap, insert_size,
#                                                     file_prefix, read_length, primer_suite_id, platform, illumina_index_id)
#                                             VALUES (%s, %s, %s, %s, %s, '%s', '%s',
#                                                     '%s', %s, '%s', '%s', '%s', %s,
#                                                     '%s', %s, %s, '%s', %s);
#         """ % (run_key_id, self.run_id, content_row.lane, dataset_id, project_id, content_row.tubelabel, content_row.barcode,
#                content_row.adaptor, dna_region_id, content_row.amp_operator, content_row.seq_operator, content_row.overlap, content_row.insert_size,
#                                                     file_prefix, content_row.read_length, primer_suite_id, self.runobj.platform, illumina_index_id)

#         self.utils.print_both("insert run_info sql = %s" % my_sql)

        cursor_info = self.my_conn.execute_no_fetch(my_sql)
        self.utils.print_both("insert run_info: %s" % cursor_info)

    def insert_primer(self):
        pass

    def del_sequence_pdr_info_by_project_dataset(self, projects = "", datasets = "", primer_suite = ""):
        my_sql1 = """DELETE FROM %s
                    USING %s JOIN run_info_ill USING (run_info_ill_id)
                    JOIN run USING(run_id)
                    JOIN project using(project_id)
                    JOIN dataset using(dataset_id)
                    JOIN primer_suite using(primer_suite_id)
                    WHERE primer_suite = "%s"
                    AND run = "%s"
                """ % (self.table_names["sequence_pdr_info_table_name"], self.table_names["sequence_pdr_info_table_name"],  primer_suite, self.rundate)
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
        cursor_info = self.my_conn.execute_no_fetch(my_sql)

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
        cursor_info = self.my_conn.execute_no_fetch(my_sql)

    def del_sequence_uniq_info(self):
        my_sql = """DELETE FROM sequence_uniq_info_ill
                    USING sequence_uniq_info_ill
                    LEFT JOIN %s USING(%s_id)
                    WHERE %s_id is NULL;""" % (self.table_names["sequence_pdr_info_table_name"], self.table_names["sequence_table_name"], self.table_names["sequence_pdr_info_table_name"])
        cursor_info = self.my_conn.execute_no_fetch(my_sql)

    def del_sequences(self):
        my_sql = """DELETE FROM %s
                    USING %s
                    LEFT JOIN %s USING(%s_id)
                    WHERE %s_id IS NULL;
                """ % (self.table_names["sequence_table_name"], self.table_names["sequence_table_name"], self.table_names["sequence_table_name"], self.table_names["sequence_pdr_info_table_name"], self.table_names["sequence_pdr_info_table_name"])
        cursor_info = self.my_conn.execute_no_fetch(my_sql)

    def count_sequence_pdr_info(self):
        results = {}
        primer_suites = self.get_primer_suite_name()
        lane          = self.get_lane().pop()

        if (self.db_server == "vamps2"):
            join_add = """ JOIN dataset using(dataset_id)
                       JOIN run_info_ill USING(dataset_id) """
        elif (self.db_server == "env454"):
            join_add = """ JOIN run_info_ill USING(run_info_ill_id) """

        for primer_suite in primer_suites:
            primer_suite_lane = primer_suite + ", lane " + lane
            my_sql = """SELECT count(%s_id)
                        FROM %s
                          %s
                          JOIN run USING(run_id)
                          JOIN primer_suite using(primer_suite_id)
                        WHERE run = '%s'
                          AND lane = %s
                          AND primer_suite = '%s';
                          """ % (self.table_names["sequence_pdr_info_table_name"], self.table_names["sequence_pdr_info_table_name"], join_add, self.rundate, lane, primer_suite)
            res    = self.my_conn.execute_fetch_select(my_sql)
            try:
                if (int(res[0][0]) > 0):
                    results[primer_suite_lane] = int(res[0][0])
#                     results.append(int(res[0][0]))
            except Exception:
                self.utils.print_both("Unexpected error from 'count_sequence_pdr_info':", sys.exc_info()[0])
                raise
        return results
#
#     def count_sequence_pdr_info_ill(self):
#         results = {}
#         primer_suites = self.get_primer_suite_name()
#         lane          = self.get_lane().pop()
#         for primer_suite in primer_suites:
#             primer_suite_lane = primer_suite + ", lane " + lane
#             my_sql = """SELECT count(%s_id)
#                         FROM %s
#                           JOIN run_info_ill USING(run_info_ill_id)
#                           JOIN run USING(run_id)
#                           JOIN primer_suite using(primer_suite_id)
#                         WHERE run = '%s'
#                           AND lane = %s
#                           AND primer_suite = '%s';
#                           """ % (self.table_names["sequence_pdr_info_table_name"], self.table_names["sequence_pdr_info_table_name"], self.rundate, lane, primer_suite)
#             res    = self.my_conn.execute_fetch_select(my_sql)
#             try:
#                 if (int(res[0][0]) > 0):
#                     results[primer_suite_lane] = int(res[0][0])
# #                     results.append(int(res[0][0]))
#             except Exception:
#                 self.utils.print_both("Unexpected error from 'count_sequence_pdr_info_ill':", sys.exc_info()[0])
#                 raise
#         return results
# #             int(res[0][0])
#
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
#         print self.suffix_used
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
        file_seq_db_counts   = self.count_sequence_pdr_info()
        file_seq_orig_count = self.count_seq_from_files_grep()

        for pr_suite, file_seq_db_count in file_seq_db_counts.items():
            if (file_seq_orig_count == file_seq_db_count):
                self.utils.print_both("All sequences from files made it to %s for %s %s: %s == %s\n" % (self.db_server, self.rundate, pr_suite, file_seq_orig_count, file_seq_db_count))
            else:
                self.utils.print_both("Warning: Amount of sequences from files not equal to the one in the db for %s %s: %s != %s\n" % (self.rundate, pr_suite, file_seq_orig_count, file_seq_db_count))
#
#                 ("Oops, amount of sequences from files not equal to the one in the db for %.\nIn file: %s != in db: %s\n==============" % (pr_suite, file_seq_orig_count, file_seq_db_count))

    def put_seq_statistics_in_file(self, filename, seq_in_file):
#        if os.path.exists(file_full):
#            os.remove(file_full)
        self.utils.write_seq_frequencies_in_file(self.unique_file_counts, filename, seq_in_file)

    def insert_taxonomy(self):
        # TODO: mv to Taxonomy?
        self.taxonomy.get_taxonomy_from_gast(self.gast_dict)
        if (self.db_server == "vamps2"):
            self.taxonomy.insert_split_taxonomy()
        elif (self.db_server == "env454"):
            self.taxonomy.insert_whole_taxonomy()
            self.taxonomy.get_taxonomy_id_dict()      
    
    def insert_pdr_info(self, run_info_ill_id):
        all_insert_pdr_info_vals = self.seq.prepare_pdr_info_values(run_info_ill_id, self.all_dataset_run_info_dict, self.db_server)

        group_vals = self.utils.grouper(all_insert_pdr_info_vals, len(all_insert_pdr_info_vals))
        sequence_table_name = self.table_names["sequence_table_name"]
        if (self.db_server == "vamps2"):
            fields = "dataset_id, %s_id, seq_count, classifier_id" % sequence_table_name
        elif (self.db_server == "env454"):
            fields = "run_info_ill_id, %s_id, seq_count" % sequence_table_name
        table_name = self.table_names["sequence_pdr_info_table_name"]
        query_tmpl = self.my_conn.make_sql_for_groups(table_name, fields)
        
        logger.debug("insert sequence_pdr_info:")
        self.my_conn.run_groups(group_vals, query_tmpl)

    def insert_sequence_uniq_info(self):
        if (self.db_server == "vamps2"):
            self.insert_silva_taxonomy_info_per_seq()
            self.seq.insert_sequence_uniq_info2()
        elif (self.db_server == "env454"):
            self.seq.insert_sequence_uniq_info_ill(self.gast_dict)

    def insert_silva_taxonomy_info_per_seq(self):

        for fasta_id, gast_entry in self.gast_dict.items():

            curr_seq          = self.seq.fasta_dict[fasta_id]
            sequence_id       = self.seq.seq_id_dict[curr_seq]
            silva_taxonomy_id = self.taxonomy.silva_taxonomy_id_per_taxonomy_dict[gast_entry[0]]
            gast_distance     = gast_entry[1]
            rank_id           = self.taxonomy.all_rank_w_id[gast_entry[2]]
            refssu_id         = 0
            refssu_count      = 0
        # self.silva_taxonomy_info_per_seq_list = [[8559950L, 2436599, '0.03900', 0, 0, 83],...
#         (taxonomy, gast_distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = self.gast_dict
            vals = "(%s,  %s,  '%s',  '%s',  %s,  '%s')" % (sequence_id, silva_taxonomy_id, gast_distance, refssu_id, refssu_count, rank_id)
            self.silva_taxonomy_info_per_seq_list.append(vals)
        fields = "sequence_id, silva_taxonomy_id, gast_distance, refssu_id, refssu_count, rank_id"
        query_tmpl = self.my_conn.make_sql_for_groups("silva_taxonomy_info_per_seq", fields)
        group_vals = self.utils.grouper(self.silva_taxonomy_info_per_seq_list, len(self.silva_taxonomy_info_per_seq_list))
        logger.debug("insert silva_taxonomy_info_per_seq:")
        self.my_conn.run_groups(group_vals, query_tmpl)

class Taxonomy:
    def __init__(self, my_conn):

        self.utils        = PipelneUtils()
        self.my_conn      = my_conn
        self.taxa_content = set()
        self.ranks        = ['domain', 'phylum', 'klass', 'order', 'family', 'genus', 'species', 'strain']
        self.taxa_by_rank = []
        self.all_rank_w_id                       = set()
        self.uniqued_taxa_by_rank_dict           = {}
        self.uniqued_taxa_by_rank_w_id_dict      = {}
        self.tax_id_dict                         = {}
        self.taxa_list_w_empty_ranks_dict        = defaultdict(list)
        self.taxa_list_w_empty_ranks_ids_dict    = defaultdict(list)
        self.silva_taxonomy_rank_list_w_ids_dict = defaultdict(list)
        self.silva_taxonomy_ids_dict             = defaultdict(list)
        self.silva_taxonomy_id_per_taxonomy_dict = defaultdict(list)
        self.get_all_rank_w_id()

    def get_taxonomy_from_gast(self, gast_dict):
        self.taxa_content = set(v[0] for v in gast_dict.values())

    def get_taxonomy_id_dict(self):
        my_sql = "SELECT %s, %s FROM %s;" % ("taxonomy_id", "taxonomy", "taxonomy")
        res        = self.my_conn.execute_fetch_select(my_sql)
        one_tax_id_dict = dict((y, int(x)) for x, y in res)
        self.tax_id_dict.update(one_tax_id_dict)

    def insert_whole_taxonomy(self):
        val_tmpl   = "('%s')"
        all_taxonomy = set([val_tmpl % taxonomy.rstrip() for taxonomy in self.taxa_content])
        group_vals = self.utils.grouper(all_taxonomy, len(all_taxonomy))
        query_tmpl = self.my_conn.make_sql_for_groups("taxonomy", "taxonomy")
        logger.debug("insert taxonomy:")
        self.my_conn.run_groups(group_vals, query_tmpl)

    def insert_split_taxonomy(self):
        self.parse_taxonomy()
        self.get_taxa_by_rank()
        self.make_uniqued_taxa_by_rank_dict()
#         if (args.do_not_insert == False):
        self.insert_taxa()
        self.silva_taxonomy()
#         if (args.do_not_insert == False):
        self.insert_silva_taxonomy()
        self.get_silva_taxonomy_ids()
        self.make_silva_taxonomy_id_per_taxonomy_dict()

    def parse_taxonomy(self):
        self.taxa_list_dict = {taxon_string: taxon_string.split(";") for taxon_string in self.taxa_content}
        self.taxa_list_w_empty_ranks_dict = {taxonomy: tax_list + [""] * (len(self.ranks) - len(tax_list)) for taxonomy, tax_list in self.taxa_list_dict.items()}

    def get_taxa_by_rank(self):
        self.taxa_by_rank = zip(*self.taxa_list_w_empty_ranks_dict.values())

    def make_uniqued_taxa_by_rank_dict(self):
        for rank in self.ranks:
            rank_num = self.ranks.index(rank)
            uniqued_taxa_by_rank = set(self.taxa_by_rank[rank_num])
            try:
                self.uniqued_taxa_by_rank_dict[rank] = uniqued_taxa_by_rank
            except:
                raise

    def insert_taxa(self):
        for rank, uniqued_taxa_by_rank in self.uniqued_taxa_by_rank_dict.items():
            insert_taxa_vals = '), ('.join(["'%s'" % key for key in uniqued_taxa_by_rank])

            shielded_rank_name = self.shield_rank_name(rank)
            rows_affected = self.my_conn.execute_insert(shielded_rank_name, shielded_rank_name, insert_taxa_vals)
#             self.utils.print_array_w_title(rows_affected, "rows affected by self.my_conn.execute_insert(%s, %s, insert_taxa_vals)" % (rank, rank))

    def shield_rank_name(self, rank):
        return "`"+rank+"`"

    def get_all_rank_w_id(self):
        all_rank_w_id = self.my_conn.get_all_name_id("rank")

        try:
            klass_id = self.utils.find_val_in_nested_list(all_rank_w_id, "klass")
        except:
            raise
        if not klass_id:
            klass_id = self.utils.find_val_in_nested_list(all_rank_w_id, "class")
        l = list(all_rank_w_id)
        l.append(("class", klass_id[0]))
        self.all_rank_w_id = dict((x, y) for x, y in set(l))
        # (('domain', 78), ('family', 82), ('genus', 83), ('klass', 80), ('NA', 87), ('order', 81), ('phylum', 79), ('species', 84), ('strain', 85), ('superkingdom', 86))



    def make_uniqued_taxa_by_rank_w_id_dict(self):
        # self.utils.print_array_w_title(self.uniqued_taxa_by_rank_dict, "===\nself.uniqued_taxa_by_rank_dict from def silva_taxonomy")

        for rank, uniqued_taxa_by_rank in self.uniqued_taxa_by_rank_dict.items():
            shielded_rank_name = self.shield_rank_name(rank)
            taxa_names         = ', '.join(["'%s'" % key for key in uniqued_taxa_by_rank])
            taxa_w_id          = self.my_conn.get_all_name_id(shielded_rank_name, rank + "_id", shielded_rank_name, 'WHERE %s in (%s)' % (shielded_rank_name, taxa_names))
            self.uniqued_taxa_by_rank_w_id_dict[rank] = taxa_w_id

    def insert_silva_taxonomy(self):
        all_insert_st_vals = []
        
        for arr in self.taxa_list_w_empty_ranks_ids_dict.values():
            insert_dat_vals = ', '.join("'%s'" % key for key in arr)
            all_insert_st_vals.append('(%s)' % insert_dat_vals)
                    
        fields = "domain_id, phylum_id, klass_id, order_id, family_id, genus_id, species_id, strain_id"
        query_tmpl = self.my_conn.make_sql_for_groups("silva_taxonomy", fields)
        group_vals = self.utils.grouper(all_insert_st_vals, len(all_insert_st_vals))
        logger.debug("insert silva_taxonomy:")
        self.my_conn.run_groups(group_vals, query_tmpl)


    def silva_taxonomy(self):
        # silva_taxonomy (domain_id, phylum_id, klass_id, order_id, family_id, genus_id, species_id, strain_id)
        self.make_uniqued_taxa_by_rank_w_id_dict()
        silva_taxonomy_list = []

        for taxonomy, tax_list in self.taxa_list_w_empty_ranks_dict.items():
            # ['Bacteria', 'Proteobacteria', 'Deltaproteobacteria', 'Desulfobacterales', 'Nitrospinaceae', 'Nitrospina', '', '']
            silva_taxonomy_sublist = []
            for rank_num, taxon in enumerate(tax_list):
                rank     = self.ranks[rank_num]
                taxon_id = int(self.utils.find_val_in_nested_list(self.uniqued_taxa_by_rank_w_id_dict[rank], taxon)[0])
                silva_taxonomy_sublist.append(taxon_id)
                # self.utils.print_array_w_title(silva_taxonomy_sublist, "===\nsilva_taxonomy_sublist from def silva_taxonomy: ")
            self.taxa_list_w_empty_ranks_ids_dict[taxonomy] = silva_taxonomy_sublist
            # self.utils.print_array_w_title(self.taxa_list_w_empty_ranks_ids_dict, "===\ntaxa_list_w_empty_ranks_ids_dict from def silva_taxonomy: ")

    def make_silva_taxonomy_rank_list_w_ids_dict(self):
        for taxonomy, silva_taxonomy_id_list in self.taxa_list_w_empty_ranks_ids_dict.items():
            rank_w_id_list = []
            for rank_num, taxon_id in enumerate(silva_taxonomy_id_list):
                rank = self.ranks[rank_num]
                t = (rank, taxon_id)
                rank_w_id_list.append(t)

            self.silva_taxonomy_rank_list_w_ids_dict[taxonomy] = rank_w_id_list
            # self.utils.print_array_w_title(self.silva_taxonomy_rank_list_w_ids_dict, "===\nsilva_taxonomy_rank_list_w_ids_dict from def make_silva_taxonomy_rank_list_w_ids_dict: ")
            """
            {'Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhodobiaceae;Rhodobium': [('domain', 2), ('phylum', 2016066), ('klass', 2085666), ('order', 2252460), ('family', 2293035), ('genus', 2303053), ('species', 1), ('strain', 2148217)], ...
            """

    def make_rank_name_id_t_id_str(self, rank_w_id_list):
        a = ""
        for t in rank_w_id_list[:-1]:
            a += t[0] + "_id = " + str(t[1]) + " AND\n"
        a += rank_w_id_list[-1][0] + "_id = " + str(rank_w_id_list[-1][1]) + "\n"
        return a

    def make_silva_taxonomy_ids_dict(self, silva_taxonomy_ids):
        for ids in silva_taxonomy_ids:
#             ids[-1] = silva_taxonomy_id, the rest are ids for each rank
            self.silva_taxonomy_ids_dict[int(ids[-1])] = [int(id) for id in ids[0:-1]]
        # self.utils.print_array_w_title(self.silva_taxonomy_ids_dict, "===\nsilva_taxonomy_ids_dict from def get_silva_taxonomy_ids: ")
        # {2436595: [2, 2016066, 2085666, 2252460, 2293035, 2303053, 1, 2148217], 2436596: [...

    def get_silva_taxonomy_ids(self):
        self.make_silva_taxonomy_rank_list_w_ids_dict()

        sql_part = ""
        for taxonomy, rank_w_id_list in self.silva_taxonomy_rank_list_w_ids_dict.items()[:-1]:
            a = self.make_rank_name_id_t_id_str(rank_w_id_list)
            sql_part += "(%s) OR " % a

        a_last = self.make_rank_name_id_t_id_str(self.silva_taxonomy_rank_list_w_ids_dict.values()[-1])
        sql_part += "(%s)" % a_last

        field_names = "domain_id, phylum_id, klass_id, order_id, family_id, genus_id, species_id, strain_id"
        table_name  = "silva_taxonomy"
        where_part  = "WHERE " + sql_part
        silva_taxonomy_ids = self.my_conn.get_all_name_id(table_name, "", field_names, where_part)

        """
        ((2436595L, 2L, 2016066L, 2085666L, 2252460L, 2293035L, 2303053L, 1L, 2148217L), ...
        """
        self.make_silva_taxonomy_ids_dict(silva_taxonomy_ids)

    def make_silva_taxonomy_id_per_taxonomy_dict(self):
        for silva_taxonomy_id, st_id_list1 in self.silva_taxonomy_ids_dict.items():
            taxon_string = self.utils.find_key_by_value_in_dict(self.taxa_list_w_empty_ranks_ids_dict.items(), st_id_list1)
            self.silva_taxonomy_id_per_taxonomy_dict[taxon_string[0]] = silva_taxonomy_id
        # self.utils.print_array_w_title(self.silva_taxonomy_id_per_taxonomy_dict, "silva_taxonomy_id_per_taxonomy_dict from silva_taxonomy_info_per_seq = ")

class Seq:
    def __init__(self, taxonomy, table_names):

        self.utils       = PipelneUtils()
        self.taxonomy    = taxonomy
        self.my_conn     = self.taxonomy.my_conn
        self.table_names = table_names
        self.seq_id_dict = {}
        self.fasta_dict  = {}

        self.sequences   = ""
        self.taxa        = ""
        self.refhvr_id   = ""
        self.the_rest    = ""

    def prepare_fasta_dict(self, filename):
        read_fasta = fastalib.ReadFasta(filename)
        self.fasta_dict = dict(zip(read_fasta.ids, read_fasta.sequences))
        read_fasta.close()

    def make_seq_upper(self, filename):
        sequences  = [seq.upper() for seq in self.fasta_dict.values()] #here we make uppercase for VAMPS compartibility
        return list(set(sequences))

    def insert_seq(self, sequences):
        val_tmpl = "(COMPRESS('%s'))"
        all_seq = set([val_tmpl % seq for seq in sequences])
        group_vals = self.utils.grouper(all_seq, len(all_seq))
        query_tmpl = self.my_conn.make_sql_for_groups(self.table_names["sequence_table_name"], self.table_names["sequence_field_name"])
        logger.debug("insert sequences:")
        self.my_conn.run_groups(group_vals, query_tmpl)

    def get_seq_id_dict(self, sequences):
#         TODO: ONCE IN CLASS
        sequence_field_name = self.table_names["sequence_field_name"]
        sequence_table_name = self.table_names["sequence_table_name"]
        id_name    = self.table_names["sequence_table_name"] + "_id"
        query_tmpl = """SELECT %s, uncompress(%s) FROM %s WHERE %s in (COMPRESS(%s))"""
        val_tmpl   = "'%s'"
        try:
            group_seq = self.utils.grouper(sequences, len(sequences))
            for group in group_seq:
                seq_part = '), COMPRESS('.join([val_tmpl % key for key in group])
                my_sql     = query_tmpl % (id_name, sequence_field_name, sequence_table_name, sequence_field_name, seq_part)
                res        = self.my_conn.execute_fetch_select(my_sql)
                one_seq_id_dict = dict((y.upper(), int(x)) for x, y in res)
                self.seq_id_dict.update(one_seq_id_dict)
        except:
            if len(sequences) == 0:
                self.utils.print_both(("ERROR: There are no sequences, please check if there are correct fasta files in the directory %s") % self.fasta_dir)
            raise
        
    def prepare_pdr_info_values(self, run_info_ill_id, all_dataset_run_info_dict, db_server):
        all_insert_pdr_info_vals = []
        for fasta_id, seq in self.fasta_dict.items():
            if (not run_info_ill_id):
                self.utils.print_both("ERROR: There is no run info yet, please check if it's uploaded to env454")
#             seq_upper = seq.upper()
            sequence_id = self.seq_id_dict[seq]
        
            seq_count = int(fasta_id.split('|')[-1].split(':')[-1])
        
            if (db_server == "vamps2"):       
                dataset_id = all_dataset_run_info_dict[run_info_ill_id]
                vals = "(%s, %s, %s, %s)" % (dataset_id, sequence_id, seq_count, C.classifier_id)
            elif (db_server == "env454"):
                vals = "(%s, %s, %s)" % (run_info_ill_id, sequence_id, seq_count)
         
            all_insert_pdr_info_vals.append(vals)
        
        return all_insert_pdr_info_vals
          
    def get_seq_id_w_silva_taxonomy_info_per_seq_id(self):
        sequence_ids_strs = [str(i) for i in self.seq_id_dict.values()]
        where_part = 'WHERE sequence_id in (%s)' % ', '.join(sequence_ids_strs)
        self.seq_id_w_silva_taxonomy_info_per_seq_id = self.my_conn.get_all_name_id("silva_taxonomy_info_per_seq", "silva_taxonomy_info_per_seq_id", "sequence_id", where_part)

    def insert_sequence_uniq_info2(self):
        self.get_seq_id_w_silva_taxonomy_info_per_seq_id()
        fields = "sequence_id, silva_taxonomy_info_per_seq_id"
        sequence_uniq_info_values = ["(%s,  %s)"  % (i1, i2) for i1, i2 in self.seq_id_w_silva_taxonomy_info_per_seq_id]
        query_tmpl = self.my_conn.make_sql_for_groups("sequence_uniq_info", fields)
        group_vals = self.utils.grouper(sequence_uniq_info_values, len(sequence_uniq_info_values))
        logger.debug("insert sequence_uniq_info_ill:")
        self.my_conn.run_groups(group_vals, query_tmpl)

    def insert_sequence_uniq_info_ill(self, gast_dict):
        all_insert_sequence_uniq_info_ill_vals = []
        for fasta_id, gast in gast_dict.items():
            (taxonomy, distance, rank, refssu_count, vote, minrank, taxa_counts, max_pcts, na_pcts, refhvr_ids) = gast
            seq  = self.fasta_dict[fasta_id]
            sequence_id = self.seq_id_dict[seq]
            rank_id = self.taxonomy.all_rank_w_id[rank]
            if taxonomy in self.taxonomy.tax_id_dict:
                taxonomy_id = self.taxonomy.tax_id_dict[taxonomy]
                vals = """(%s,  %s,  '%s',  '%s',  %s,  '%s')
                        """ % (sequence_id, taxonomy_id, distance, refssu_count, rank_id, refhvr_ids.rstrip())
                all_insert_sequence_uniq_info_ill_vals.append(vals)
        group_vals = self.utils.grouper(all_insert_sequence_uniq_info_ill_vals, len(all_insert_sequence_uniq_info_ill_vals))
        
        fields = "%s_id, taxonomy_id, gast_distance, refssu_count, rank_id, refhvr_ids" % (self.table_names["sequence_table_name"])
        query_tmpl = self.my_conn.make_sql_for_groups("sequence_uniq_info_ill", fields)
        
        logger.debug("insert sequence_uniq_info_ill:")
        self.my_conn.run_groups(group_vals, query_tmpl)


