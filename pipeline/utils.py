import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))
import subprocess
import constants as C
import time
import datetime

from contextlib import closing
import zipfile
import zlib
import collections
sys.path.append('/bioware/linux/seqinfo/bin/python_pipeline/py_mbl_sequencing_pipeline')
from pipeline.pipelinelogging import logger
from subprocess import call
import getpass
import math


def it_is_py3():
    if sys.version_info[0] < 3:
        return False
    if sys.version_info[0] >= 3:
        return True
        
if it_is_py3():
    import string
    base_complement_translator = bytes.maketrans(b"ACGTRYMK", b"TGCAYRKM")
    from itertools import zip_longest as izip_longest
else:
    from string import maketrans
    base_complement_translator = maketrans("ACGTRYMK", "TGCAYRKM")
    from itertools import izip_longest as izip_longest

# the json expected files get loaded and parsed into Unicode strings
# but the asserts won't work comparing unicode to ascii so we need change them
# to plain strings
def convert_unicode_dictionary_to_str(data):
    if isinstance(data, unicode):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(convert_unicode_dictionary_to_str, data.items()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert_unicode_dictionary_to_str, data))
    else:
        return data


########################################################
def check_for_Ns(seq):
    """Doc string here.."""
    nCount = seq.count('N')
    if( nCount > 0 ):
        return nCount
    else:
        return 0

def remove_runkey(seq,runkeys):
    """Doc string here.."""
    found_key = ''
    for key in runkeys:
        if (seq.find(key) == 0):
            found_key = key
            seq = seq[len(key):]
            break
        else:
            continue
    return found_key, seq

def find_sequence_direction( direction = '' ):
    """Doc string here.."""
    seqFwd = 0
    seqRev = 0
    if(direction == 'F'):
        seqFwd = 1
    if(direction == "R"):
        seqRev = 1

    if ( not direction or (seqFwd + seqRev) == 0 ):
        return 0
    elif ( (seqFwd + seqRev) == 2 ):
        return "B"
    elif (seqFwd):
        return "F"
    else:
        return "R"

def check_for_quality(rawseq, trimseq, quality_scores):
    """Doc string here.."""
    start = rawseq.find(trimseq)
    end   = start + len(trimseq)
    #logger.debug('IN Check_for_quality ' + str(quality_scores))
    scores = quality_scores[start:end]

    return sum(scores,0.0) / len(scores)

def revcomp(sequence):
    reversed = str(sequence[::-1])
    return reversed.translate(base_complement_translator)

def set_trim1():
    return (True,False,False,False)

def set_trim2():
    return (True,True,False,False)

def set_trim3():
    return (True,True,True,False)

def set_chim4():
    return (False,True,False,False)

def set_chim5():
    return (False,True,True,False)

def set_chim6():
    return (False,True,True,True)

def set_gast7():
    return (False,False,True,False)

def set_gast8():
    return (False,False,True,True)

def set_vamps9():
    return (False,False,False,True)

def set_all10():
    return (True,True,True,True)

options = {
        1 : set_trim1,
        2 : set_trim2,
        3 : set_trim3,
        4 : set_chim4,
        5 : set_chim5,
        6 : set_chim6,
        7 : set_gast7,
        8 : set_gast8,
        9 : set_vamps9,
        10 : set_all10,
}

def wait_for_cluster_to_finish(my_running_id_list):
    #print('My IDs',running_id_list)
    logger.debug('Max run time set to ' + str(C.cluster_max_wait) + ' seconds')
    logger.debug('These are my running qsub IDs ' + str(my_running_id_list))
    my_working_id_list = my_running_id_list

    counter =0

    time.sleep(C.cluster_initial_check_interval)

    while my_working_id_list:

        qstat_codes = get_qstat_id_list()
        if not qstat_codes['id']:
            #print('No qstat ids')
            logger.debug("id list not found: may need to increase initial_interval if you haven't seen running ids.")
            return ('SUCCESS','id list not found','',)
        if 'Eqw' in qstat_codes['code']:
            logger.debug( "Check cluster: may have error code(s), but they may not be mine!")


        got_one = False

        #print('working ids',my_working_id_list)
        if my_working_id_list[0] in qstat_codes['id']:

            got_one = True
            name = qstat_codes['name'][qstat_codes['id'].index(my_working_id_list[0])]
            user = qstat_codes['user'][qstat_codes['id'].index(my_working_id_list[0])]
            code = qstat_codes['code'][qstat_codes['id'].index(my_working_id_list[0])]


            if code == 'Eqw':
                return ('FAIL','Found Eqw code',my_working_id_list[0])
            elif code == 'qw':
                logger.debug("id is still queued: " +  str(my_working_id_list[0]) + " " + str(code))
            elif code == 'r':
                logger.debug("id is still running: " + str(my_working_id_list[0]) + " " + str(code))
            else:
                logger.debug('Unknown qstat code ' + str(code))
        else:
            my_working_id_list = my_working_id_list[1:]
            logger.debug('id finished ' + str(my_working_id_list))

        if not my_working_id_list:
            return ('SUCCESS','not my_working_id_list','')
        #if not got_one:
            #print('IN not got one',)
        #    return ('SUCCESS','not got one','')

        time.sleep(C.cluster_check_interval)
        counter = counter + C.cluster_check_interval
        if counter >= C.cluster_max_wait:
            return ('FAIL','Max Time exceeded',C.cluster_max_wait)

    return ('FAIL','Unknown','Unknown')

def get_qstat_id_list():

    # ['139239', '0.55500', 'usearch', 'avoorhis', 'r', '01/22/2012', '09:00:39', 'all.q@grendel-07.bpcservers.pr', '1']
    # 1) id
    # 2)
    # 3) name
    # 4) username
    # 5) code r=running, Ew=Error
    qstat_cmd = 'qstat'
    qstat_codes={}
    output = subprocess.check_output(qstat_cmd)
    #print(output)
    output_list = output.strip().split("\n")[2:]
    qstat_codes['id'] = [n.split()[0] for n in output_list]
    qstat_codes['name'] = [n.split()[2] for n in output_list]
    qstat_codes['user'] = [n.split()[3] for n in output_list]
    qstat_codes['code'] = [n.split()[4] for n in output_list]
    #print('Found IDs',qstat_ids)



    return qstat_codes


def find_key(dic, val):
    """return the first key of dictionary dic given the value"""
    return [k for k, v in dic.items() if v == val][0]

def mysort(uniques,names):
    """ Sorts the uniques using the uniques and names hashes:

    uniques[lane_tag][trimmed_sequence] = read_id
    names[lane_tag][read_id1] = [read_id1, read_id2, read_id3, read_id4]

    returns a list of tuples (read_id, count, sequence) highest to lowest
    """
    sorted_uniques = []

    # sorted_names should be list of ids with the highest number at the top
    sorted_names = sorted(names.items(), key=lambda x: len(x[1]), reverse=True)

    for n in sorted_names:

        seq = find_key(uniques, n[0])
        sorted_uniques.append( (n[0], len(n[1]),seq) )

    return sorted_uniques

def extract_zipped_file(run_date, outdir, filename):
    """

    """
    # check if zipped
    assert os.path.isdir(outdir)
    archivename = os.path.join(outdir,run_date+'.zip')
    if zipfile.is_zipfile(archivename):
        zf = zipfile.ZipFile(archivename, 'r')

        try:
            data = zf.read(filename)
        except KeyError:
            logger.error('ERROR: Did not find %s in zip file' % filename)
        else:
            logger.error(filename, ':')
            logger.error(repr(data))
        print
        zf.close()
    else:
        logger.error("No zipfile archive found:",archivename)

def zip_up_directory(run_date, dirPath, mode='a'):
    """
    This should be run at the end of each process to zip the files in each directory
    """
    files_to_compress = ['fa','db','names','sff','fasta','fastq']
    assert os.path.isdir(dirPath)
    zipFilePath = os.path.join(dirPath,run_date+'.zip')

    zf = zipfile.ZipFile(zipFilePath, mode)

    for (archiveDirPath, dirNames, fileNames) in os.walk(dirPath):
        for file in fileNames:
            if file.split('.')[-1] in files_to_compress:
                filePath = os.path.join(dirPath, file)
                zf.write(filePath, compress_type=zipfile.ZIP_DEFLATED)

    zipInfo = zipfile.ZipInfo(zipFilePath)

    for i in zf.infolist():
        dt = datetime.datetime(*(i.date_time))
        logger.debug("%s\tSize: %sb\tCompressed: %sb\t\tModified: %s" % (i.filename, i.file_size, i.compress_size, dt.ctime()))
        os.remove(i.filename)

    zf.close()

def write_status_to_vamps_db(site='vampsdev', id='0', status='Test', message=''):
    """
    This should be available to write status updates to vamps:vamps_upload_status.
    It is especially important for MoBeDAC uploads because the qiime site
    will 'see' and react to the message in the db.  <-- not true any longer 2014-02-01 AAV


    """
    import ConMySQL
    from pipeline.db_upload import MyConnection
    today   = str(datetime.date.today())
    if site == 'vamps':
        db_host    = 'vampsdb'
        db_name    = 'vamps'
        db_home = '/xraid2-2/vampsweb/vamps/'
    else:
        db_host    = 'bpcweb7'
        db_name    = 'vamps'
        db_home = '/xraid2-2/vampsweb/vampsdev/'
    #obj=ConMySQL.New(db_host, db_name, db_home)
    #my_conn = MyConnection(db_host, db_name)
    obj=ConMySQL.New(db_host, db_name, db_home)
    conn = obj.get_conn()
    cursor = conn.cursor()
    query = "update vamps_upload_status set status='%s', status_message='%s', date='%s' where id='%s'" % (status, message, today, id)
    try:
        cursor.execute(query)
        #print("executing",query)
    except:
        conn.rollback()
        logger.error("ERROR status update failed")
    else:
        conn.commit()

class PipelneUtils:
    def __init__(self):
        pass

    def find_val_in_nested_list(self, hey, needle):
        return [v for k, v in hey if k.lower() == needle.lower()]

    def print_array_w_title(self, message, title = 'message'):
        print(title)
        print(message)
        
    def grouper(self, iterable, obj_len, fillvalue=None):
        n = 10 ** self.magnitude(obj_len)
        args = [iter(iterable)] * n
        return izip_longest(*args, fillvalue=fillvalue)        

    def find_val_in_nested_list(self, hey, needle):
        return [v for k, v in hey if k.lower() == needle.lower()]
  
    def find_key_by_value_in_dict(self, hey, needle):
        return [k for k, v in hey if v == needle]

    def make_insert_values(self, matrix):
        all_insert_vals = ""
        
        for arr in matrix[:-1]:
            insert_dat_vals = ', '.join("'%s'" % key for key in arr)
            all_insert_vals += insert_dat_vals + "), ("
        
        all_insert_vals += ', '.join("'%s'" % key for key in matrix[-1])
        
        # self.print_array_w_title(all_insert_vals, "all_insert_vals from make_insert_values")
        return all_insert_vals


    def call_sh_script(self, script_name_w_path, where_to_run):
        try:
            call(['chmod', '0774', script_name_w_path])
            if self.is_local():
                self.print_both("call(['qsub', script_name_w_path], cwd=(where_to_run))")
                call(['bash', script_name_w_path], cwd=(where_to_run))
            else:
                call(['qsub', script_name_w_path], cwd=(where_to_run))
#             pass
        except:
            self.print_both("Problems with script_name = %s or qsub" % (script_name_w_path))
            raise

    def make_users_email(self):
        username = getpass.getuser()
        return username + "@mbl.edu"

    def open_write_close(self, file_name, text):
        ini_file = open(file_name, "w")
        ini_file.write(text)
        ini_file.close()

    def create_job_array_script(self, command_line, dir_to_run, files_list):
        files_string         = " ".join(files_list)
        files_list_size         = len(files_list)
        command_file_name = os.path.basename(command_line.split(" ")[0])
        script_file_name  = command_file_name + "_" + self.runobj.run + "_" + self.runobj.lane_name + ".sh"
        script_file_name_full = os.path.join(dir_to_run, script_file_name)
        log_file_name     = script_file_name + ".sge_script.sh.log"
        email_mbl         = self.make_users_email()
        text = (
                '''#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N %s
# Giving the name of the output log file
#$ -o %s
# Combining output/error messages into one file
#$ -j y
# Send mail to these users
#$ -M %s
# Send mail at job end (e); -m as sends abort, suspend.
#$ -m as
#$ -t 1-%s
# Now the script will iterate %s times.

  file_list=(%s)

  i=$(expr $SGE_TASK_ID - 1)
#   echo "i = $i"
  # . /etc/profile.d/modules.sh
  # . /xraid/bioware/bioware-loader.sh
  . /xraid/bioware/Modules/etc/profile.modules
  module load bioware

  echo "%s ${file_list[$i]}"
  %s ${file_list[$i]}
''' % (script_file_name, log_file_name, email_mbl, files_list_size, files_list_size, files_string, command_line, command_line)
# ''' % (script_file_name, log_file_name, email_mbl, files_list_size, files_list_size, files_string, command_line)
                )
        self.open_write_close(script_file_name_full, text)
        return script_file_name
    
    def read_file(self, source_name):
        with open(source_name, "r") as sources:
            return sources.readlines()

    def get_vsearch_version(self):
        import commands
        return commands.getstatusoutput('vsearch --help | head -1')

    def flatten_list_of_lists(self, mylist):
        return [item for sublist in mylist for item in sublist]


    def write_seq_frequencies_in_file(self, out_file, fa_file_name, seq_in_file):
        try:
            with open(out_file, "a") as myfile:
                myfile.write(str(fa_file_name) + ": " + str(seq_in_file) + "\n")
        except Exception:
            logger.error(Exception)

    def is_local(self):
        logger.debug(os.uname()[1])
        dev_comps = ['ashipunova.mbl.edu', "as-macbook.home", "as-macbook.local", "Ashipunova.local", "Annas-MacBook-new.local", "Annas-MacBook.local"]
        if os.uname()[1] in dev_comps:
            return True
        else:
            return False

    def is_vamps(self):
        logger.debug(os.uname()[1])
        dev_comps = ['bpcweb8','bpcweb7','bpcweb7.bpcservers.private', 'bpcweb8.bpcservers.private', 'vampsdev', 'vampsdb']
        if os.uname()[1] in dev_comps:
            return True
        else:
            return False

    def find_in_nested_dict(self, nested_dict, to_find):
        reverse_linked_q = list()
        reverse_linked_q.append((list(), nested_dict))
        while reverse_linked_q:
            this_key_chain, this_v = reverse_linked_q.pop()
            # finish search if found the mime type
            if this_v == to_find:
                return this_key_chain
            # not found. keep searching
            # queue dicts for checking / ignore anything that's not a dict
            try:
                items = this_v.items()
            except AttributeError:
                continue  # this was not a nested dict. ignore it
            for k, v in items:
                reverse_linked_q.append((this_key_chain + [k], v))
        # if we haven't returned by this point, we've exhausted all the contents

        return False



    def check_if_array_job_is_done(self, job_name):
        cluster_done = False
        check_qstat_cmd_line = "qstat -r | grep %s | wc -l" % job_name
        logger.debug("check_qstat_cmd_line = %s" % check_qstat_cmd_line)
        try:
            p = subprocess.Popen(check_qstat_cmd_line, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            num_proc = int(output)
            logger.debug("qstat is running %s '%s' processes" % (num_proc, job_name))
    #         pprint(p)

            if (num_proc == 0):
                cluster_done = True
    #         print("cluster_done from check_if_cluster_is_done = %s" % cluster_done)
        except:
            logger.error("%s can be done only on a cluster." % job_name)
            raise
        return cluster_done

    def run_until_done_on_cluster(self, job_name):
        start = time.time()
        time_before = self.get_time_now()
        logger.debug("time_before = %s" % time_before)
        logger.debug("Waiting for the cluster...")
        while True:
            if self.is_local():
                time.sleep(1)
            else:
                time.sleep(120)
            cluster_done = self.check_if_array_job_is_done(job_name)
            logger.debug("cluster_done = %s" % cluster_done)
            if (cluster_done):
                break

        elapsed = (time.time() - start)
        logger.debug("Cluster is done with %s in: %s" % (job_name, elapsed))

    def get_time_now(self):
        """date and hour only!"""
        return time.strftime("%m/%d/%Y %H:%M", time.localtime())
# '2009-01-05 22'

    def print_both(self, message):
        print(message)
        logger.debug(message)
        
    def magnitude(self, x):
        return int(math.log10(x))        

class Dirs:
    """get input dir from args, create all other dirs
input_dir - directory with fastq or sff files
Output path example: /xraid2-2/g454/run_new_pipeline/illumina/miseq/20121025/analysis/gast
dir_prefix is <vamps_user>_<random_number> for vamps user uplads
and run_date otherwise
example of initiation all directories in pipelie_ui
example of getting all directory name in illumina_files
"""
    def __init__(self, is_user_upload, dir_prefix, platform, lane_name = '', site = ''):
        self.utils             = PipelneUtils()
        self.output_dir_name   = None
        self.get_path(is_user_upload, dir_prefix, platform, lane_name, site)
        self.analysis_dir      = os.path.join(self.output_dir,   C.subdirs['analysis_dir'])
        self.gast_dir          = os.path.join(self.analysis_dir, C.subdirs['gast_dir'])
        self.reads_overlap_dir = os.path.join(self.analysis_dir, C.subdirs['reads_overlap_dir'])
        self.vamps_upload_dir  = os.path.join(self.analysis_dir, C.subdirs['vamps_upload_dir'])
        self.chimera_dir       = os.path.join(self.analysis_dir, C.subdirs['chimera_dir'])
        self.trimming_dir      = os.path.join(self.analysis_dir, C.subdirs['trimming_dir'])

        self.unique_file_counts = os.path.join(self.analysis_dir, "unique_file_counts")


    def check_and_make_dir(self, dir_name):
#        if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except OSError:
            if os.path.isdir(dir_name):
                logger.error("\nDirectory %s already exists."  % (dir_name))
#                 confirm_msg = "Do you want to continue? (Yes / No) "
#                 answer = raw_input(confirm_msg)
#                 if answer != 'Yes':
#                     sys.exit("There was an error in the directory " + dir_name + " creation - Exiting.")
#                 elif answer == 'Yes':
                pass
            else:
            # There was an error on creation, so make sure we know about it
                raise
        return dir_name

    def check_dir(self, dir_name):
        if os.path.isdir(dir_name):
            return dir_name
        else:
            return self.check_and_make_dir(dir_name)

    def get_path(self, is_user_upload, dir_prefix, platform, lane_name, site):
        """
        dir_prefix is <vamps_user>_<random_number> for vamps user uploads
        and run_date otherwise
        """
        if is_user_upload:
            if site == 'vamps':
                root_dir  = C.output_root_vamps_users
            elif site == 'new_vamps':
                root_dir  = C.output_root_newvamps_users
            else:
                root_dir  = C.output_root_vampsdev_users

            self.output_dir = os.path.join(root_dir, dir_prefix)

        else:
            id_number = dir_prefix
#            runobj['run']
            root_dir  = C.output_root_mbl

            if self.utils.is_local():
                root_dir  = C.output_root_mbl_local

            self.output_dir = os.path.join(root_dir, platform, id_number)
            if (lane_name != ''):
                self.output_dir = os.path.join(root_dir, platform, id_number, lane_name)


    def check_and_make_output_dir(self):
        self.check_and_make_dir(self.output_dir)

    def create_all_output_dirs(self):
        self.check_and_make_dir(self.analysis_dir)
        self.check_and_make_dir(self.gast_dir)
        self.check_and_make_dir(self.reads_overlap_dir)
        self.check_and_make_dir(self.vamps_upload_dir)
        self.check_and_make_dir(self.chimera_dir)
        self.check_and_make_dir(self.trimming_dir)

    def create_gast_name_dirs(self, name_iterator):
        gast_name_dirs = []
#         gast_name_dirs1 = [self.check_dir(os.path.join(self.gast_dir, key)) for key in name_iterator]
        for key in name_iterator:
            #gast_name_dirs.append(self.check_and_make_dir(os.path.join(self.gast_dir, key)))
            gast_name_dirs.append(self.check_dir(os.path.join(self.gast_dir, key)))
        return gast_name_dirs

    def delete_file(self, filename):
        try:
            os.remove(filename)
            logger.debug("DELETE %s" % (filename))
        except OSError:
            pass


    def get_all_files(self, walk_dir_name, ext = ""):
        files = {}
        for dirname, dirnames, filenames in os.walk(walk_dir_name, followlinks=True):
            if ext:
                filenames = [f for f in filenames if f.endswith(ext)]

            for file_name in filenames:
                full_name = os.path.join(dirname, file_name)
                (file_base, file_extension) = os.path.splitext(os.path.join(dirname, file_name))
                files[full_name] = (dirname, file_base, file_extension)
    #        print("len(files) = %s" % len(files))
        return files

    def get_all_files_by_ext(self, walk_dir_name, extension):
        return [file for file in os.listdir(walk_dir_name) if file.endswith(extension)]

    def chmod_all(self, dir_name):
      try:
        call(['chmod', '-R', 'ug+w', dir_name])
      except Exception:
        logger.error("call(['chmod', '-R', 'ug+w', %s]) didn't work: \n" % (dir_name))
        logger.error(Exception)
        pass


if __name__=='__main__':
    logger.debug("GTTCAAAGAYTCGATGATTCAC")
    logger.debug(revcomp("GTTCAAAGAYTCGATGATTCAC"))

