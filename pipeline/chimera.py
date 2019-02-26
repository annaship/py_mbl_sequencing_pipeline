import subprocess
import sys, os
import re
import time
from pipeline.pipelinelogging import logger
from pipeline.utils import Dirs, PipelneUtils
from pipeline.utils import *
from pprint import pprint
from collections import defaultdict, namedtuple
# sys.path.append("/xraid/bioware/linux/seqinfo/bin")
# sys.path.append("/Users/ashipunova/bin/illumina-utils")
# sys.path.append("/Users/ashipunova/bin/illumina-utils/illumina-utils/scripts")
# sys.path.append("/bioware/merens-illumina-utils")

import IlluminaUtils.lib.fastqlib as fq
import IlluminaUtils.lib.fastalib as fa

import pipeline.constants as C
import json


class Chimera:
    """ Define here """
    def __init__(self, runobj = None):
        self.utils      = PipelneUtils()
        self.runobj     = runobj
        self.run_keys   = self.runobj.run_keys
        self.rundate    = self.runobj.run
        try:
            self.use_cluster = self.runobj.use_cluster
        except:
            self.use_cluster = True
        self.chg_suffix         = ".chg"
        self.chimeras_suffix    = ".chimeras"
        self.ref_suffix         = ".db"
        self.denovo_suffix      = ".txt"
        self.nonchimeric_suffix = "." + C.nonchimeric_suffix #".nonchimeric.fa"
        self.chimeric_suffix    = ".chimeric.fa"
        self.base_suffix        = "unique" + self.chimeras_suffix
        
        self.cluster_slots = {
          "grendel": [12, 8],
          "cricket": [40],
          "cluster5": [32]          
        }

        try:
            if self.runobj.lane_name:
                lane_name = self.runobj.lane_name
            else:
                lane_name = ''
        except:
            lane_name = ''

        if self.runobj.vamps_user_upload:
            os.environ['SGE_ROOT'] ='/opt/sge'
            os.environ['SGE_CELL'] ='grendel'
            path                   = os.environ['PATH']
            os.environ['PATH']     = '/opt/sge/bin/lx24-amd64:'+path
            site       = self.runobj.site
            dir_prefix = self.runobj.user + '_' + self.runobj.run
            self.dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, lane_name = lane_name, site = site)
            self.idx_keys = convert_unicode_dictionary_to_str(json.loads(open(self.runobj.trim_status_file_name,"r").read()))["new_lane_keys"]
            self.analysis_dir = self.dirs.check_dir(self.dirs.analysis_dir)
            self.indir  = self.dirs.check_dir(self.dirs.trimming_dir)
            self.outdir = self.dirs.check_dir(self.dirs.chimera_dir)
            self.gast_dir = self.dirs.check_dir(self.dirs.gast_dir)

        else:
            site = ''
            dir_prefix = self.runobj.run
            self.dirs = Dirs(self.runobj.vamps_user_upload, dir_prefix, self.runobj.platform, lane_name = lane_name, site = site)
            self.indir  = self.dirs.check_dir(self.dirs.reads_overlap_dir)
            self.outdir = self.dirs.check_dir(self.dirs.chimera_dir)


#         self.usearch_cmd = C.usearch_cmd
        self.usearch_cmd = C.usearch6_cmd
        if self.utils.is_local():
            self.usearch_cmd = C.usearch6_cmd_local
        #self.abskew      = C.chimera_checking_abskew
        self.refdb       = C.chimera_checking_refdb
        if self.utils.is_local():
            self.refdb_local = C.chimera_checking_refdb_local
        self.its_refdb   = C.chimera_checking_its_refdb
        self.input_file_names  = self.make_chimera_input_illumina_file_names()
#         self.output_file_names = self.make_chimera_output_illumina_file_names(self.input_file_names)

    def get_ref_db(self, dna_region):
        ref_db = ''
        if dna_region.upper() == 'ITS':
            ref_db = C.chimera_checking_its_refdb
            logger.debug("got an ITS dna region so using refdb: " + ref_db)
        else:
            ref_db = C.chimera_checking_refdb
            if self.utils.is_local():
                ref_db = C.chimera_checking_refdb_local
            logger.debug("using standard refdb: " + ref_db)
        return ref_db

    def make_chimera_input_illumina_file_names(self):
        input_file_names = {}

        for idx_key in self.run_keys:
            file_name = idx_key + "_" + C.filtered_suffix + ".unique"

            if os.path.exists(os.path.join(self.indir, file_name)):
                input_file_names[idx_key] = file_name

        return input_file_names

    def get_current_dirname(self, in_or_out = ""):
        if in_or_out == "":
            cur_dirname    = self.indir
        else:
            cur_dirname    = self.outdir
        return cur_dirname

    def is_chimera_check_file(self, filename):
        return filename.endswith((self.chimeras_suffix + self.denovo_suffix, self.chimeras_suffix + self.ref_suffix, self.chimeric_suffix, self.nonchimeric_suffix))

    def get_current_filenames(self, cur_dirname):
        cur_file_names = []
        if cur_dirname == self.indir:
            cur_file_names = self.input_file_names.values()
        elif cur_dirname == self.outdir:
            cur_file_names = self.get_chimera_file_names(self.outdir)
        return cur_file_names

    def get_chimera_file_names(self, cur_dirname):
        cur_file_names = []
        for dirname, dirnames, filenames in os.walk(cur_dirname):
            cur_file_names = [filename for filename in filenames if (self.is_chimera_check_file(filename))]
        return cur_file_names

    def read_file(self, source_name):
        with open(source_name, "r") as sources:
            return sources.readlines()

    def illumina_sed(self, lines, target_name, regex, replace, uppercase):
        with open(target_name, "w") as target:
            for line in lines:
                if line.startswith(">"):
                    line1 = regex.sub(replace, line)
                else:
                    if (uppercase):
                        line1 = line.upper()
                    else:
                        line1 = line
                target.write(line1)


    def call_illumina_sed(self, from_to):
        """
            from_to = from_frequency_to_size or from_size_to_frequency
        """
        sed_from_to = namedtuple('sed_from_to', 'find, replace, cur_dirname, cur_file_names, change_from_suffix, change_to_suffix, uppercase')

        from_frequency_to_size = sed_from_to(
        find               = "frequency:",
        replace            = ";size=",
        cur_dirname        = self.indir,
        cur_file_names     = self.get_current_filenames(self.indir),
        change_from_suffix = "",
        change_to_suffix   = self.chg_suffix,
        uppercase          = True
        )

        from_size_to_frequency = sed_from_to(
        find               = ";size=",
        replace            = "frequency:",
        cur_dirname        = self.outdir,
        cur_file_names     = self.get_chimera_file_names(self.outdir),
        change_from_suffix = "",
        change_to_suffix   = "",
        uppercase          = False
        )

        if (from_to == "from_frequency_to_size"):
            tuple_name = from_frequency_to_size
        elif (from_to == "from_size_to_frequency"):
            tuple_name = from_size_to_frequency

        regex          = re.compile(r"%s" % tuple_name.find)
#         logger.debug("find = %s, replace = %s" % (find, replace))
        if (not tuple_name.cur_file_names) and (tuple_name == from_frequency_to_size):
            self.utils.print_both('ERROR: Did not find uniqued files ("%s") in %s, please check if the previous step has finished. Exiting.\n' % (C.filtered_suffix + ".unique" , self.indir))
            sys.exit()
        for cur_file_name in tuple_name.cur_file_names:
            file_name = os.path.join(tuple_name.cur_dirname, cur_file_name)
            source_name = file_name + tuple_name.change_from_suffix
            target_name = file_name + tuple_name.change_to_suffix
            lines = self.read_file(source_name)
            self.illumina_sed(lines, target_name, regex, tuple_name.replace, tuple_name.uppercase)

    def illumina_freq_to_size_in_chg(self):
        find1    = "frequency:"
        replace1 = ";size="
        regex1   = re.compile(r"%s" % find1)

#         logger.debug("cur_file_names: ")
#         pprint(cur_file_names)
        cur_dirname        = self.get_current_dirname()
        cur_file_names     = self.get_current_filenames(cur_dirname)
        change_from_suffix = ""
        change_to_suffix   = self.chg_suffix
#         logger.debug("find = %s, replace = %s" % (find, replace))

        for cur_file_name in cur_file_names:
            file_name = os.path.join(cur_dirname, cur_file_name)
            lines = self.utils.read_file(file_name + change_from_suffix)
            with open(file_name + change_to_suffix, "w") as target:
                for line in lines:
                    if line.startswith(">"):
                        line1 = regex1.sub(replace1, line)
                    else:
                        line1 = line.upper()
#                     logger.debug(line1)
                    target.write(line1)


    def illumina_size_to_freq_in_chimer(self):
        find1           = ";size="
        replace1        = "frequency:"
        regex1          = re.compile(r"%s" % find1)

        cur_file_names = self.get_chimera_file_names(self.outdir)

        for file_chim in cur_file_names:
            file_chim_path = os.path.join(self.outdir, file_chim)
            lines = self.utils.read_file(file_chim_path)
            with open(file_chim_path, "w") as target:
                for line in lines:
                    line1 = regex1.sub(replace1, line)
                    target.write(line1)

    def illumina_rm_size_files(self):
        for idx_key in self.input_file_names:
            file_name = os.path.join(self.indir, self.input_file_names[idx_key] + self.chg_suffix)
            if os.path.exists(file_name):
                pass
                # os.remove(file_name)

    def check_if_chimera_dir_empty(self):
        if not os.listdir(self.outdir):
            self.utils.print_both('ERROR: Did not find files in %s, something is wrong. First check if you ran the command on a cluster. Exiting.\n' % self.outdir)
            sys.exit()

    def check_if_cluster_is_done(self, time_before):
        cluster_done = False
        check_qstat_cmd_line = "qstat | grep \"%s\" | grep chimera_ch | wc -l" % time_before
#         check_qstat_cmd_line = "qstat | grep vsearch"

        self.utils.print_both("check_qstat_cmd_line = %s" % check_qstat_cmd_line)

        try:
            p = subprocess.Popen(check_qstat_cmd_line, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            num_proc = int(output)
            self.utils.print_both("qstat is running %s 'vsearch' processes" % num_proc)
    #         pprint(p)

            if (num_proc == 0):
                cluster_done = True
    #         logger.debug("cluster_done from check_if_cluster_is_done = %s" % cluster_done)
        except:
            self.utils.print_both("Chimera checking can be done only on a cluster.")
            raise

        return cluster_done

    def create_chimera_cmd(self, ref_db):
        """
        /usr/local/bin/vsearch
        -uchime_denovo
        /Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/miseq/20150223/lane_1_B/analysis/reads_overlap/TGACCA_NNNNCGACG_1_MERGED-MAX-MISMATCH-3.unique.chg
        -uchimeout
        /Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/miseq/20150223/lane_1_B/analysis/chimera/TGACCA_NNNNCGACG_1_MERGED-MAX-MISMATCH-3.unique.chimeras.txt
        -chimeras
        /Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/miseq/20150223/lane_1_B/analysis/chimera/TGACCA_NNNNCGACG_1_MERGED-MAX-MISMATCH-3.unique.chimeras.txt.chimeric.fa
        -notrunclabels
        ---
        /usr/local/bin/vsearch
        -uchime_ref
        /Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/miseq/20150223/lane_1_B/analysis/reads_overlap/TGACCA_NNNNCGACG_1_MERGED-MAX-MISMATCH-3.unique.chg
        -uchimeout
        /Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/miseq/20150223/lane_1_B/analysis/chimera/TGACCA_NNNNCGACG_1_MERGED-MAX-MISMATCH-3.unique.chimeras.db
        -chimeras
        /Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test/miseq/20150223/lane_1_B/analysis/chimera/TGACCA_NNNNCGACG_1_MERGED-MAX-MISMATCH-3.unique.chimeras.db.chimeric.fa
        -notrunclabels
        -strand
        plus
        -db
        /groups/g454/blastdbs/rRNA16S.gold.fasta

        """
        command_line = []

        ref_or_novo_options = {self.denovo_suffix: "-uchime_denovo", self.ref_suffix: "-uchime_ref"}
        for suff, opt in ref_or_novo_options.items():
            input_file_name  = self.indir  + "/$filename_base" + self.chg_suffix
            output_file_name = self.outdir + "/$filename_base" + self.chimeras_suffix + suff

            ref_add = ""
            if (opt == "-uchime_ref"):
                ref_add = "-strand plus -db %s" % ref_db

            uchime_cmd = """%s %s %s -uchimeout %s -chimeras %s%s -notrunclabels %s
            """ % (self.usearch_cmd, opt, input_file_name, output_file_name, output_file_name, self.chimeric_suffix, ref_add)
            logger.debug("UUU = uchime_cmd = %s" % uchime_cmd)
            logger.debug("+++")
            command_line.append(uchime_cmd)

        return command_line


    def create_chimera_cmd_old(self, input_file_name, output_file_name, ref_or_novo, ref_db = ""):
        """
        http://www.drive5.com/usearch/manual/uchime_denovo.html
        from usearch -help
        Chimera detection (UCHIME ref. db. mode):
          usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
            [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]

        Chimera detection (UCHIME de novo mode):
          usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
             [-uchimeout results.uch] [-uchimealns results.alns]
          Input is estimated amplicons with integer abundances specified using ";size=N".
        usearch -uchime_denovo amplicons.fasta -uchimeout results.uchime
        """

        uchime_cmd_append = ""
        db_cmd_append     = ""
        dir_cmd_append    = ""

        if (ref_or_novo == "denovo"):
            uchime_cmd_append = " -uchime_denovo "
            output_file_name  = output_file_name + self.chimeras_suffix + self.denovo_suffix
        elif (ref_or_novo == "ref"):
            uchime_cmd_append = " -uchime_ref "
            output_file_name  = output_file_name + self.chimeras_suffix + self.ref_suffix
            db_cmd_append     = " -db " + ref_db
            dir_cmd_append    = " -strand plus"
        else:
            self.utils.print_both("Error: Incorrect method, should be \"denovo\" or \"ref\"")
        self.utils.print_both("output_file_name = %s" % output_file_name)


        uchime_cmd = C.clusterize_cmd
        if self.utils.is_local():
            uchime_cmd = ""
        uchime_cmd += " "
        uchime_cmd += self.usearch_cmd
        logger.debug("self.usearch_cmd FROM create_chimera_cmd = %s" % (uchime_cmd))

        uchime_cmd += uchime_cmd_append + input_file_name
        logger.debug("uchime_cmd_append FROM create_chimera_cmd = %s" % (uchime_cmd_append))


        uchime_cmd += db_cmd_append

        logger.debug("db_cmd_append FROM create_chimera_cmd = %s" % (db_cmd_append))

        uchime_cmd += " -uchimeout " + output_file_name

        """if we need nonchimeric for denovo and db separate we might create them here"""
        uchime_cmd += " -nonchimeras "
        uchime_cmd += (output_file_name + self.nonchimeric_suffix)

        uchime_cmd += " -chimeras " + (output_file_name + self.chimeric_suffix)
        uchime_cmd += dir_cmd_append

        uchime_cmd += " -notrunclabels"


        logger.debug("uchime_cmd FROM create_chimera_cmd = %s" % (uchime_cmd))
        return uchime_cmd

    def get_sge_cluster_name(self):
      # import subprocess
      result = subprocess.run(['qstat', '-F'], stdout = subprocess.PIPE)
      a1 = result.stdout.decode('utf-8').split()
      for line in a1:
          if (line.find("hostname") != -1): #qf:hostname=grendel-01.bpcservers.private
                  return line.split("=")[1].split("-")[0]

    def get_sge_slot_number(self): # doesn't work on cricket because: 	qc:slots=12 and qc:slots=8
      result = subprocess.run(['qstat', '-F', 'slots'], stdout = subprocess.PIPE)
      a1 = result.stdout.decode('utf-8').split()
      slots = []
      for line in a1:
         if line.startswith('qc:slots'):
             slots.append(int(line.split("=")[-1]))
      slots_uniq = set(slots)
      return max(slots_uniq)

    # TODO: temp! take from util. change illumina-files to use util, too
    #   create_job_array_script(self, command_line, dir_to_run, files_list, runobj)
    # feb 25 2019 removed, because didn't work on grendel:
    #  Use the allslots pe and all available slots on that cluster
    # #$ -pe allslots %s
    def create_job_array_script(self, script_file_name_base, command_line, dir_to_run, files_list):
        # sge_slot_number = self.get_sge_slot_number()
        sge_cluster_name = self.get_sge_cluster_name()
        sge_slot_number  = self.cluster_slots[sge_cluster_name][0]
        logger.debug("sge_slot_number FROM create_job_array_script = %s" % (sge_slot_number))

        files_string    = " ".join(files_list)
        files_list_size = len(files_list)
#         command_file_name = os.path.basename(command_line.split(" ")[0])
        script_file_name  = script_file_name_base + "_" + self.runobj.run + "_" + self.runobj.lane_name + ".sh"
        script_file_name_full = os.path.join(dir_to_run, script_file_name)
        log_file_name     = script_file_name + ".sge_script.sh.log"
        email_mbl         = self.utils.make_users_email()
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
# max_running_tasks
#$ -tc 15
-# Use the allslots pe and all available slots on that cluster
#$ -pe allslots %s
#$ -t 1-%s
# Now the script will iterate %s times.

  file_list=(%s)

  i=$(expr $SGE_TASK_ID - 1)
  echo "i = $i"
  . /bioware/root/Modules/etc/profile.modules
  module load bioware
  module load vsearch

  INFILE=${file_list[$i]}

  filename=$(basename $INFILE)
  echo "INFILE = $INFILE"
  filename_base="${filename%%.*}"
  echo "filename_base = $filename_base"
  echo "%s"
  echo "%s"
  %s
  %s
''' % (script_file_name, log_file_name, email_mbl, sge_slot_number, files_list_size, files_list_size, files_string, command_line[0], command_line[1], command_line[0], command_line[1])
# ''' % (script_file_name, log_file_name, email_mbl, files_list_size, files_list_size, files_string, command_line)
                )
        self.utils.open_write_close(script_file_name_full, text)
        return script_file_name

    def create_not_SGE_script(self, script_file_name_base, command_line, dir_to_run, files_list):

        files_string = " ".join(files_list)
        script_file_name = script_file_name_base + "_" + self.runobj.run + "_" + self.runobj.lane_name + ".sh"
        script_file_name_full = os.path.join(dir_to_run, script_file_name)
        text = (
                '''#!/bin/bash

    file_list=(%s)

    . /bioware/root/Modules/etc/profile.modules
    module load bioware
    module load vsearch

    n=0
    for INFILE in "${file_list[@]}"
    do      
    n=$[n + 1]
    echo $n
    echo "INFILE = $INFILE"
    filename=$(basename $INFILE)
    filename_base="${filename%.*}"
    echo "filename_base = $filename_base"

    echo "%s"
    echo "%s"
    %s
    %s
    done
    ''' % (files_string, command_line[0], command_line[1], command_line[0], command_line[1])
            # ''' % (script_file_name, log_file_name, email_mbl, files_list_size, files_list_size, files_string, command_line)
        )
        self.utils.open_write_close(script_file_name_full, text)
        return script_file_name

    def chimera_checking(self):
        chimera_region_found = False

        file_list   = self.dirs.get_all_files_by_ext(self.indir, self.chg_suffix)
        logger.debug("FFF = file_list = %s" % (file_list))

#         TODO: method
        dna_region = list(set([self.runobj.samples[idx_key].dna_region for idx_key in self.input_file_names]))[0]
        if dna_region in C.regions_to_chimera_check:
            chimera_region_found = True
        else:
            logger.debug('region not checked: ' +  dna_region)
        ref_db = self.get_ref_db(dna_region)
        command_line = self.create_chimera_cmd(ref_db)
        sh_script_file_name = self.create_job_array_script("chimera_checking", command_line, self.indir, file_list)
        script_file_name_full = os.path.join(self.indir, sh_script_file_name)
        self.utils.call_sh_script(script_file_name_full, self.indir)
        self.utils.print_both("self.dirs.chmod_all(%s)" % (self.indir))
        self.dirs.chmod_all(self.indir)
        logger.debug('sh_script_file_name: ' +  sh_script_file_name)
        if not chimera_region_found:
            return ('NOREGION', 'No regions found that need checking', '')
        else:
            return ("The vsearch commands were created")

    def get_chimeric_ids(self):
        ids = set()
        chimera_file_names = self.get_chimera_file_names(self.outdir)
        file_ratio = self.check_chimeric_stats()

        for file_name in chimera_file_names:
#             logger.debug("from get_chimeric_ids: file_name = %s" % file_name)
            if file_name.endswith(self.chimeric_suffix):
                both_or_denovo = self.get_chimeras_suffix(file_ratio, file_name)
#                 TODO: run ones for each file_base = ".".join(file_name.split(".")[0:3]) (for txt and db)
                if file_name.endswith(both_or_denovo):
                    file_name_path = os.path.join(self.outdir, file_name)
                    self.utils.print_both("Get ids from %s" % file_name_path)
                    read_fasta     = fa.ReadFasta(file_name_path)
                    ids.update(set(read_fasta.ids))
        return ids


    def get_chimeras_suffix(self, file_ratio, file_name):
        """ use only de-novo (.txt) chimeric if
            check_chimeric_stats shows
            ratio ref to de-novo > 3
            e.g.
            if denovo_only:
                chimeric_suffix = self.chimeras_suffix + self.denovo_suffix + self.chimeric_suffix
            if no:
                chimeras_suffix = self.chimeric_suffix

            if file_name.endswith(chimeric_suffix):
            ...
                #     first_name, last_name = get_name()

        """
#         for file_basename in file_ratio:
        (percent_ref, ratio) = file_ratio[".".join(file_name.split(".")[0:3])]

        chimeric_fa_suffix = ""
#         logger.debug("percent_ref = %s, ratio = %s" % (percent_ref, ratio))
#         if (percent_ref > 15) and (ratio > 2):
        if ratio > 3:
            chimeric_fa_suffix = self.chimeras_suffix + self.denovo_suffix + self.chimeric_suffix
        else:
            chimeric_fa_suffix = self.chimeric_suffix
        return chimeric_fa_suffix

    def move_out_chimeric(self):
        chimeric_ids = self.get_chimeric_ids()
        for idx_key in self.input_file_names:
            fasta_file_path    = os.path.join(self.indir, self.input_file_names[idx_key])
            read_fasta         = fa.ReadFasta(fasta_file_path)
            read_fasta.close()

            non_chimeric_file  = fasta_file_path + self.nonchimeric_suffix
            non_chimeric_fasta = fa.FastaOutput(non_chimeric_file)

            fasta              = fa.SequenceSource(fasta_file_path, lazy_init = False)
            while fasta.next():
                if not fasta.id in chimeric_ids:
                    non_chimeric_fasta.store(fasta, store_frequencies = False)
            non_chimeric_fasta.close()


    def check_chimeric_stats(self):
        all_lines_suffix      = self.denovo_suffix # ".txt" or ".db, doesn't matter"
        chimera_ref_suffix    = self.ref_suffix + self.chimeric_suffix #".db.chimeric.fa"
        chimera_denovo_suffix = self.denovo_suffix + self.chimeric_suffix # ".txt.chimeric.fa"
        filenames             = self.get_basenames(self.get_current_filenames(self.outdir))
        file_ratio            = {}
        for file_basename in filenames:
            # logger.debug(file_basename)
            all_lines      = 0
            ref_lines      = 0
            denovo_lines   = 0
            ratio          = 0
            percent_ref    = 0
            percent_denovo = 0

            all_lines_file_name    = os.path.join(self.outdir, file_basename + all_lines_suffix)
            ref_lines_file_name    = os.path.join(self.outdir, file_basename + chimera_ref_suffix)
            denovo_lines_file_name = os.path.join(self.outdir, file_basename + chimera_denovo_suffix)

            all_lines    = int(self.wccount(all_lines_file_name) or 0)
            ref_lines    = int(self.get_fa_lines_count(ref_lines_file_name) or 0)
            denovo_lines = int(self.get_fa_lines_count(denovo_lines_file_name) or 0)

            # denovo_lines = int(denovo_lines or 0)
            if (ref_lines == 0) or (all_lines == 0):
                file_ratio[file_basename] = (0, 0)
                continue
            else:
                percent_ref = self.percent_count(all_lines, ref_lines)

            if (denovo_lines == 0):
                file_ratio[file_basename] = (percent_ref, percent_ref) #use ref instead of ratio, because we are actually looking for a huge difference between ref and denovo (ref > 15 and denovo = 0)
                continue

            if (denovo_lines > 0):
                ratio          = self.count_ratio(ref_lines, denovo_lines)
                percent_denovo = self.percent_count(all_lines, denovo_lines)
            file_ratio[file_basename] = (percent_ref, ratio)
            # percent_ref = int(percent_ref or 0)
            if (percent_ref > 15):
                self.utils.print_both("=" * 50)

                self.utils.print_both(file_basename)
                # logger.debug("all_lines_file_name = %s, ref_lines_file_name = %s, denovo_lines_file_name = %s" % (all_lines_file_name, ref_lines_file_name, denovo_lines_file_name))
                self.utils.print_both("all_lines = %s, ref_lines = %s, denovo_lines = %s" % (all_lines, ref_lines, denovo_lines))
                self.utils.print_both("ratio = %s" % ratio)
                self.utils.print_both("percent_ref = %s, percent_denovo = %s" % (percent_ref, percent_denovo))
        return file_ratio


    def get_basenames(self, filenames):
        file_basenames = set()
        for f in filenames:
            file_basename = ".".join(f.split(".")[0:3])
            if file_basename.endswith(self.base_suffix):
                file_basenames.add(file_basename)

        return file_basenames

    def wccount(self, filename):
        return subprocess.check_output(['wc', '-l', filename]).split()[0]

    def count_ratio(self, ref_num, denovo_num):
        try:
            return float(ref_num or 0) / float(denovo_num or 0)
        except ZeroDivisionError:
            # logger.debug("There is no denovo chimeras to count ratio.")
            pass

    def get_fa_lines_count(self, file_name):
        # todo: use fastalib to get cnt?
        # return fa.SequenceSource(file_name, lazy_init = False).total_seq
        try:
            file_open = open(file_name)
            return len([l for l in file_open.readlines() if l.startswith('>')])
        except IOError:
            e = sys.exc_info()[1]
            self.utils.print_both(e)
            return 0
            # logger.error("%s\nThere is no such file: %s" % (e, file_name))

    def percent_count(self, all_lines, chimeric_count):
        try:
            return float(chimeric_count or 0) * 100 / float(all_lines or 0)
        except ZeroDivisionError:
            # logger.error("There is no denovo chimeras to count ratio.")
            pass




    """
    -----------------------------------------------------------------------------
        For 454.
        not tested
    """
    def chimera_denovo(self):
        chimera_region_found = False
        output = {}
        cluster_id_list = []



        for idx_key in self.idx_keys:
            input_file_name  = os.path.join(self.indir,  idx_key +'.abund.fa')
            if os.path.isfile(input_file_name):
                output_file_name = os.path.join(self.outdir, idx_key +'.chimera.denovo')
                #open(output_file_name, 'a').close()  # make sure file exists
                log_file = os.path.join(self.outdir,idx_key+".denovo.log")

                dna_region       = self.runobj.samples[idx_key].dna_region
                logger.debug("dna_region = %s" % dna_region)
                if self.runobj.vamps_user_upload:
                    # VAMPS users can chimera check regardless of region chosen
                    chimera_region_found = True
                else:
                    if dna_region in C.regions_to_chimera_check:
                        chimera_region_found = True
                    else:
                        logger.debug('region not checked: ' +  dna_region)
                        continue


                self.utils.print_both("input_file_name = %s \noutput_file_name = %s" % (input_file_name, output_file_name))

    #             uchime_cmd = C.clusterize_cmd
    #             uchime_cmd += " "
    #             uchime_cmd += self.usearch_cmd
    #             uchime_cmd += " --uchime "
    #             uchime_cmd += input_file_name
    #             uchime_cmd += " --uchimeout "
    #             uchime_cmd += output_file_name
    #             uchime_cmd += " --abskew "
    #             uchime_cmd += self.abskew
                uchime_cmd=''
                if self.use_cluster:
                    uchime_cmd += C.clusterize_cmd
                    uchime_cmd += " "
                    uchime_cmd += " -log "
                    uchime_cmd += log_file
                    uchime_cmd += " "
                uchime_cmd += self.usearch_cmd
                uchime_cmd += " -uchime_denovo "
                uchime_cmd += input_file_name
                uchime_cmd += " -uchimeout "
                uchime_cmd += output_file_name

                logger.debug("uchime_denovo_cmd = %s" % (uchime_cmd))

                try:
                    logger.info("chimera denovo command: " + str(uchime_cmd))
    #                 subprocess.Popen(uchime_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                    self.utils.print_both("chimera denovo command: " + str(uchime_cmd))
                    #output[idx_key] = subprocess.Popen(uchime_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    output[idx_key] = subprocess.check_output(uchime_cmd, shell=True)
                    self.utils.print_both("chimera denovo result: " + str(output[idx_key]))
                    #self.utils.print_both("output[idx_key] = %s" % output[idx_key])
                    #if idx_key in output and len(output[idx_key].split()) > 1:
                    #self.utils.print_both(output[idx_key].split()[2])
                    items = output[idx_key].split()
                    if len(items) > 2:
                        cluster_id_list.append(items[2])




                except OSError:
                    e = sys.exc_info()[1]
                    self.utils.print_both("Error: Problems with this command: %s" % (uchime_cmd))
                    if self.utils.is_local():
                        print >>sys.stderr, "Error: Execution of %s failed: %s" % (uchime_cmd, e)
                    else:
                        print >>sys.stderr, "Error: Execution of %s failed: %s" % (uchime_cmd, e)
                        self.utils.print_both("Error: Execution of %s failed: %s" % (uchime_cmd, e))
                        raise

# ???
        if not chimera_region_found:
            return ('NOREGION', 'No regions found that need checking', '')

        # ???
#         for idx_key in output:
#             if len(output[idx_key]) > 50 or len(output[idx_key]) < 40:
#                 return ('ERROR','uchime ref may have broken or empty', idx_key)

        # finally
        self.utils.print_both('Finished Chimera Denovo')
        if cluster_id_list:
            return ('SUCCESS', 'uchime ref seems to have been submitted successfully', cluster_id_list)
        else:
            return ('ERROR', 'uchime ref returned no cluster IDs', cluster_id_list)

    def chimera_reference(self):

        chimera_region_found = False
        output = {}
        cluster_id_list = []
        for idx_key in self.run_keys:

            dna_region  = self.runobj.samples[idx_key].dna_region
            if self.runobj.vamps_user_upload:
                # VAMPS users can chimera check regardless of region chosen
                chimera_region_found = True
            else:
                if dna_region in C.regions_to_chimera_check:
                    chimera_region_found = True
                else:
                    logger.debug('region not checked: ' + dna_region)
                    continue


            input_file_name  = os.path.join(self.indir,  idx_key +'.abund.fa')
            output_file_name    = os.path.join(self.outdir,idx_key+".chimera.ref")
            #open(output_file_name, 'a').close()  # make sure file exists
            log_file = os.path.join(self.outdir,idx_key+".ref.log")
            logger.debug("OUT FILE NAME: " + output_file_name)

            #out_file_name = self.prefix[idx_key] + ".chimeras.db"
            input_file_name  = os.path.join(self.indir,  idx_key +'.abund.fa')
            if os.path.isfile(input_file_name):
                output_file_name    = os.path.join(self.outdir,idx_key+".chimera.ref")
                #open(output_file_name, 'a').close()  # make sure file exists
                log_file = os.path.join(self.outdir,idx_key+".ref.log")
                logger.debug("OUT FILE NAME: " + output_file_name)
                # which ref db to use?
                ref_db = ''
                if dna_region.upper() == 'ITS':
                    logger.debug("got an ITS dna region so using refdb: " + self.its_refdb)
                    ref_db = self.its_refdb
                else:
                    logger.debug("using standard refdb: " + self.refdb)
                    ref_db = self.refdb

                uchime_cmd=''
                if self.use_cluster:
                    uchime_cmd = C.clusterize_cmd
                    uchime_cmd += " "
                    uchime_cmd += " -log "
                    uchime_cmd += log_file
                    uchime_cmd += " "
                uchime_cmd += self.usearch_cmd
                uchime_cmd += " -uchime_ref "
                uchime_cmd += input_file_name
                uchime_cmd += " -uchimeout "
                uchime_cmd += output_file_name
                uchime_cmd += " -db "
                uchime_cmd += ref_db
                uchime_cmd += " -strand "
                uchime_cmd += "plus"

                logger.debug("uchime_ref_cmd = %s" % (uchime_cmd))

                try:

                    logger.info("vsearch version: " % (self.utils.get_vsearch_version))
                    logger.info("chimera reference command: " + str(uchime_cmd))
                    output[idx_key] = subprocess.check_output(uchime_cmd, shell=True)
                    #logger.debug('outsplit',output[idx_key].split()[2])
                    cluster_id_list.append(output[idx_key].split()[2])
                    #logger.debug('Have %d bytes in output' % len(output))
                    #logger.debug('ref',idx_key,output,len(output))
                    if len(output[idx_key]) < 50 and len(output[idx_key]) > 40:
                        logger.debug(idx_key + " uchime ref seems to have been submitted successfully")
                    else:
                        if self.use_cluster:
                            print >>sys.stderr, "Error: uchime ref may be broke"
                            self.utils.print_both("Error: uchime ref may be broke")

                except OSError:
                    e = sys.exc_info()[1]
                    print >>sys.stderr, "Error: Execution of chimera_reference failed: %s" % (uchime_cmd, e)
                    self.utils.print_both("Error: Execution of chimera_reference failed: %s" % (uchime_cmd, e))
                    raise

        if not chimera_region_found:
            return ('NOREGION','No regions found that need checking','')

        for idx_key in output:
            if (len(output[idx_key]) > 50 or len(output[idx_key]) < 40) and self.use_cluster:
                return ('ERROR','uchime ref may have broken or empty',idx_key)
        self.utils.print_both('Finished Chimera Reference')
        return ('SUCCESS','uchime ref seems to have been submitted successfully',cluster_id_list)


    def write_chimeras_to_deleted_file(self):

        for idx_key in self.run_keys:
            # open  deleted file and append chimera to it
            # open and read both chimeras files: chimeras.db and chimeras.txt

            # hash to remove dupes
            chimera_deleted = {}
            denovo_file = os.path.join(self.outdir, idx_key +'.chimera.denovo')
            ref_file = os.path.join(self.outdir,idx_key+".chimera.ref")
            # deleted file is in trimming dir for vampsuser
            deleted_file = os.path.join(self.indir, idx_key+".deleted.txt")
            for file in [denovo_file, ref_file]:
                if os.path.isfile(file):
                    fh = open(file,"r")
                    # make a list of chimera deleted read_ids
                    for line in fh.readlines():
                        lst = line.strip().split()
                        id = lst[1].split(';')[0]
                        chimera_yesno = lst[-1]
                        if(chimera_yesno) == 'Y':
                            chimera_deleted[id] = 'chimera'
            # open to append as trimming deletions are already there
            fh_del = open(deleted_file,"a")
            for id in chimera_deleted:
                fh_del.write(id+"\tChimera\n")
            fh_del.close()

# # http://drive5.com/uchime/uchime_quickref.pdf
# # The --uchimeout file is a tab-separated file with the following 17 fields.
# # Field Name Description
# # 1 Score Value >= 0.0, high score means more likely to be a chimera.
# # 2 Query Sequence label
# # 3 Parent A Sequence label
# # 4 Parent B Sequence label
# # 5 IdQM %id between query and model made from (A, crossover, B)
# # 6 IdQA %id between query and parent A.
# # 7 IdQB %id between query and parent B
# # 8 IdAB %id between parents (A and B).
# # 9 IdQT %id between query and closest reference sequence / candidate parent.
# # 10 LY Yes votes on left
# # 11 LN No votes on left
# # 12 LA Abstain votes on left
# # 13 RY Yes votes on right
# # 14 RN No votes on right
# # 15 RA Abstain votes on right
# # 16 Div Divergence ratio, i.e. IdQM - IdQT
# # 17 YN Y (yes) or N (no) classification as a chimera. Set to Y if score >= threshold
