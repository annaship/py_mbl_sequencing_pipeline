import sys
import fastalib as u
import MySQLdb
from os import listdir, walk
from os.path import isfile, join

from pipeline.pipelinelogging import logger
#import logging
import constants as C

class MyConnection:
  """
  Connection to env454
  By default takes parameters from "conf.txt", host = "newbpcdb2"
  if different: change conf.txt and use my_conn = MyConnection(config_file_name, server_name)
  Conf.txt has lines as: bpcdb2:bpcdb2:3306:my_password                                                                                                            
  """
  def __init__(self, file_name="db_conn.conf", server_name="newbpcdb2_ill"):
    self.file_name   = file_name
    self.server_name = server_name
    self.conn        = None
    self.cursor      = None
    
    try:
      content = [line.strip() for line in open(self.file_name).readlines()]
      for line in content:
        fields = line.split(':')
        if fields[0] == self.server_name:
          print "server_name = " + str(self.server_name)
          print "=" * 40
          # conn = MySQLdb.connect (host = str(fields[1]), port = int(fields[2]), user = "ashipunova", passwd = str(fields[3]), db = "env454")
          # self.conn = MySQLdb.connect (host = str(fields[1]), port = int(fields[2]), user = "ashipunova", passwd = str(fields[3]), db = "env454")
          self.conn = MySQLdb.connect (host = str(fields[1]), port = int(fields[2]), user = "ashipunova", passwd = str(fields[3]), db = str(fields[4]))
          self.cursor = self.conn.cursor()       
    except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      raise
    except:                       # catch everything
      print "Unexpected:"         # handle unexpected exceptions
      print sys.exc_info()[0]     # info about curr exception (type,value,traceback)
      raise                       # re-throw caught exception   


class dbUpload:
    """db upload methods"""
    Name = "dbUpload"
    def __init__(self, run = None):

        self.run 	 = run
        self.outdir  = run.output_dir
        try:
            self.basedir = run.basedir
        except:
            self.basedir = self.outdir
        self.rundate = self.run.run_date
        self.use_cluster = 1
        self.fasta_dir        = self.run.input_dir + "/fasta/" 
        self.filenames   = []
#        get_fasta_file_names(fasta_file_path)
#        self.fasta       = u.SequenceSource(fasta_file_path) 

        
#        os.environ['SGE_ROOT']='/usr/local/sge'
#        os.environ['SGE_CELL']='grendel'
#        path = os.environ['PATH']
#        os.environ['PATH'] = path + ':/usr/local/sge/bin/lx24-amd64'
        #First step is to check for (or create via mothur)
        # a uniques fasta file and a names file 
        # one for each dataset.
        # If we are here from a vamps gast process
        # then there should be just one dataset to gast
        # but if MBL pipe then many datasets are prbably involved.
        self.refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
   
    def get_fasta_file_names(self, fasta_dir):
      for (dirpath, dirname, files) in walk(fasta_dir):
        return files
        
        
#    def 
#            logger.info("Finished clustergast")

     

