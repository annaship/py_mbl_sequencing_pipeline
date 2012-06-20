import subprocess
import sys, os,stat
import time
import shutil
from pipeline.pipelinelogging import logger
import constants as C

class Vamps:
    """Uploads date to the VAMPS (or vampsdev) database"""
    Name = "VAMPS"
    def __init__(self, run = None):

        self.run 	 = run
        self.outdir  = run.output_dir
        self.basedir = run.basedir
        self.rundate = self.run.run_date
        self.use_cluster = 1
        self.project = run.project
        self.dataset = run.dataset
        
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
        self.refdb_dir = '/xraid2-2/vampsweb/blastdbs/'
        # 1) clustergast
        # 2) gast cleanup
        # 3) gast2tax
        print "vampsupload 1- init"
    

    def taxonomy(self, lane_keys):
        """
        fill vamps_data_cube, vamps_junk_data_cube and vamps_taxonomy tables
        """
        print "vampsupload 3- tax"
        pass
        
    def sequences(self, lane_keys):
        """
        fill vamps_sequences table
        """
        print "vampsupload 4- seqs"
        pass  
        
    def exports(self, lane_keys):
        """
        fill vamps_exports table
        """
        print "vampsupload 5- exports"
        pass
        
    def projects(self, lane_keys):
        """
        fill vamps_projects_datasets table
        """
        print "vampsupload 2- projects"
        pass
        
    def info(self, lane_keys):
        """
        fill vamps_project_info table
        """
        print "vampsupload 6- info"
        pass

        
