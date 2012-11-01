#!/bioware/python/bin/python

import os
import MySQLdb

class New(object):

    def __init__(self,host,dbname,home):
        self.db_host = host
        self.db_name = dbname
        self.db_home = home
	    #my $db_info = shift;
	    #my $version = "mysql5";
        #print self.db_host, self.db_name, self.db_home
        
           
        userconf = os.popen("cat "+home+"/.dbconf").readlines()
        
        
        self.db_user = userconf[0].strip()
        self.db_pass = userconf[1].strip()
        
        self.db = MySQLdb.connect(host=self.db_host, db=self.db_name, user=self.db_user, passwd=self.db_pass )
        self.cursor = self.db.cursor()
    
    def get_conn(self):
        return self.db
        
    def get_cursor(self):
        return self.cursor
        
    def get_db_user(self):
        return self.db_user 
        
    def get_db_host(self):
        return self.db_host
        
    def get_db_name(self):
        return self.db_name 





