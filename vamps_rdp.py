#!/usr/bin/env python

##!/bioware/python/bin/python



##!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

import os
from stat import * # ST_SIZE etc
import sys
import shutil
import types
import time
import random
import csv
from time import sleep
import argparse
sys.path.append( '/bioware/python/lib/python2.7/site-packages/' )

import datetime
import subprocess
import logging
import ConMySQL



   
def run_fasta2tax(myobject):
    """
    
    """   
    user = myobject['user']
    project = myobject['project']
    run = myobject['runcode']
    site = myobject['site']
    output_dir = myobject['output_dir']
    
    outSeqTable = 'vamps_upload_trimseq'
    
    
    
    if site == 'vamps':

        # for vamps:
        db_hostname_vamps = "vampsdb";
        db_hostname_user = "bpcdb2";
        dbName_vamps = 'vamps';
        dbName_user = 'vamps_user_uploads';
        web_user = "vampshttpd";
        cluster_path= "/xraid2-2/vampsweb/vamps";
        
   
    else:

        # for vampsdev
        db_hostname_vamps = "vampsdev.mbl.edu"
        db_hostname_user = "bpcdb2"
        dbName_vamps = 'vamps'
        dbName_user = 'vampsdev_user_uploads'
        web_user = "vampsdevhttpd"
        cluster_path= "/xraid2-2/vampsweb/vampsdev"
    
    
    
    
    if output_dir and os.path.exists(output_dir):
        print "files path exists:",output_dir
        rundir = output_dir
        #gast_input_source = 'files'
        #file_base = output_dir
        # This may be a mobedac upload and we should try to use the files here
        # rather than look to the database for data

    else:
        runcode = user+'_'+run+'_rdp'
        output_dir = os.path.join(cluster_path,'tmp',runcode)
        rundir = output_dir 
        print "Files path doesn't exist: attempting to get data from database"
        print "Creating directory",output_dir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
    
        
    print db_hostname_vamps, dbName_vamps, cluster_path
    obj=ConMySQL.New(db_hostname_vamps, dbName_vamps, cluster_path)
    vampsconn = obj.get_conn()
    vamps_cursor = vampsconn.cursor()
    obj=ConMySQL.New(db_hostname_user, dbName_user, cluster_path)
    userconn = obj.get_conn()
    user_cursor = userconn.cursor()
    
    
    #myconn_user = MyConnection(host=db_hostname_user,db=dbName_user)
    #myconn_vamps = MyConnection(host=db_hostname_vamps,db=dbName_vamps)
    dbconf = cluster_path+"/.dbconf"
    

    if os.path.exists(dbconf):
        res = subprocess.check_output("/bin/cat "+dbconf, shell=True)
        db_user = res.strip().split()[0]
        db_password = res.strip().split()[1]
    else:
        sys.exit("Unable to connect to the database, please contact your database administrator for access privileges")

    
    apps_directory = os.path.join(cluster_path,'apps')
    got_atleast_one = 0;
    
    q0 = "select distinct dataset from vamps_projects_datasets_pipe where project='"+project+"' and dataset != '' and dataset != 'NoKey' "
    
    #dsrows = myconn_vamps.execute_fetch_select(q0)
    vamps_cursor.execute(q0)
    dsrows = vamps_cursor.fetchall()
    if not dsrows:
        print "No datasets found using query:", q0
        sys.exit()
    for ds in dsrows:                
        ds  = ds[0]
        
        where = "WHERE deleted = 0 AND user = '"+user+"' AND project = '"+project+"' AND dataset='"+ds+"' "
        ffilename = os.path.join(output_dir,project+'--'+ds+".fa")
    
    
        # RUN IN FOREGROUND (wait until return to start fasta2tax)
        db2fasta_cmd = cluster_path + "/db2fasta_vamps -id read_id -seq sequence -sql \"SELECT sequence, read_id FROM "+outSeqTable+" "+where+" \" -o "+ffilename+" -db user -site "+site+" -v"     
        print  db2fasta_cmd
    
        db2fasta_result = subprocess.call(db2fasta_cmd, shell=True)
    
        #print "Sleeping for 20\n";
        sleep(2);
           
           
        if not os.path.exists(ffilename):
            print  "Unable to locate input fasta file: "+ffilename
            sys.exit()
        else:
            pass
                      
            # run fasta2tax.pl on $ffilename
            fasta2tax_cmd = apps_directory+"/fasta2tax.pl "
            fasta2tax_cmd += " --user="+user
            fasta2tax_cmd += " --inputfile="+ffilename
            fasta2tax_cmd += " --project="+project
            fasta2tax_cmd += " --dataset="+ds
            fasta2tax_cmd += " --path-to-apps="+apps_directory
            fasta2tax_cmd += " --database="+dbName_vamps
            fasta2tax_cmd += " --table1=vamps_data_cube_uploads";
            fasta2tax_cmd += " --table2=vamps_junk_data_cube_pipe";
            fasta2tax_cmd += " --db_hostname="+db_hostname_vamps
            fasta2tax_cmd += " --db-user="+db_user
            fasta2tax_cmd += " --db-password="+db_password
            #fasta2tax_cmd += " & ";
               
            print  fasta2tax_cmd
            fasta2tax_result = subprocess.call(fasta2tax_cmd, shell=True)

            print "updating vamps_projects_datasets_pipe"
            updatePDP = "UPDATE vamps_projects_datasets_pipe set has_tax='1' where project='"+project+"' and dataset='"+ds+"'";
            #myconn_vamps.execute_no_fetch(updatePDP)
            #myconn_user.execute_no_fetch(updatePDP)
            vamps_cursor.execute(updatePDP)
            user_cursor.execute(updatePDP)
            
 
    
    #######################################
    #
    # insert ignore vamps_upload_info
    #
    #######################################
    print "updating vamps_upload_info"
    updateInfo = "UPDATE vamps_upload_info set has_tax='1' where user='"+user+"' and project_name='"+project+"' ";
    #myconn_vamps.execute_no_fetch(updateInfo)
    #myconn_user.execute_no_fetch(updateInfo)
    vamps_cursor.execute(updateInfo)
    user_cursor.execute(updateInfo)   
    
    #######################################
    #
    # insert ignore vamps_users
    #
    #######################################
    print "updating vamps_users"
    insert_projects = "INSERT IGNORE INTO vamps_users (user, project) values ('"+user+"', '"+project+"' )";
    #myconn_vamps.execute_no_fetch(insert_projects)
    vamps_cursor.execute(insert_projects)   
    
    vampsconn.commit()
    userconn.commit()
    vamps_cursor.close()
    user_cursor.close()
    userconn.close()
    vampsconn.close()

        
if __name__ == '__main__':
    
    
    # DEFAULTS
    site = 'vampsdev'
    
    

    data_object = {}
    

    
    myusage = """usage: vamps_rdp.py  [options]
         
         
         
         where
            
            -site       vamps or vampsdev        
            
            -r,   code
            
            -u, --user       Needed for otu naming and db access.
            
            -project
            
            
    
    
    """
    parser = argparse.ArgumentParser(description="" ,usage=myusage)                 
    
    
                                                     
    parser.add_argument("-site",                    required=True,  action="store",   dest = "site", 
                                                        help="""database hostname: vamps or vampsdev
                                                        [default: vampsdev]""")  
    parser.add_argument("-r", "--runcode",          required=True,  action="store",   dest = "runcode", 
                                                        help="")  
    parser.add_argument("-u", "--user",             required=True,  action="store",   dest = "user", 
                                                        help="user name")         
    parser.add_argument("-p", "--project",          required=True,  action='store', dest = "project", 
                                                        help="")   
    parser.add_argument("-out", "--output_dir",     required=False,  action="store",   dest = "output_dir", default='',
                                                        help = '')                                         
    
    #parser.add_argument("-cl", "--use_cluster",     required=False,  action="store",   dest = "use_cluster", default=True,
    #                                                    help = '')                                                
                                        
                                                
                                                
    args = parser.parse_args()
    
    
    
    # fill command line object
    data_object['datetime']     = str(datetime.date.today())
    data_object['runcode']      = args.runcode
    data_object['site']         = args.site
    data_object['user']         = args.user
    data_object['output_dir']   = args.output_dir
    data_object['project']      = args.project[:1].capitalize() + args.project[1:]
    
    
    
#    data_object['use_cluster']      = args.use_cluster
#    if data_object['use_cluster'] == 'True' or data_object['use_cluster'] ==  'true':
#        data_object['use_cluster'] = True
#    else:
#        data_object['use_cluster'] = False
#     if data_object['site'] == 'vamps':
#         db_host_vamps   = 'vampsdb'
#         db_host_user    = 'bpcdb2'
#         db_name_vamps   = 'vamps'
#         db_name_user    = 'vamps_user_uploads'
#         db_home         = '/xraid2-2/vampsweb/vamps/'
#     else:
#         db_host_vamps   = 'vampsdev'
#         db_host_user    = 'bpcdb2'
#         db_name_vamps   = 'vamps'
#         db_name_user    = 'vampsdev_user_uploads'
#         db_home         = '/xraid2-2/vampsweb/vampsdev/'
#     
#     
#     vamps_obj=ConMySQL.New(db_host_vamps, db_name_vamps, db_home)
#     data_object['vamps_cursor'] = vamps_obj.get_cursor()
#     user_obj=ConMySQL.New(db_host_user, db_name_user, db_home)
#     data_object['user_cursor'] = user_obj.get_cursor()
    
    run_fasta2tax(data_object)
        
    