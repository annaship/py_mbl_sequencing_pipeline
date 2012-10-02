#!/usr/bin/env python

#!/bioware/python/bin/python

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
import random
import csv
from time import sleep
#sys.path.append( '/bioware/python/lib/python2.7/site-packages/' )
import MySQLdb
import ConMySQL
import datetime
from pipeline.pipelinelogging import logger

def get_metadata(file_base):  
    fh = open(file_base+'/metafile_meta_clean','r')
    metadata = {}
    for line in fh:
        items = line.strip().split("\t")
        #metadata[items[0]] = {'lanekey':'1:'+items[0],'direction':items[1],'dna_region':items[2],'project':items[3],'dataset':items[4]}
        metadata['1_'+items[0]] = {  'key':items[0],        'project':items[1],             'dataset':items[2], 
                                'dna_region':items[3],          'domain':items[4],		 		'direction':items[5],
                                'project_title':items[6],       'project_description':items[7], 'dataset_description':items[8], 
                                'env_sample_source_id':items[9]
                                }
    fh.close()
    return metadata
    
    
def gather_and_store_info(myobject):
    """
    
    """
    vamps_cursor    = myobject['vamps_cursor']
    user_cursor     = myobject['user_cursor']
    user            = myobject['user']
    runcode         = myobject['runcode']
    datetime        = str(myobject['datetime'])
    loadtype        = myobject['type']
    upload_source   = myobject['upload_source']
    file_base       = myobject['file_base']
    
    ptitle = " "
    pdescription = " "
    # this is the original count of seqs in the raw sequence file
    # or trimed if uploaded as trimmed
    # do we need it?
    seq_count         = myobject['seq_count']
    if loadtype == 'raw':
        # this is purely for env_sample_source (But useful for anything else we need to pull from the original metadata.csv file)
        metadata_obj = get_metadata(file_base)
        env_sources = {}
        ddescriptions = {}
        pdescriptions = {}
        ptitles = {}
        for r in metadata_obj:
            
            project = metadata_obj[r]['project'][0].capitalize() + metadata_obj[r]['project'][1:]
            #ddescriptions[metadata_obj[r]['project']] = {}
            env_sources[project] = metadata_obj[r]['env_sample_source_id']
            pdescriptions[project] = metadata_obj[r]['project_description']
            ptitles[project] = metadata_obj[r]['project_title']
            if project in ddescriptions:
                ddescriptions[project][metadata_obj[r]['dataset']] = metadata_obj[r]['dataset_description']
            else:
                ddescriptions[project] = {}
    vamps_cursor.execute("select first_name, last_name, email, institution from vamps_auth where user='"+user+"' ")
    (first,last,email,institution) = vamps_cursor.fetchone()
    contact = first+" "+last
    trimseq_select1 ="select project, count(*) as knt from vamps_upload_trimseq where user='"+user+"' and run='"+runcode+"'  \
                            and deleted='0' group by project"
    user_cursor.execute(trimseq_select1)
    #print trimseq_select
    project_rows = user_cursor.fetchall()
    
    
    for p in project_rows:
        project = p[0][0].capitalize() + p[0][1:]
        pcount = str(p[1]) 
        if loadtype == 'trimmed':
            #from CL
            env_source_id = myobject['env_source_id']
            pdescription = ''
            ptitle      = ''
        else:
            env_source_id = env_sources[project] 
            pdescription = pdescriptions[project]
            ptitle      = ptitles[project]
            
            
        q0 = "select vamps_sample_source_id from vamps_sample_source where vamps_sample_source_id='"+env_source_id+"'";
        vamps_cursor.execute(q0)
        envid = vamps_cursor.fetchone()
        if not envid:
            env_source_id = '100'   # unknown
        
        
        info_insert = "insert ignore into vamps_upload_info \
					(project_name, title, description, contact, email, institution,user,edits, \
					upload_date,upload_function,has_tax,seq_count,env_source_id,project_source) \
					values(   '"+project+"',  '"+ptitle+"',  '"+pdescription+"', \
					            '"+contact+"', 	'"+email+"',   '"+institution+"', \
					            '"+user+"',     '"+user+"',   '"+datetime+"', \
					            '"+loadtype+"', '0', 	  '"+pcount+"', \
					            '"+env_source_id+"',  '"+upload_source+"')"
        #parallel tables: update both
        vamps_cursor.execute(info_insert)
        user_cursor.execute(info_insert)
        
        trimseq_select2 = "SELECT dataset, count(*) as dataset_count, entry_date \
						FROM vamps_upload_trimseq \
						WHERE project='"+project+"' and user='"+user+"' and run='"+runcode+"' and deleted='0' \
						GROUP BY dataset"
        user_cursor.execute(trimseq_select2)
						
        dataset_rows = user_cursor.fetchall()
        
        for d in dataset_rows:
            dataset = d[0]
            dcount = str(d[1])
            if loadtype == 'trimmed':
                ddescription = dataset
            else:
                ddescription = ddescriptions[project][dataset]
            #ddescription = " "
		    #
            # now update vamps_projects_datasets_pipe
            #
            pd_pipe_insert = "INSERT ignore INTO vamps_projects_datasets_pipe \
                            (project, dataset, dataset_count, has_tax, date_trimmed, dataset_info) \
                           VALUES ('"+project+"','"+dataset+"','"+dcount+"','0','"+datetime+"','"+ddescription+"')"
            # parallel tables: update both
            vamps_cursor.execute(pd_pipe_insert )	
            user_cursor.execute(pd_pipe_insert )
            
	    #
		# now update vamps_users
	    #
        vamps_cursor.execute("INSERT ignore INTO vamps_users (user, project)  values ('"+user+"','"+project+"')" )	
		
		
    vamps_cursor.close()
    user_cursor.close()
    
if __name__ == '__main__':
    import argparse
    
    # DEFAULTS
    site = 'vampsdev'
    
    user = ''  
    upload_source = 'unknown'
    data_object = {}
    

    
    myusage = """usage: otu_matrix2db.py -im matrixFile -it taxFile [options]
         
         Put user created otus into the vamps database. The OTUs must be in the
         form of a matrix file and a taxonomy file.
         
         where
            -i, --infile The name of the matrix file.  [required]
            
            -t, --upload_type    upload type.   [required]
            
            -p, --project    The name of the sequence file.   [optional]
            
            -d, --dataset     Has to be either usearch, crop or slp. 
                                [default: unknown]
            -s, --source       Percent similarity: generally 3, 6 or 10 percent
                                [default: unknown]    
            -r, --site            vamps or vampsdev.
                                [default: vampsdev]
            -code            upload_id.
                                [default: random number]
            -u, --user       Needed for otu naming and db access.
                             Will be retrieved from .dbconf if not supplied
            
    
    
    """
    parser = argparse.ArgumentParser(description="Gets information (contact,email counts...) from just uploaded and trimmed \
                                                sequences and put it into the appropriate tables on vamps." ,usage=myusage)
    parser.add_argument('-file_base',       required=True, action="store",   dest = "file_base", 
                                                    help = '')                                                 
    parser.add_argument("-t","--upload_type",        required=True,  action="store",   dest = "type", 
                                                    help="raw or trimmed")                                   
    
    
                                                     
    parser.add_argument("-site", "--site",               required=True,  action="store",   dest = "site", 
                                                    help="""database hostname: vamps or vampsdev
                                                        [default: vampsdev]""")  
    parser.add_argument("-rc", "--runcode",      required=True,  action="store",   dest = "runcode", 
                                                    help="""""")  
 
    parser.add_argument("-u", "--user",         required=True,  action="store",   dest = "user", 
                                                    help="user name")  
    
## optional       
    parser.add_argument("-upload_source",       required=False,  action="store",   dest = "upload_source", 
                                                    help="") 
                                          
    parser.add_argument("-count", "--count",      required=True,  action="store",   dest = "seq_count", 
                                                    help=" ")                                                  
    parser.add_argument("-env_id", "--env_id",      required=False,  action="store",   dest = "env_source_id", 
                                                    help=" ")                                                    
    
    logger.info("Starting vamps_info.py")
    args = parser.parse_args()
    
    data_object['datetime'] = str(datetime.date.today())
    data_object['type'] = args.type
    data_object['runcode'] = args.runcode
    data_object['site'] =  args.site
    data_object['user'] =  args.user
    data_object['seq_count']  = args.seq_count
    data_object['file_base']  = args.file_base
    data_object['env_source_id']  = args.env_source_id
    if args.upload_source:
        upload_source = args.upload_source
    data_object['upload_source'] = upload_source
    
    
    if data_object['site'] == 'vamps':
        db_host_vamps   = 'vampsdb'
        db_host_user    = 'bpcdb2'
        db_name_vamps   = 'vamps'
        db_name_user    = 'vamps_user_uploads'
        db_home         = '/xraid2-2/vampsweb/vamps/'
    else:
        db_host_vamps   = 'bpcweb7'
        #db_host_vamps   = 'BPCWeb7'
        db_host_user    = 'bpcdb2'
        db_name_vamps   = 'vamps'
        db_name_user    = 'vampsdev_user_uploads'
        db_home         = '/xraid2-2/vampsweb/vampsdev/'
    
    #print 'vamps',db_host_vamps, db_name_vamps, db_home
    #print 'user',db_host_user, db_name_user, db_home
    vamps_obj = ConMySQL.New(db_host_vamps, db_name_vamps, db_home)
    data_object['vamps_cursor'] = vamps_obj.get_cursor()
    
    user_obj = ConMySQL.New(db_host_user, db_name_user, db_home)
    data_object['user_cursor'] = user_obj.get_cursor()
    
    
    gather_and_store_info(data_object)
    logger.info("Finished:info "+args.runcode)
    
    data_object['vamps_cursor'].close()
    data_object['user_cursor'].close()
    