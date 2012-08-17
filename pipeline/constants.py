# constants.py
# for mbl sequncing pipeline
############################
#
#  comments '#' and empty lines are ignored
#
#
############ VALIDATION ############################################################################
#"run","data_owner","run_key","lane","dataset","project","tubelabel","barcode","adaptor",
#"dna_region","amp_operator","seq_operator","barcode_index","overlap","insert_size","file_prefix",
#"read_length","primer_suite","first_name","last_name","email","institution","project_title",
#"project_description","funding","env_sample_source","dataset_description"
csv_header_list = {
'illumina' :    ["run",          "data_owner",       "run_key",      "lane",         "project",  "dataset",      
                            "tubelabel",    "barcode",          "adaptor",      "dna_region",   "amp_operator", "seq_operator",     
                            "barcode_index","overlap",          "insert_size",  "file_prefix",  "read_length",  "primer_suite", 
                            "first_name",   "last_name",        "email",        "institution",  "project_title","project_description",
                            "funding",      "env_sample_source","dataset_description" ],
                            
'454' :         [ "run",          "data_owner",       "run_key",      "lane",         "project",  "dataset",
                            "tubelabel",    "barcode",          "adaptor",      "dna_region",   "amp_operator",     "seq_operator",
                            
                            "empcr_operator","direction",       "enzyme",       "concentration", "quant_method",    "domain",               
                            "primer_suite",   "pool",                            
                            
                            "first_name",       "last_name",    "email",        "institution",  "project_title",    "project_description",  
                            "funding",          "env_sample_source","dataset_description" ],

'ion_torrent' : [ "run",          "data_owner",       "run_key",      "lane",         "project",  "dataset",
                            "tubelabel",    "barcode",          "adaptor",      "dna_region",   "amp_operator",     "seq_operator",
                            
                            "empcr_operator","direction",       "enzyme",       "concentration", "quant_method",    "domain",               
                            "primer_suite",   "pool",                            
                            
                            "first_name",       "last_name",    "email",        "institution",  "project_title",    "project_description",  
                            "funding",          "env_sample_source","dataset_description" ]
}

known_platforms = ('illumina','454','ion_torrent','vamps')
primer_suites    = ["bacterialv6suite","Bacterial v6 Suite","bacterial_v6_suite","archaealv6suite","eukaryalv9suite"]
dna_regions      = ["v3", "v3v1", "v3v5", "v3v6", "v4", "v4v5", "v4v6", "v5v3", "v5v4", "v6", "v6a", 
                    "v6v4", "v6v4a", "v6_dutch", "v9", "v9v6"]


VALIDATE_STEP           = "validate"
TRIM_STEP               = "trim"
CHIMERA_STEP            = "chimera"
GAST_STEP               = "gast"
VAMPSUPLOAD             = "vampsupload"
CLUSTER_STEP            = "cluster"
DELETE_RUN              = "delete"
ENV454UPLOAD            = "env454upload"
ENV454RUN_INFO_UPLOAD   = "env454run_info_upload"
STATUS_STEP             = 'status'
CLEAN_STEP              = 'clean'
ILLUMINA_FILES_STEP     = 'illumina_files'

existing_steps = [VALIDATE_STEP, TRIM_STEP, CHIMERA_STEP, GAST_STEP, CLUSTER_STEP, ILLUMINA_FILES_STEP, ENV454RUN_INFO_UPLOAD, ENV454UPLOAD, VAMPSUPLOAD, STATUS_STEP, CLEAN_STEP]

########### RUN CONFIG #############################################################################
input_file_formats = ['sff', 'fasta', 'fasta-mbl', 'fastq', 'fastq-illumina', 'fastq-sanger']


                            
############# defaults for TRIMMING ################################################################
minimumLength   = 50
maximumLength   = ''
minAvgQual      = 30
maxN            = 0

##### DEFAULT PIPELINE RUN ITEMS ##################################################################
# these should contain all the items that are possible for each platform with defaults
#general_run_items = ['baseoutputdir' , 'run' , 'platform','configPath']
# VAMPS as a platform
# if you add another item here be sure to check that
# these other places have been updated
# runconfig.py:initializeFromDictionary() make sure the run object gets the general_config item
# metadata.py:convert_csv_to_ini()        make sure that the item gets written to
#    the NEW ini file

pipeline_run_items = {
'vamps' : { 'dna_region':'v6',
            'domain':'bacteria',
            'project':'test_project',
            'dataset':'test_dataset',
            'from_fasta':False,
            'fasta_file':'',
            'load_vamps_database':True,
            'envsource':'100',
            'use_cluster':True,
            'cluster_nodes':100,
            'site':'vampsdev',
            'user':'',
            'baseoutputdir':'output',
            'require_distal':True,
            'commandline':False,
            'config_file_type':'ini',
            'platform':'vamps'
            },
'illumina' : {'input_file_format':'fastq',
                'compressed':True,
                'database_host':'vampsdev',
                'database_name':'test',
                'platform':'illumina',
                'use_cluster':True,
                'csvPath':'',
                'site':'vampsdev',
                'load_vamps_database':False,
                'anchor_file':'',
                'baseoutputdir':'output',
                'input_dir':'.',
                'primer_file':'',
                'require_distal':True
			},
'454' : {   'input_file_format':'sff',
			'input_file_suffix':'sff',
			'platform':'454',
			'use_cluster':True,
			'csvPath':'',
			'input_dir':'.',
			'baseoutputdir':'output'
		},
'ion_torrent' : {'platform':'ion_torrent',
                'csvPath':'',
                'input_dir':'.',
                'baseoutputdir':'output'
                }
}
# this is the maximum distance from the end of the sequence where script
# will accept a distal primer (if found)
distal_from_end  = 12

complementary_bases = {'A':'T',
                       'T':'A',
                       'C':'G',
                       'G':'C'}

#          for anchor trimming: 
# lowenshtein distance
max_divergence = 0.9

# anchor locations and types
# 0: anchor begin -- where to start looking for an internal anchor (if distal, doesn't really matter)
# 1: anchor end -- where to stop looking for an internal anchor, set to -1 for distal trimming
# 2: minimum length -- delete if sequence is shorter than this length
# 3: trim type is it looking internal = inside for an anchor or distal = at the end
trim_lengths = {
    'v3'    : {'start' : 110, 'end' : -1,  'length' : 110,  'trim_type' : "distal"},
    'v4'    : {'start' : 112, 'end' : -1,  'length' : 112,  'trim_type' : "distal"},
    'v6'    : {'start' : 50,  'end' : -1,  'length' : 50,   'trim_type' : "distal"},
    'v9'    : {'start' : 70,  'end' : -1,  'length' : 70,   'trim_type' : "distal"},
    'v6v4'  : {'start' : -420, 'end' : -525, 'length' : 400,  'trim_type' : "internal"},
    'v6v4a' : {'start' : -325, 'end' : -425, 'length' : 325,  'trim_type' : "internal"},
    'v3v5'  : {'start' : 375, 'end' : 450, 'length' : 350,  'trim_type' : "internal"}
    }
################################################################################################     
    
############# defaults for CHIMERA checking #####################
# if its not in this list chimera checking will be skipped
regions_to_chimera_check = ['v6v4','v3v5','v4v5','v4v6','v3v1','v5v3']
cluster_max_wait                = 1*60*60  # 1 hour
cluster_check_interval          = 2
cluster_initial_check_interval  = 10
################################################################################################  

############# defaults for GAST ################################################################ 
usearch6_cmd        = '/bioware/uclust/usearch6.0.192_i86linux32'
usearch5_cmd        = '/bioware/uclust/usearch5.0.151_i86linux32'
usearch64           = '/bioware/uclust/usearch6.0.192_i86linux64'
fastasampler_cmd    = '/bioware/seqinfo/bin/fastasampler'
calcnodes_cmd       = '/bioware/seqinfo/bin/calcnodes'
mysqlimport_cmd     = '/usr/bin/mysqlimport'
qsub_cmd            = '/usr/local/sge/bin/lx24-amd64/qsub'
clusterize_cmd      = '/bioware/seqinfo/bin/clusterize'
#mothur_cmd          = '/bioware/mothur/mothur'
fastaunique_cmd     = '/bioware/seqinfo/bin/fastaunique'
ref_database_dir    = '/xraid2-2/g454/blastdbs/gast_distributions/'
tax_database_dir    = '/xraid2-2/g454/blastdbs/gast_distributions/'
analysis_dir        = 'analysis'
gast_dir            = 'analysis/gast'
illumina_reads_dir  = 'analysis/perfect_reads'

max_accepts = 10
max_rejects = 0
pctid_threshold = 0.70
majority = 66
cluster_nodes = 0
use_full_length = 0
ignore_terminal_gaps = 0
ignore_all_gaps = 0
max_distance = {'default': 0.30, 'v6': 0.30, 'v6a': 0.30, 'v6v4': 0.25, 'v3v5': 0.25}
#cluster wait
maxwaittime = 1000  # seconds
sleeptime = 5      # seconds
refdbs = {'unknown':'refssu_all',
        'v1v3':'refv1v3',
        'v1v3a':'refv1v3a',
        'v3'    :'refv3',
        'v3a'   :'refv3a',
        'v3v5'  :'refv3v5',
        'v4'    :'refv4',
        'v4v5'  :'refv4v5',
        'v4v6'  :'refv4v6',
        'v6v4'  :'refv4v6',
        'v4v6a' :'refv4v6a',
        'v6v4a' :'refv4v6a',
        'v5'    :'refv5',
        'v6'    :'refv6',
        'v6a'   :'refv6a',
        'v9'    :'refv9',
        'ITS'   :'refITS'       
        }
    

    
########### VAMPS UPLOAD ###########################################################################  
ranks = ('domain','phylum','class','orderx','family','genus','species','strain')
database_tables = {
'vamps_user_uploads': {
            'tax_dc_tbl'      : 'vamps_data_cube_uploads',
            'tax_summed_tbl'  : 'vamps_junk_data_cube_pipe',
            'tax_tbl'         : 'vamps_taxonomy_pipe',
            'sequences_tbl'   : 'vamps_sequences_pipe',
            'export_tbl'      : 'vamps_export_pipe',
            'info_tbl'        : 'vamps_upload_info',
            'datasets_tbl'    : 'vamps_projects_datasets_pipe',
            'users_tbl'       : 'vamps_users'
            },
'vamps_mbl_origin': {
            'tax_dc_tbl'      : 'vamps_data_cube',
            'tax_summed_tbl'  : 'vamps_junk_data_cube',
            'tax_tbl'         : 'vamps_taxonomy',
            'sequences_tbl'   : 'vamps_sequences',
            'export_tbl'      : 'vamps_export',
            'info_tbl'        : 'vamps_projects_info',
            'datasets_tbl'    : 'vamps_projects_datasets',
            'users_tbl'       : 'vamps_users'
            }

}



####################################################################################################
  
    
