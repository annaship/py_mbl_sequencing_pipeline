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
#"project_description","funding","env_sample_source_id","dataset_description"
csv_header_list = {
'illumina' :    ["run",    "data_owner",    "run_key",    "lane",    "dataset",    "project",    "tubelabel",    "barcode",   
                            "adaptor",    "dna_region",    "amp_operator",    "seq_operator",    "barcode_index",    "overlap",    "insert_size",    
                            "read_length",    "primer_suite",    "first_name",    "last_name",    "email",    "institution",    
                            "project_title",    "project_description",    "funding",    "env_sample_source_id",    "dataset_description"],
                            
'454' :         [ "run",          "data_owner",       "run_key",      "lane",         "project",  "dataset",
                            "tubelabel",    "barcode",          "adaptor",      "dna_region",   "amp_operator",     "seq_operator",
                            
                            "empcr_operator","direction",       "enzyme",       "concentration", "quant_method",    "domain",               
                            "primer_suite",   "pool",                            
                            
                            "first_name",       "last_name",    "email",        "institution",  "project_title",    "project_description",  
                            "funding",          "env_sample_source_id","dataset_description" ],

'ion_torrent' : [ "run",          "data_owner",       "run_key",      "lane",         "project",  "dataset",
                            "tubelabel",    "barcode",          "adaptor",      "dna_region",   "amp_operator",     "seq_operator",
                            
                            "empcr_operator","direction",       "enzyme",       "concentration", "quant_method",    "domain",               
                            "primer_suite",   "pool",                            
                            
                            "first_name",       "last_name",    "email",        "institution",  "project_title",    "project_description",  
                            "funding",          "env_sample_source_id","dataset_description" ]
}

known_platforms = ('illumina','454','ion_torrent','vamps')
#primer_suites   = ["bacterialv6suite","bacterial v6 suite","bacterial_v6_suite","archaeal v6 suite","archaealv6suite","eukaryalv9suite","bacterial v4-v5 suite"]
primer_suites   = ["Archaeal V6 Suite", "Archaeal V6-V4 Suite", "Bacterial V3 Suite", "Bacterial V3-V1 Suite", 
                   "Bacterial V3-V5 Suite", "Bacterial V4-V5 Suite", "Bacterial V4-V6 Suite", "Bacterial V5-V3 Suite", 
                   "Bacterial V6 Suite", "Bacterial V6-V4 Suite", "CDSIII", "Eukaryal V9 Suite", "eukv9_1380", 
                   "eukv9_1389", "Fungal ITS1 Suite", "HMP V3-V1 Suite", "HMP V5-V3 Suite", "hmpv3v1", "hmpv5v3", 
                   "Relman", "ti_v3v6", "ti_v6", "topo", "v6v4", "v6_dutch", "Vibrio V4"]
dna_regions     = ["v3", "v3v1", "v3v5", "v3v6", "v4", "v4v5", "v4v6", "v5v3", "v5v4", "v6", "v6a", 
                    "v6v4", "v6v4a", "v6_dutch", "v9", "v9v6", "its1"]


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
ILLUMINA_FILES_DEM_STEP = 'illumina_files_demultiplex_only'
ILLUMINA_CHIMERA_ONLY_STEP = 'illumina_chimera_only'
ILLUMINA_CHIMERA_AFTER_CLUSTER = 'illumina_chimera_after_cluster'

existing_steps = [VALIDATE_STEP, TRIM_STEP, CHIMERA_STEP, GAST_STEP, CLUSTER_STEP, ILLUMINA_CHIMERA_ONLY_STEP, ILLUMINA_CHIMERA_AFTER_CLUSTER, ILLUMINA_FILES_DEM_STEP, ILLUMINA_FILES_STEP, ENV454RUN_INFO_UPLOAD, ENV454UPLOAD, VAMPSUPLOAD, STATUS_STEP, CLEAN_STEP]

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
                'require_distal':True,
                'do_perfect':True,
                'lane_name':''
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
regions_to_chimera_check = ['v6v4','v3v5','v4v5','v4v6','v3v1','v1v3','v5v3']
cluster_max_wait                = 1*60*60  # 1 hour
cluster_check_interval          = 2
cluster_initial_check_interval  = 10

################################################################################################  

############# directories ################################################################ 
"""
Output path example: /xraid2-2/g454/run_new_pipeline/illumina/miseq/20121025/analysis/gast
"""
#output data
#subdirs = ['analysis_dir', 'gast_dir', 'reads_overlap_dir', 'vamps_upload_dir', 'chimera_dir', 'trimming_dir']

#under rundate
subdirs = {'analysis_dir' : 'analysis',
#under analysis
'gast_dir'          : 'gast',
'reads_overlap_dir' : 'reads_overlap',
'vamps_upload_dir'  : 'vamps_upload',
'chimera_dir'       : 'chimera',
'trimming_dir'      : 'trimming'}

#root:
#VAMPS users
output_root_vampsdev_users = '/xraid2-2/vampsweb/vampsdev/tmp/'
output_root_vamps_users = '/xraid2-2/vampsweb/vamps/tmp/'

#MBL
output_root_mbl         = '/xraid2-2/g454/run_new_pipeline/'

output_root_mbl_local   = '/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test' 
# cmd_path_local          = "/Users/ashipunova/bin/illumina-utils/"
cmd_path_local          = "/Users/ashipunova/bin/illumina-utils/illumina-utils/scripts/"

############# defaults for illumina ################################################################ 
perfect_overlap_cmd       = "analyze-illumina-v6-overlaps"
perfect_overlap_cmd_local = cmd_path_local + perfect_overlap_cmd    
partial_overlap_cmd       = "merge-illumina-pairs"
partial_overlap_cmd_local = cmd_path_local + partial_overlap_cmd          
filter_mismatch_cmd       = "filter-merged-reads" 
filter_mismatch_cmd_local = cmd_path_local + filter_mismatch_cmd          

################################################################################################  
filtered_suffix    = 'MERGED-MAX-MISMATCH-3' #result of filter-merged-reads
nonchimeric_suffix = "nonchimeric.fa"
############# defaults for GAST ################################################################
usearch_cmd        = 'usearch'
usearch6_cmd       = 'usearch6'
# usearch6_cmd        = '/bioware/uclust/usearch6'
#usearch5_cmd        = '/bioware/uclust/usearch5.2.32_i86linux32'
#usearch64           = '/bioware/uclust/usearch6.0.192_i86linux64'
fastasampler_cmd    = '/bioware/seqinfo/bin/fastasampler'
calcnodes_cmd       = '/bioware/seqinfo/bin/calcnodes'
mysqlimport_cmd     = '/usr/bin/mysqlimport'
qsub_cmd            = '/usr/local/sge/bin/lx24-amd64/qsub'
clusterize_cmd      = '/bioware/seqinfo/bin/clusterize'
#mothur_cmd          = '/bioware/mothur/mothur'
#fastaunique_cmd     = '/bioware/seqinfo/bin/fastaunique'
fastaunique_cmd     = '/bioware/seqinfo/bin/fastaunique'
#local commands
usearch6_cmd_local     = cmd_path_local + 'usearch'
fastasampler_cmd_local = cmd_path_local + 'fastasampler'
calcnodes_cmd_local    = cmd_path_local + 'calcnodes'
mysqlimport_cmd_local  = '/usr/local/mysql/bin/mysqlimport'
fastaunique_cmd_local  = cmd_path_local + 'fastaunique'
ref_database_dir_local = cmd_path_local
tax_database_dir_local = cmd_path_local

ref_database_dir    = '/xraid2-2/g454/blastdbs/gast_distributions/'
tax_database_dir    = '/xraid2-2/g454/blastdbs/gast_distributions/'
vamps_ref_database_dir    = '/xraid2-2/vampsweb/blastdbs/'
vamps_tax_database_dir    = '/xraid2-2/vampsweb/blastdbs/'

max_accepts = 10
max_rejects = 0
pctid_threshold = 0.70
majority = 66
cluster_nodes = 100
use_full_length = 0
ignore_terminal_gaps = 0
ignore_all_gaps = 0
max_gast_distance = {'default': 0.30, 'v6': 0.30, 'v6a': 0.30, 'v6v4': 0.25, 'v3v5': 0.25}
#cluster wait
maxwaittime = 50000  # seconds
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
        'its1'   :'refits1'       
        }

""" Use '_6' for usearch6 """    
chimera_checking_refdb_6     = '/xraid2-2/g454/blastdbs/rRNA16S.gold.udb'
chimera_checking_refdb       = '/xraid2-2/g454/blastdbs/rRNA16S.gold.fasta'
chimera_checking_its_refdb_6 = '/xraid2-2/g454/blastdbs/fungalITS.udb'
chimera_checking_its_refdb   = '/xraid2-2/g454/blastdbs/fungalITS.fa'
chimera_checking_abskew      = '1.9'    
    
########### VAMPS UPLOAD ###########################################################################  
ranks = ('domain','phylum','class','orderx','family','genus','species','strain')
domains = ('Archaea','Bacteria','Eukarya','Organelle','Unknown')
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
  
    
