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
from collections import defaultdict

csv_header_list = {
'illumina' :    ["run", "data_owner", "run_key", "lane", "dataset", "project", "tubelabel", "barcode",
 "adaptor", "dna_region", "amp_operator", "seq_operator", "barcode_index", "overlap", "insert_size",
 "read_length", "platform", "primer_suite", "first_name", "last_name", "email", "institution",
 "project_title", "project_description", "funding", "env_sample_source_id", "dataset_description"],
'miseq' :    ["run", "data_owner", "run_key", "lane", "dataset", "project", "tubelabel", "barcode",
 "adaptor", "dna_region", "amp_operator", "seq_operator", "barcode_index", "overlap", "insert_size",
 "read_length", "platform", "primer_suite", "first_name", "last_name", "email", "institution",
 "project_title", "project_description", "funding", "env_sample_source_id", "dataset_description"],
'hiseq' :    ["run", "data_owner", "run_key", "lane", "dataset", "project", "tubelabel", "barcode",
 "adaptor", "dna_region", "amp_operator", "seq_operator", "barcode_index", "overlap", "insert_size",
 "read_length", "platform", "primer_suite", "first_name", "last_name", "email", "institution",
 "project_title", "project_description", "funding", "env_sample_source_id", "dataset_description"],
'nextseq' :    ["run", "data_owner", "run_key", "lane", "dataset", "project", "tubelabel", "barcode",
 "adaptor", "dna_region", "amp_operator", "seq_operator", "barcode_index", "overlap", "insert_size",
 "read_length", "platform", "primer_suite", "first_name", "last_name", "email", "institution",
 "project_title", "project_description", "funding", "env_sample_source_id", "dataset_description"],

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

known_platforms = ('illumina','454','ion_torrent','vamps','hiseq','miseq','nextseq')
illumina_list = ['hiseq', 'miseq', 'nextseq']
#primer_suites   = ["bacterialv6suite","bacterial v6 suite","bacterial_v6_suite","archaeal v6 suite","archaealv6suite","eukaryalv9suite","bacterial v4-v5 suite"]
# todo: take from db!
primer_suites   = ["archaeal v6 suite", "archaeal v4 suite", "archaeal v6mod suite", "archaeal v6-v4 suite", "bacterial v3 suite", "bacterial v3-v1 suite",
                   "bacterial v3-v5 suite", "archaeal v4-v5 suite", "bacterial v4-v5 suite", "bacterial v4-v6 suite", "bacterial v4 suite", "bacterial v5-v3 suite",
                   "bacterial v6 suite", "bacterial v6-v4 suite", "cdsiii", "eukaryal v9 suite", "eukv9_1380",
                   "eukv9_1389", "fungal its1 suite", "hmp v3-v1 suite", "hmp v5-v3 suite", "hmpv3v1", "hmpv5v3",
                   "relman", "ti_v3v6", "ti_v6", "topo", "v6v4", "v6_dutch", "vibrio v4", "eukaryal v4 suite", "eukaryal hssu suite", "eukaryal hlsu suite"]
dna_regions     = ["v3", "v3v1", "v3v5", "v3v6", "v4", "v4v5", "v4v6", "v5v3", "v5v4", "v6", "v6a",
                    "v6v4", "v6v4a", "v6_dutch", "v9", "v9v6", "its1", "hlsu", "hssu"]

#K = [G,T]
#M = [A,C]
#R = [A,G]
#S = [G,C]
#W = [A,T]
#Y = [C,T]
#todo: add to env454 "combined_primer"
primers_dict = defaultdict(dict)
primers_dict["archaeal v6mod suite"]["proximal_primer"]  = "AATTGGCGGGGGAGCAC"
primers_dict["archaeal v6mod suite"]["distal_primer"]    = "GCCATGCACC[A,T]CCTCT"
primers_dict["archaeal v4-v5 suite"]["proximal_primer"]  = "G[C,T][C,T]TAAA..[A,G][C,T][C,T][C,T]GTAGC"
primers_dict["archaeal v4-v5 suite"]["distal_primer"]    = "CCGGCGTTGA.TCCAATT"
primers_dict["bacterial v4-v5 suite"]["proximal_primer"] = "CCAGCAGC[C,T]GCGGTAA."
primers_dict["bacterial v4-v5 suite"]["distal_primer"]   = "CCGTC[A,T]ATT[C,T].TTT[G,A]A.T"
primers_dict["fungal its1 suite"]["proximal_primer"]     = "GTAAAAGTCGTAACAAGGTTTC"
primers_dict["fungal its1 suite"]["distal_primer"]       = "GTTCAAAGA[C,T]TCGATGATTCAC"

# hard coded in analyze-illumina-v6-overlaps
primers_dict["archaeal v6 suite"]["proximal_primer"]  = "A.TCAACGCCGG"
primers_dict["archaeal v6 suite"]["distal_primer"]    = "G[A,T]GGT[G,A]"
primers_dict["bacterial v6 suite"]["proximal_primer"] = "AGGTG."
primers_dict["bacterial v6 suite"]["distal_primer"]   = "[A,G]AACCT[CT]A.C"
primers_dict["eukaryal v4 suite"]["proximal_primer"] = "CCAGCA[C,G]C[C,T]GCGGTAATTCC"
# primers_dict["eukaryal v4 suite"]["distal_primer"]   = "ACTTTCGTTCTTGAT[C,T][A,G]A"
primers_dict["eukaryal v4 suite"]["distal_primer"]   = "ACTTTCGTTCTTGAT[C,T][A,G][A,G]"

primers_dict["archaeal v4 suite"]["proximal_primer"]  = "GTGTG[CT]CAGC[AC]GCCGCGGTAA"
primers_dict["archaeal v4 suite"]["distal_primer"]    = "CCGGACTAC[ACGT][ACG]GGGT[AT]TCTAAT"
primers_dict["bacterial v4 suite"]["proximal_primer"] = "GTGTG[CT]CAGC[AC]GCCGCGGTAA"
primers_dict["bacterial v4 suite"]["distal_primer"]   = "CCGGACTAC[ACGT][ACG]GGGT[AT]TCTAAT"

primers_dict["eukaryal hssu suite"]["proximal_primer"] = "GCGGTAATTCCAGCTCCA"
primers_dict["eukaryal hssu suite"]["distal_primer"]   = "GATCAGTGAAAACATCCCTGG"
primers_dict["eukaryal hlsu suite"]["proximal_primer"] = "GGT[AG]TCGGAGA[AG]GGTGAGAATCC"
primers_dict["eukaryal hlsu suite"]["distal_primer"]   = "TCAGACTCCTTGGTCCGTGTTTCT"




VALIDATE_STEP           = "validate"
TRIM_STEP               = "trim"
CHIMERA_STEP            = "chimera"
GAST_STEP               = "gast"
VAMPSUPLOAD             = "vampsupload"
NEW_VAMPS               = "new_vamps"
CLUSTER_STEP            = "cluster"
DELETE_RUN              = "delete"
ENV454UPLOAD            = "env454upload"
ENV454UPLOAD_NO_SEQ    = "env454upload_no_seq"
ENV454RUN_INFO_UPLOAD   = "env454run_info_upload"
STATUS_STEP             = 'status'
CLEAN_STEP              = 'clean'
ILLUMINA_FILES_STEP     = 'illumina_files'
ILLUMINA_FILES_DEM_STEP = 'illumina_files_demultiplex_only'
ILLUMINA_CHIMERA_ONLY_STEP = 'illumina_chimera_only'
ILLUMINA_CHIMERA_AFTER_CLUSTER = 'illumina_chimera_after_cluster'

existing_steps = [VALIDATE_STEP, TRIM_STEP, CHIMERA_STEP, GAST_STEP, CLUSTER_STEP, ILLUMINA_CHIMERA_ONLY_STEP, ILLUMINA_CHIMERA_AFTER_CLUSTER, ILLUMINA_FILES_DEM_STEP, ILLUMINA_FILES_STEP, ENV454RUN_INFO_UPLOAD, ENV454UPLOAD, ENV454UPLOAD_NO_SEQ, VAMPSUPLOAD, NEW_VAMPS, STATUS_STEP, CLEAN_STEP]

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
'miseq' : {'input_file_format':'fastq',
                'compressed':True,
                'database_host':'vampsdev',
                'database_name':'test',
                'platform':'miseq',
                'use_cluster':True,
                'csvPath':'',
                'site':'vampsdev',
                'load_vamps_database':False,
                'anchor_file':'',
                'baseoutputdir':'output',
                'input_dir':'.',
                'primer_file':'',
                'require_distal':True,
                'do_perfect': False,
                'lane_name':'1'
			},
'nextseq' : {'input_file_format':'fastq',
                'compressed':True,
                'database_host':'vampsdev',
                'database_name':'test',
                'platform':'nextseq',
                'use_cluster':True,
                'csvPath':'',
                'site':'vampsdev',
                'load_vamps_database':False,
                'anchor_file':'',
                'baseoutputdir':'output',
                'input_dir':'.',
                'primer_file':'',
                'require_distal':True,
                'do_perfect': False,
                'lane_name':'1'
			},
'hiseq' : {'input_file_format':'fastq',
                'compressed':True,
                'database_host':'vampsdev',
                'database_name':'test',
                'platform':'hiseq',
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
regions_to_chimera_check = ['v6v4','v3v5','v4v5','v4v6','v3v1','v1v3','v5v3', 'ITS1', 'v4', 'HSSU', 'HLSU']
cluster_max_wait                = 1*60*60  # 1 hour
cluster_check_interval          = 2
cluster_initial_check_interval  = 10

################################################################################################

############# directories ################################################################
"""
Output path example: /groups/g454/run_new_pipeline/illumina/miseq/20121025/analysis/gast
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
output_root_vampsdev_users      = '/groups/vampsweb/vampsdev/tmp/'
output_root_vamps_users         = '/groups/vampsweb/vamps/tmp/'
#New VAMPS
output_root_newvamps_users      = '/groups/vampsweb/vampsdev_user_data/'
new_vamps_database = 'vamps2'
#MBL
output_root_mbl         = '/groups/g454/run_new_pipeline/'

output_root_mbl_local   = '/Users/ashipunova/BPC/py_mbl_sequencing_pipeline/test'
# cmd_path_local          = "/Users/ashipunova/bin/illumina-utils/", "/Users/ashipunova/bin/illumina-utils/illumina-utils/scripts/"
cmd_path_local          = "/usr/local/bin/" 
cmd_path_vamps          = "/groups/vampsweb/seqinfobin/"
############# defaults for illumina ################################################################
# perfect_overlap_cmd       = "iu-analyze-v6-complete-overlaps"
overlap_command           = "iu-merge-pairs"
perfect_overlap_cmd       = overlap_command
perfect_overlap_cmd_local = perfect_overlap_cmd
partial_overlap_cmd       = overlap_command
partial_overlap_cmd_local = partial_overlap_cmd
trim_primers_cmd          = "iu-trim-V6-primers"
filter_mismatch_cmd       = "iu-filter-merged-reads"
filter_mismatch_cmd_local = filter_mismatch_cmd

################################################################################################
trimming_length    = 251
filtered_suffix    = "MERGED-MAX-MISMATCH-3" #result of filter-merged-reads
nonchimeric_suffix = "nonchimeric.fa"
unique_suffix      = "unique"
############# defaults for GAST ################################################################
# usearch_cmd        = '/bioware/usearch/5.2.236/x86/usearch'     #usearch5 32bit
# usearch6_cmd       = '/bioware/usearch/6.0.217/x86/usearch'     #usearch6 32bit
# usearch64_cmd      = '/bioware/usearch/7.0.1090/x86_64/usearch'  #usearch6 64bit, for non-parallel execution ONLY
# usearch6_cmd_local     = cmd_path_local + 'usearch'
# use the whole path here please
usearch_cmd        = '/bioware/seqinfo/bin/vsearch'
usearch6_cmd       = usearch_cmd
usearch64_cmd      = usearch_cmd
usearch6_cmd_local = cmd_path_local + 'vsearch'
py_pipeline_base = '/bioware/seqinfo/bin/python_pipeline/py_mbl_sequencing_pipeline/pipeline/'

fastasampler_cmd    = '/bioware/seqinfo/bin/fastasampler'
calcnodes_cmd       = '/bioware/seqinfo/bin/calcnodes'
mysqlimport_cmd     = '/usr/bin/mysqlimport'
qsub_cmd            = '/usr/local/sge/bin/lx24-amd64/qsub'
clusterize_cmd      = '/bioware/seqinfo/bin/clusterize'
#mothur_cmd          = '/bioware/mothur/mothur'
#fastaunique_cmd     = '/bioware/seqinfo/bin/fastaunique'
fastaunique_cmd     = '/bioware/seqinfo/bin/fastaunique'
#local commands
fastasampler_cmd_local = cmd_path_local + 'fastasampler'
calcnodes_cmd_local    = cmd_path_local + 'calcnodes'
mysqlimport_cmd_local  = '/usr/local/mysql/bin/mysqlimport'
fastaunique_cmd_local  = cmd_path_local + 'fastaunique'
ref_database_dir_local = cmd_path_local
tax_database_dir_local = cmd_path_local
# VAMPS commands
fastasampler_cmd_vamps  = cmd_path_vamps + 'fastasampler'
fastaunique_cmd_vamps   = cmd_path_vamps + 'fastaunique'
calcnodes_cmd_vamps     = cmd_path_vamps + 'calcnodes'
vsearch_cmd_vamps       = cmd_path_vamps + 'vsearch'
tophit_cmd_vamps        = cmd_path_vamps + 'clustergast_tophit'
clusterize_cmd_vamps    = cmd_path_vamps + 'clusterize_vamps'
py_pipeline_base_vamps = '/groups/vampsweb/py_mbl_sequencing_pipeline/pipeline/'


#####
ref_database_dir    = '/groups/g454/blastdbs/gast_distributions/'
tax_database_dir    = '/groups/g454/blastdbs/gast_distributions/'
# these are symlinks to above
vamps_ref_database_dir    = '/groups/vampsweb/blastdbs/gast_distributions/'
vamps_tax_database_dir    = '/groups/vampsweb/blastdbs/gast_distributions/'
#vamps_ref_database_dir    = '/groups/vampsweb/blastdbs/'
#vamps_tax_database_dir    = '/groups/vampsweb/blastdbs/'

max_accepts         = 10
max_rejects         = 0
pctid_threshold     = 0.70
majority            = 66
cluster_nodes       = 100
use_full_length     = 0
ignore_terminal_gaps= 0
ignore_all_gaps     = 0
max_gast_distance   = {'default': 0.30, 'v6': 0.30, 'v6a': 0.30, 'v6v4': 0.25, 'v3v5': 0.25}
#cluster wait
maxwaittime         = 50000  # 50,000 seconds == 833 minutes == 13.9 hours
sleeptime           = 3      # seconds
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
# chimera_checking_refdb_6     = '/groups/g454/blastdbs/rRNA16S.gold.udb'
chimera_checking_refdb_6     = '/groups/g454/blastdbs/rRNA16S.gold.fasta'
chimera_checking_refdb       = '/groups/g454/blastdbs/rRNA16S.gold.fasta'
chimera_checking_its_refdb_6 = '/groups/g454/blastdbs/fungalITS.udb'
chimera_checking_its_refdb   = '/groups/g454/blastdbs/fungalITS.fa'
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


