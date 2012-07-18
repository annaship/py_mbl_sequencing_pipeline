# constants.py
# for mbl sequncing pipeline
############################
#
#  comments '#' and empty lines are ignored
#
#
############ VALIDATION ############################################################################
csv_header_list_illumina = ["run","data_owner","run_key","lane","dataset","project","tubelabel","barcode","adaptor","dna_region",
                            "amp_operator","seq_operator","barcode_index","overlap","insert_size","file_prefix","read_length",
                            "primer_suite","first_name","last_name","email","institution","project_title","project_description",
                            "funding","env_sample_source","dataset_description" ]
csv_header_list_454 = []                        
primer_suites    = ["bacterialv6suite","bacterial_v6_suite","archaealv6suite","eukaryalv9suite"]
dna_regions      = ["v3", "v3v1", "v3v5", "v3v6", "v4", "v4v5", "v4v6", "v5v3", "v5v4", "v6", "v6a", 
                    "v6v4", "v6v4a", "v6_dutch", "v9", "v9v6"]



########### RUN CONFIG #############################################################################
input_file_formats = ['sff', 'fasta', 'fasta-mbl', 'fastq', 'fastq-illumina', 'fastq-sanger']


                            
############# defaults for TRIMMING ################################################################
minimumLength   = 50
maximumLength   = ''
minAvgQual      = 30
maxN            = 0



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
usearch = '/bioware/uclust/usearch'
fastasampler = '/bioware/seqinfo/bin/fastasampler'
calcnodes = '/bioware/seqinfo/bin/calcnodes'
mysqlimport = '/usr/bin/mysqlimport'
#qsub = '/usr/local/sge/bin/lx24-amd64/qsub'
clusterize = '/bioware/seqinfo/bin/clusterize'
max_accepts = 10
max_rejects = 0
pctid_threshold = 0.70
majority = 66
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
################################################################################################      

  
    
