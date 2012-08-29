
import sys, os, stat
import subprocess
import pipeline.galaxy.fastq_trimmer_by_quality as galaxytrimmer
#from Bio import SeqIO

class IlluminaFiltering:
    """reads an Illumina fastq file and outputs a new filtered file
         Filtering can be on chastity, Ns, mate exists, B-tails, etc.
         Based partially on minoche_filtering_pipeline.pl, but with more controls
         /bioware/illumina_scripts/minoche_filtering_pipeline.pl
         
    Usage:  illumina_filtering.pl -chastity -N -Btail -len min_length -trim trim_length -clip clip_bases -unmated mate.fastq -in input.fastq -failed failed.fastq -out output.fastq -outmate outmate.fastq

      ex:  illumina_filtering.pl -chastity -in PhiX_TTAGGC_L008_R1_001.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.fastq 
           illumina_filtering.pl -chastity -N -in PhiX_TTAGGC_L008_R1_001.fastq -failed PhiX_TTAGGC_L008_R1_001.failed.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.fastq 
           illumina_filtering.pl -unmated PhiX_TTAGGC_L008_R2_001.chaste.N.fastq -outmate PhiX_TTAGGC_L008_R2_001.chaste.N.mated.fastq -in PhiX_TTAGGC_L008_R1_001.chaste.N.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.N.mated.fastq 

 Options:  
           -in          input fastq file
           -out         output fastq file of reads passing quality filters
           -failed      output fastq file of reads failing quality filters
           -unmated     input mate fastq file to compare for removing reads that no longer have a mate
           -outmate     output fastq file for mated reads from -unmated input file

           -chastity    remove reads that do not pass the Illumina chastity filter
           -N           remove reads that contain an ambiguous base
           -Btail       remove reads that contain a Btail
           -unmated     remove reads that no longer have a mate in the mate fastq file
           -len         remove reads that are shorter than minimum length
           -f50         remove reads that have >=34 qual scores < 30 in the first 50 nt (>2/3s)
           -trim        trim remaining reads to trim length
           -clip        clip N bases from the start of remaining reads (for R2)
           -Nx          ignore ambiguous base at position X (seems to be in the v6 primer)
    """
    Name = "IlluminaFiltering"
    def __init__(self, run_object ):
        self.runobj         = run_object
        self.indir          = self.runobj.input_dir
        self.outdir         = os.path.join(self.runobj.output_dir,"trim")
        

    def btails_filter(  self, infile=None,    outfile=None,   qual_score=0,   filter_first50=False,   filter_Ns=False,
                            filter_Nx=0,    length=0,       trim=0, clip=0, failed_fastq=False,     format='sanger',
                            wsize=1,        wstep=1,        trim_ends='53', agg_action='min',       exc_count=0,    
                            score_comp='>=',   keep_zero_length=False
                        ):
        """
        infile=None,            outfile=None, 
        format='sanger',        wsize=1,        wstep=1,            trim_ends='53', 
        agg_action='min',       exc_count=0,    score_comp='>=',    qual_score=0,   
        filter_first50=False,   filter_Ns=False,filter_Nx=0,        failed_fastq=False,
        length=0,               trim=0,         clip=0,             keep_zero_length=False
        """
        
        print "Running illumina Filter"
        in_filename     = infile
        in_filepath = os.path.join(self.indir,in_filename)
        try:
            filebase    = in_filename.split('/')[1].split('.')[0]
        except:
            filebase    = in_filename.split('.')[0]
        out_filename    = filebase+".filtered.fastq"
        out_filepath    = os.path.join(self.outdir, out_filename)
        
        #report_filename = out_filepath + ".btail.report";
        #error_filename = out_filepath + ".btail.errors";
        #print "Running Galaxy fastq_trimmer_by_quality to remove B tails...this can take several hours on a large file...\n";
        #print "Output from fastq_trimmer_by_quality:\n >> report_filename";
        
                    
        galaxytrimmer.trim_by_quality(  infile=in_filepath,         outfile=out_filepath,   qual_score=qual_score,  filter_first50=filter_first50, 
                                        filter_Ns=filter_Ns,        filter_Nx=filter_Nx,    length=length,          trim=trim, clip=clip,
                                        failed_fastq=failed_fastq,  format=format,          wsize=wsize,            wstep=wstep,  
                                        trim_ends=trim_ends,        agg_action=agg_action,  exc_count=exc_count,    score_comp=score_comp,   
                                        keep_zero_length=keep_zero_length
                                    )
         
        #btail_cmd = "python /bioware/galaxy/tools/fastq/fastq_trimmer_by_quality.py -s 1 -t 1 -e 53 -a mean -x 0 -c '>=' -q 3 "+in_filepath+" "+out_filepath+" >> "+report_filename+" 2>"+error_filename
        #subprocess.call(btail_cmd,shell=True)
    
        return out_filename    
        
