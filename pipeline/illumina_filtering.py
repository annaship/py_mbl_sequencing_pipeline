#!/bin/env python

# this has been adaped for the MBL pipeline by Andy Voorhis
# the original is from Galaxy
# written by Dan Blankenberg
# and located here:  /bioware/galaxy/tools/fastq/fastq_trimmer_by_quality.py

import sys, os, stat
import subprocess
from pipeline.pipelinelogging import logger
from pipeline.galaxy.fastq import fastqReader, fastqWriter

def mean( score_list ):
    return float( sum( score_list ) ) / float( len( score_list ) )

ACTION_METHODS = { 'min':min, 'max':max, 'sum':sum, 'mean':mean }

class IlluminaFiltering:
    """From Sue's original perl script called illumina_filtering:
    
    reads an Illumina fastq file and outputs a new filtered file
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
        self.outdir         = os.path.join(self.runobj.output_dir,"illumina_filtered")
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
            
    def compare( self, aggregated_value, operator, threshold_value ):
        if operator == '>':
            return aggregated_value > threshold_value
        elif operator == '>=':
            return aggregated_value >= threshold_value
        elif operator == '==':
            return aggregated_value == threshold_value
        elif operator == '<':
            return aggregated_value < threshold_value
        elif operator == '<=':
            return aggregated_value <= threshold_value
        elif operator == '!=':
            return aggregated_value != threshold_value

    def exclude( self, value_list, exclude_indexes ):
        rval = []
        for i, val in enumerate( value_list ):
            if i not in exclude_indexes:
                rval.append( val )
        return rval

    def exclude_and_compare( self, aggregate_action, aggregate_list, operator, threshold_value, exclude_indexes = None ):
        if not aggregate_list or self.compare( aggregate_action( aggregate_list ), operator, threshold_value ):
            return True
        if exclude_indexes:
            for exclude_index in exclude_indexes:
                excluded_list = self.exclude( aggregate_list, exclude_index )
                if not excluded_list or self.compare( aggregate_action( excluded_list ), operator, threshold_value ):
                    return True
        return False
    
    def trim_by_quality(self, infile=None, 
                        format='sanger',        wsize=1,        wstep=1,            trim_ends='53', 
                        agg_action='min',       exc_count=0,    score_comp='>=',    qual_score=0,   
                        filter_first50=False,   filter_Ns=False,filter_Nx=0,        failed_fastq=False,
                        length=0,               trim=0,         clip=0,             keep_zero_length=False):
        #format
        window_size         = wsize
        window_step         = wstep
        #trim_ends
        aggregation_action  = agg_action
        exclude_count       = exc_count
        score_comparison    = score_comp
        quality_score       = qual_score
        filter_length       = length
        trim_length         = trim
        clip_length         = clip
        if not infile:
            sys.exit( "illumina_fastq_trimmer: Need to specify an input file" )
        
        if window_size < 1:
            sys.exit( 'illumina_fastq_trimmer: You must specify a strictly positive window size' )
        
        if window_step < 1:
            sys.exit( 'illumina_fastq_trimmer: You must specify a strictly positive step size' )
            
        print "\nRunning illumina Filtering"
        
        in_filepath = os.path.join(self.indir,infile)
        try:
            filebase    = infile.split('/')[1].split('.')[0]
        except:
            filebase    = infile.split('.')[0]
            
        out_filename    = filebase+".filtered.fastq"
        out_filepath    = os.path.join(self.outdir, out_filename)
        
        
        
            
        #determine an exhaustive list of window indexes that can be excluded from aggregation
        exclude_window_indexes = []
        last_exclude_indexes = []
        for exclude_count in range( min( exclude_count, window_size ) ):
            if last_exclude_indexes:
                new_exclude_indexes = []
                for exclude_list in last_exclude_indexes:
                    for window_index in range( window_size ):
                        if window_index not in exclude_list:
                            new_exclude = sorted( exclude_list + [ window_index ] )
                            if new_exclude not in exclude_window_indexes + new_exclude_indexes:
                                new_exclude_indexes.append( new_exclude )
                exclude_window_indexes += new_exclude_indexes
                last_exclude_indexes = new_exclude_indexes
            else:
                for window_index in range( window_size ):
                    last_exclude_indexes.append( [ window_index ] )
                exclude_window_indexes = list( last_exclude_indexes )
        out = fastqWriter( open( out_filepath, 'wb' ), format = format )
        action = ACTION_METHODS[ aggregation_action ]
        if failed_fastq:
            fail = fastqWriter( open( out_filepath+'.failed', 'wb' ), format = format )
        num_reads = None
        num_reads_excluded = 0
        count_of_unchaste = 0
        count_of_trimmed  = 0
        count_of_first50  = 0
        count_of_Ns  = 0
        if self.runobj.compressed:
            import gzip
            try:
                logger.info( "illumina_filtering: opening compressed file: "+in_filepath)
                fp = gzip.open( in_filepath )
            except:
                logger.info( "illumina_filtering: opening uncompressed file: "+in_filepath)
                fp = open( in_filepath )
        else:
            logger.info(  "illumina_filtering: opening uncompressed file: "+in_filepath)
            fp = open( in_filepath )
        for num_reads, fastq_read in enumerate( fastqReader( fp, format = format ) ):
            ############################################################################################
            # Put chastity code here
            #print fastq_read.identifier
            seq = fastq_read.get_sequence()
            
            desc_items = fastq_read.identifier.split(':')
                
            if desc_items[7] == 'Y':
                count_of_unchaste += 1
                #print 'failed chastity'
                if failed_fastq:
                    fail.write( fastq_read )
                continue
            
            # Filter reads with ambiguous bases
            if filter_Ns:                
                countN = seq.count('N')
                if countN > 1 or (countN == 1 and seq[filter_Nx-1:filter_Nx] != 'N'):
                    #print 'failed Ns',infile
                    count_of_Ns += 1
                    if failed_fastq:
                        fail.write( fastq_read )
                    continue
            
            
            
            # Filter reads below first 50 base quality
            if filter_first50:                
                first50 = 50
                first50_maxQ = 30
                first50_maxQ_count = 34
                
                quals = fastq_read.get_decimal_quality_scores()[:first50]
                count_lt30 = 0
                
                for q in quals:
                    if q < first50_maxQ:
                        count_lt30 += 1
                if count_lt30 >= first50_maxQ_count:
                    #print 'failed first50'
                    if failed_fastq:
                        fail.write( fastq_read )
                    count_of_first50 += 1
                    continue
                    
            ##### END CHASTITY #####################
            ############################################################################################
            ##### START Btails CODE ################
            quality_list = fastq_read.get_decimal_quality_scores()
            
            for trim_end in trim_ends:
                
                
                if trim_end == '5':
                    lwindow_position = 0 #left position of window
                    while True:
                        if lwindow_position >= len( quality_list ):
                            fastq_read.sequence = ''
                            fastq_read.quality = ''
                            break
                        if self.exclude_and_compare( action, quality_list[ lwindow_position:lwindow_position + window_size ], score_comparison, quality_score, exclude_window_indexes ):
                            fastq_read = fastq_read.slice( lwindow_position, None )
                            break
                        lwindow_position += window_step
                else:
                    rwindow_position = len( quality_list ) #right position of window
                    while True:
                        lwindow_position = rwindow_position - window_size #left position of window
                        if rwindow_position <= 0 or lwindow_position < 0:
                            fastq_read.sequence = ''
                            fastq_read.quality = ''
                            break
                        if self.exclude_and_compare( action, quality_list[ lwindow_position:rwindow_position ], score_comparison, quality_score, exclude_window_indexes ):
                            fastq_read = fastq_read.slice( None, rwindow_position )
                            break
                        rwindow_position -= window_step
                        
            ######## END Btails CODE ###############################            
            ############################################################################################
            # put  length/trim/clip code here
            quality_list = fastq_read.get_decimal_quality_scores()
            
            if filter_length:
                if len(quality_list) < filter_length:
                    print 'failed length'
                    if failed_fastq:
                        fail.write( fastq_read )
                    continue
            
            # Trim initial bases -- remove first 10 bases from read 2   
            if clip_length:
                # remove from the front:
                fastq_read = fastq_read.slice( clip_length, None )
                count_of_trimmed += 1
                
            # Trim to max length -- read 2 trim to 90.
            if trim_length:
                if len(quality_list) > trim_length:
                    # remove from the end:
                    fastq_read = fastq_read.slice( None, len(fastq_read.get_sequence()) - trim_length )
                    count_of_trimmed += 1
            
            
            if keep_zero_length or len( fastq_read ):
                out.write( fastq_read )
            else:
                num_reads_excluded += 1
        out.close()
        if failed_fastq:
            fail.close()
        print "file:",infile
        print 'count_of_trimmed             (for length):', count_of_trimmed
        print 'count_of_first50 (avg first50 quals < 34):', count_of_first50
        print "count_of_unchaste             ('Y' in id):", count_of_unchaste
        print 'count_of_Ns                (reads with N):', count_of_Ns
        if num_reads is None:
            print "No valid FASTQ reads could be processed."
        else:
            print "%i FASTQ reads were processed." % ( num_reads + 1 )
        if num_reads_excluded:
            print "%i reads of zero length were excluded from the output." % num_reads_excluded

        return out_filename
        
        
if __name__ == "__main__": trim_by_quality()
