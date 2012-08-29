#!/bin/env python

# this has been adaped for the MBL pipeline by Andy Voorhis
# the original is from Galaxy
# written by Dan Blankenberg
# and located here:  /bioware/galaxy/tools/fastq/fastq_trimmer_by_quality.py

import sys
#from optparse import OptionParser
from pipeline.galaxy.fastq import fastqReader, fastqWriter

def mean( score_list ):
    return float( sum( score_list ) ) / float( len( score_list ) )

ACTION_METHODS = { 'min':min, 'max':max, 'sum':sum, 'mean':mean }

def compare( aggregated_value, operator, threshold_value ):
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

def exclude( value_list, exclude_indexes ):
    rval = []
    for i, val in enumerate( value_list ):
        if i not in exclude_indexes:
            rval.append( val )
    return rval

def exclude_and_compare( aggregate_action, aggregate_list, operator, threshold_value, exclude_indexes = None ):
    if not aggregate_list or compare( aggregate_action( aggregate_list ), operator, threshold_value ):
        return True
    if exclude_indexes:
        for exclude_index in exclude_indexes:
            excluded_list = exclude( aggregate_list, exclude_index )
            if not excluded_list or compare( aggregate_action( excluded_list ), operator, threshold_value ):
                return True
    return False
    
def trim_by_quality(infile=None,            outfile=None, 
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
        
    if not infile or not outfile:
        sys.exit( "Need to specify an input file and an output file" )
    
    if window_size < 1:
        sys.exit( 'You must specify a strictly positive window size' )
    
    if window_step < 1:
        sys.exit( 'You must specify a strictly positive step size' )
        
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
    out = fastqWriter( open( outfile, 'wb' ), format = format )
    action = ACTION_METHODS[ aggregation_action ]
    if failed_fastq:
        fail = fastqWriter( open( outfile+'.failed', 'wb' ), format = format )
    num_reads = None
    num_reads_excluded = 0
    count_of_unchaste = 0
    count_of_trimmed  = 0
    count_of_first50  = 0
    count_of_Ns  = 0
    for num_reads, fastq_read in enumerate( fastqReader( open( infile ), format = format ) ):
        
        ############################################################################################
        # Put chastity code here
        #print fastq_read.identifier
        seq = fastq_read.get_sequence()
        #print '1',seq
        #print '2',fastq_read
        #print '3',fastq_read.get_decimal_quality_scores()
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
                
        ##### END CHASTITY ###################
        ############################################################################################
        # Btails code
        quality_list = fastq_read.get_decimal_quality_scores()
        for trim_end in trim_ends:
            
            
            if trim_end == '5':
                lwindow_position = 0 #left position of window
                while True:
                    if lwindow_position >= len( quality_list ):
                        fastq_read.sequence = ''
                        fastq_read.quality = ''
                        break
                    if exclude_and_compare( action, quality_list[ lwindow_position:lwindow_position + window_size ], score_comparison, quality_score, exclude_window_indexes ):
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
                    if exclude_and_compare( action, quality_list[ lwindow_position:rwindow_position ], score_comparison, quality_score, exclude_window_indexes ):
                        fastq_read = fastq_read.slice( None, rwindow_position )
                        break
                    rwindow_position -= window_step
                    
        ############################################################################################
        # put  length/trim/clip code here
        #quality_list = fastq_read.get_decimal_quality_scores()
        if filter_length:
            if len(quality_list) < filter_length:
                print 'failed length'
                if failed_fastq:
                    fail.write( fastq_read )
                continue
        # Trim initial bases -- remove first 10 bases from read 2   
        if clip_length:
            # remove from the front:
            #seq = seq[clip_length:]
            fastq_read = fastq_read.slice( clip_length, None )
            #record.seq = seq
            #record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][clip_length:]
            count_of_trimmed += 1
            
        # Trim to max length -- read 2 trim to 90.
        if trim_length:
            if len(quality_list) > trim_length:
                # remove from the back side:
                #seq = seq[:trim_length]
                fastq_read = fastq_read.slice( None, trim_length )
                #record.seq = seq
                #record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][:trim_length]
                count_of_trimmed += 1
        
        
        
        if keep_zero_length or len( fastq_read ):
            out.write( fastq_read )
        else:
            num_reads_excluded += 1
    out.close()
    print 'count_of_trimmed', count_of_trimmed
    print 'count_of_first50', count_of_first50
    print 'count_of_unchaste', count_of_unchaste
    print 'count_of_Ns', count_of_Ns
    if num_reads is None:
        print "No valid FASTQ reads could be processed."
    else:
        print "%i FASTQ reads were processed." % ( num_reads + 1 )
    if num_reads_excluded:
        print "%i reads of zero length were excluded from the output." % num_reads_excluded


if __name__ == "__main__": trim_by_quality()
