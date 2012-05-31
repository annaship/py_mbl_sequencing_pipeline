
import os, sys
import pipeline.anchortrimming_mbl as anchortrim
from pipeline.pipelinelogging import logger
import logging

#from pipeline.utils import revcomp

import constants as C
try:
    import Levenshtein
except:
    print '''
    You need Levenshtein module installed to run this software.

    Here is a fast implementation of Levenshtein distance for Python:

        http://code.google.com/p/pylevenshtein/

'''
    sys.exit(-1)

def count_keys(hash):
    return len([x for x in hash])
        
def trim_stop_seq( stop_seqs, seq, trim_type, start, end ):
    
    for anchor in stop_seqs:
        anchor_length = len(anchor)
        logger.debug( anchor + " " + str(start) + " " + str(end)  + " " + str(len(seq)))
        for pos in range(start,end):
            seq_window = seq[pos:pos+anchor_length]

            
            dist = abs( Levenshtein.ratio( anchor,     seq_window ) )
            #dist2 = abs( Levenshtein.ratio( seq_window, anchor )     )
            if dist == 1.0:
                # perfect match
                # do I trim off before or after anchor?
                #before
                return anchor,seq[pos+anchor_length:],seq[:pos]
                #after (include anchor in trimmed seq):
                #return anchor,seq[pos:],seq[:pos+anchor_length]
            if dist >= C.max_divergence:
                pass
            #print pos,seq_window,dist1,dist2
    
    return '','',seq

# the anchor strings were already reversed back at the start of trim_run    
# this now just calls the real worker so that we can more easily test the trim_anchor_helper()
# without having to create anchor data structures in the test code
def trim_anchor(anchor_name, expanded_anchor_sequences, anchor_definition, sequence):
    return trim_anchor_helper(anchor_name, expanded_anchor_sequences, anchor_definition['freedom'], anchor_definition['length'], anchor_definition['start'], sequence)

# The real worker function    
def trim_anchor_helper(anchor_name, expanded_anchor_sequences, freedom, length, start, sequence):
    exact = ''
    exactTrimmedOff = ''    

    logger.debug( 'looking for anchor: ' + anchor_name + " start: " + str(start) + " length: " + str(length))
    max_divergence  = C.max_divergence
    logger.debug('anchor_list: ' + str(expanded_anchor_sequences))
    list_of_tuples = anchortrim.generate_tuples(start, freedom, length, list_of_tuples = [], reversed_read=False)
    logger.debug('anchor tuples: ' + str(list_of_tuples))
    anchor, location = anchortrim.find_best_distance(sequence, expanded_anchor_sequences, max_divergence, list_of_tuples)
    
    if anchor and location:
        logger.debug( 'anchor: ' + anchor + ' loc tuple: ' + str(location))
        trimmed_sequence = sequence[:location[1]] # same thing here for the reversed == False
        exact = anchor
        exactTrimmedOff = sequence[location[1]:]
    else:
        logger.debug( 'no anchor location found' )
        trimmed_sequence = sequence
    return exact, exactTrimmedOff, trimmed_sequence

def trim_distal_primer(primers_list, seq):
    """Doc string here"""
    d_primer_found = ""
    trimmedPortion = ""
    loc = 0
    seq_len = len(seq)
    for p in primers_list:
        (primer_found, trimmedPortion, ret_seq) = do_actual_distal_trim(p, seq, allow_offset = True)
        if primer_found:
            return primer_found, trimmedPortion, ret_seq
        else:
            truncLength = len(p) + 1 #add one just in case!
            while truncLength >= 5:
                short_p = p[:truncLength]
                (primer_found, trimmedPortion, ret_seq) = do_actual_distal_trim(short_p, seq, allow_offset = False)
                if primer_found:
                    return primer_found, trimmedPortion, ret_seq
                # prepare for next try                    
                truncLength = truncLength - 1
    return '', '', seq

# do the real attempt to do the distal trim
def do_actual_distal_trim(p, seq, allow_offset):
    seq_len = len(seq)
    # find the furthest RIGHT p in seq
    p_index = seq.rfind(p)
    if p_index != -1:
        # found whole exact primer
        primer_found = p
        # keep EVERYTHING we trimmed off 
        trimmedPortion = seq[p_index:] 
        len_of_trimmed = len(trimmedPortion)
        # was this found at exactly the end of the sequence OR
        # if we have more than just the primer
        # then is the 'extra' after the primer less or equal than
        # what we allow
        if len_of_trimmed == len(p) or (allow_offset and (seq_len-(p_index+len(p)) <= C.distal_from_end)):
            seq = seq[:p_index]            
            return primer_found, trimmedPortion, seq
        else:
            seq = seq[:p_index]            
            return primer_found, trimmedPortion, seq
    return ('','','')
  
def trim_fuzzy_distal(anchors_list, seq, trim_type, start, end):
    """Doc string here.."""
    max_distance = 3
    best_distance = max_distance + 1
    found_fuzzy = 0
    fuzzy_match = ""
    for anchor in anchors_list:
        anchor_length = len(anchor)
        for pos in range(start,end):
        
            seq_window = seq[pos:anchor_length]
            
            dist = 0
            
            #dist1 = abs( Levenshtein.ratio( anchor,     seq_window ) )
           # dist2 = abs( Levenshtein.ratio( seq_window, anchor )     )
            dist1 = abs( levenshtein( anchor,     seq_window ) )
            dist2 = abs( levenshtein( seq_window, anchor )     )
            if dist1 >= dist2:  dist = dist1
            else:               dist = dist2
            
            if (dist <= max_distance) and (dist < best_distance) and (seq_window[:2] == anchor[:2]):
                if seq_window[-3:] != anchor[-3:]:
                    
                    # check for deletion
                    if(seq_window[-4:][:3] == anchor[-3:]):                        
                        seq_window.strip()
                        logger.debug( "Fuzzy with deletion " + seq_window)
                    # check for insertion
                    elif(seq_window[-3:] == anchor[-4:][:3]):                        
                        seq_window = seq_window + anchor[-1:]
                        logger.debug("fuzzy with insertion " + seq_window)
                        
                # Found a fuzzy match within tolerances, so store it
                found_fuzzy = 1;
                best_distance = dist;
                best_position = pos;
                fuzzy_match = seq_window;        
                if dist == 0: 
                    found_exact = 1
                    break
    fuzz_right = ''
    if found_fuzzy:
        fuzzy_right = seq
        if( trim_type == 'internal'):
            seq = seq[:best_position + len(fuzzy_match)]    
        else:
            seq = seq[:best_position]
            
        fuzzy_right = fuzz_right[len(seq):]
    return fuzz_right, best_distance, seq, fuzzy_match
        
def levenshtein(s1, s2):
    """ Use the Levenshtein.py module instead """
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if not s1:
        return len(s2)
 
    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
 
    return previous_row[-1]

def trim_proximal_primer(primer_list, seq):
    """Doc string here.."""
    p_primer_found = "";
    offset         = 0;
    for primer in primer_list:
        #print primer
        prim_len = len(primer)
        idx = seq.find(primer)
        if(idx != -1):
            if(idx == 0):
                p_primer_found = primer
                seq = seq[prim_len:]
                break
            elif(idx == 1):
                # sometimes get a stray base between runkey and proximal primer
                # kind of close??
                p_primer_found = primer
                offset = idx
                seq = seq[idx + prim_len:]
                break
    return p_primer_found, offset, seq
   
    
def expand(seq):
    """Takes a single ambiguous primer sequence and expands it to remove ambiguities.
    Adds them to a hash (prevents dupes) and returns the hash.
    """
    expandedPrims={}  # fully clean distals - hash to prevent dupes
    workingPrimers=[]    # still cleaning holder
    bases = ['A','C','G','T']
    
    seq = seq.replace('.','N')
        #print '1',primer
    workingPrimers.append(seq.upper())
    
  
    while (len(workingPrimers) > 0):
        # this is the only place items are removed from working primers
        d = workingPrimers.pop() 
        
        # remove all commas
        # don't think we needed this d = d.replace(',','')
        
#     # For each N (blast doesn't like them) expand to 4 distals, one for each base
        if ('N' in d ):    
            for b in bases:      
                workingPrimers.append(d.replace('N',b,1))

#     # For R, Y, W, S, M, K expand to the pair of bases
        elif ('R' in d):
            workingPrimers.append(d.replace('R','A',1))
            workingPrimers.append(d.replace('R','G',1))
            #print '1',d.replace('R','G',1)

        elif ('Y' in d):
            workingPrimers.append(d.replace('Y','C',1))
            workingPrimers.append(d.replace('Y','T',1))        

        elif ('W' in d):
            workingPrimers.append(d.replace('W','A',1))
            workingPrimers.append(d.replace('W','T',1))     

        elif ('S' in d):
            workingPrimers.append(d.replace('S','G',1))
            workingPrimers.append(d.replace('S','C',1))     

        elif ('M' in d):
            workingPrimers.append(d.replace('M','A',1))
            workingPrimers.append(d.replace('M','C',1))     

        elif ('K' in d):
            workingPrimers.append(d.replace('K','G',1))
            workingPrimers.append(d.replace('K','T',1))     

#     # For each [CT] or [AGT] ... expand to each base
        elif ('[' in d):
            
            if( d.find('[') + 3 == d.find(']') ):
                # ie [AC]
                base1 = d[d.find('[')+1:d.find('[')+2]
                base2 = d[d.find('[')+2:d.find('[')+3]
                replace = '['+base1+base2+']'
                workingPrimers.append(d.replace(replace,base1,1))
                workingPrimers.append(d.replace(replace,base2,1)) 
            
            if( d.find('[') + 4 == d.find(']') ):
                # ie: [CTA]
                base1 = d[d.find('[')+1:d.find('[')+2]
                base2 = d[d.find('[')+2:d.find('[')+3]
                base3 = d[d.find('[')+3:d.find('[')+4]
                replace = '['+base1+base2+base3+']'
                workingPrimers.append(d.replace(replace,base1,1))
                workingPrimers.append(d.replace(replace,base2,1))       
                workingPrimers.append(d.replace(replace,base3,1))  

#     # expand \?
        elif ('?' in d):
            if(d.find('?') == 0):
                workingPrimers.append(d[1:])
            else:
                preceder_plus = d[d.find('?')-1:d.find('?')+1] 
                # the preceing base exists: remove '?' only
                workingPrimers.append(d.replace('?','',1))
                # preceding base doesn't exist: 'remove '?' and preceding base
                workingPrimers.append(d.replace(preceder_plus,'',1))

#     # next expand + to 1 or 2
        elif ('+' in d):
            preceder = d[d.find('+')-1:d.find('+')]
            # the preceing base exists once: remove '+' only
            workingPrimers.append(d.replace('+','',1))
            # the preceing base exists twice: change '+' to preceding base
            workingPrimers.append(d.replace('+',preceder,1))

#     # expand * to 0,1,2
        elif ('*' in d):
            if(d.find('*') == 0):
                workingPrimers.append(d[1:])
            else:
                preceder = d[d.find('*')-1:d.find('*')]
                preceder_plus = d[d.find('*')-1:d.find('*')+1] 
                #print preceder,preceder_plus
                # the preceding base doesn't exist': remove '*' and preceding base
                workingPrimers.append(d.replace(preceder_plus,'',1))
                # the preceding base exists once: remove '*' only
                workingPrimers.append(d.replace('*','',1))
                # the preceding base exists twice: change '*' to preceding base
                workingPrimers.append(d.replace('*',preceder,1))

#     # For repeat bases, e.g., C{5,8} becomes 4 new primers: Cx5, Cx6, Cx7, Cx8
        elif('{' in d):
            if(d.find('}') == d.find('{') + 4 ):
                repeatBase = d[d.find('{')-1:d.find('{')]
                toReplace = d[d.find('{')-1:d.find('}')+1]
                minCount = d[ d.find('{') + 1:d.find('{') + 2 ]
                maxCount = d[ d.find('{') + 3:d.find('{') + 4 ]
                for i in range (int(minCount),int(maxCount) + 1 ):
                    homopolymer = i * repeatBase
                    workingPrimers.append(d.replace(toReplace, homopolymer))
                
        # If it made it through everything else, move it to the final set
        # Use hash so can filter out duplicates
        else:
            if(d):
                expandedPrims[d] = 1

    return expandedPrims.keys()  

    
def get_anchor_list(run, anchor_name, adtnl_anchors_list):
    raw_sequence = run.anchors[anchor_name]['sequence']
    expanded_seqs = expand(raw_sequence)

    for a in adtnl_anchors_list:
        expanded_seqs = expanded_seqs + expand(a)

    return expanded_seqs


if __name__=='__main__':
    import unittest
    from utils import revcomp
    
    class TestAll(unittest.TestCase):
    
        # utility to replace 'x' with new_x and 'y' with new_y
        def replace_template_xy(self, orig_str_array, new_x, new_y):
            return  [ty.replace("y", new_y) for ty in [tx.replace("x", new_x) for tx in orig_str_array]]
        
        def test_expand(self):
            #test expand of .
            self.assertEqual(set(expand("A.A")), set(["AAA", "ACA", "AGA", "ATA"]))

            #test expand of N
            self.assertEqual(set(expand("ANA")), set(["AAA", "ACA", "AGA", "ATA"]))

            # test expand of N and R
            self.simple_replace_test("R", "A", "G")
            self.simple_replace_test("Y", "C", "T")
            self.simple_replace_test("W", "A", "T")
            self.simple_replace_test("S", "C", "G")
            self.simple_replace_test("M", "C", "A")
            self.simple_replace_test("K", "G", "T")
            
            # test group stuff
            # "N[CT]" => "AC", "CC", "GC", "TC", "AT", "CT", "GT", "TT"
            self.assertEqual(set(expand("N[CT]")),set(["AC", "CC", "GC", "TC", "AT", "CT", "GT", "TT"]))
            # "N[CTA]" => "AC", "CC", "GC", "TC", "AT", "CT", "GT", "TT"   , "AA", "CA", "GA", "TA"
            self.assertEqual(set(expand("N[CTA]")),set(["AC", "CC", "GC", "TC", "AT", "CT", "GT", "TT"   , "AA", "CA", "GA", "TA"]))
            
            # test optional ?
            self.assertEqual(set(expand("NC?")),set(["AC", "CC", "GC", "TC", "A", "C", "G", "T"]))
            
            # test optional +
            self.assertEqual(set(expand("NC+")),set(["ACC", "CCC", "GCC", "TCC", "AC", "CC", "GC", "TC"]))

            # test optional *
            self.assertEqual(set(expand("NC*")),set(["ACC", "CCC", "GCC", "TCC", "AC", "CC", "GC", "TC",  "A", "C", "G", "T"]))

            # test optional {1,2}
            self.assertEqual(set(expand("C{1,2}N")),set(["CCA", "CCC", "CCG", "CCT",      "CA", "CC", "CG", "CT"]))
        
        def test_trim_proximal(self):
            # test not found
            seq = "AATAT"
            exactLeft, offset, trimmed_seq =  trim_proximal_primer(["TAT"], seq)
            self.assertEqual((exactLeft, offset, trimmed_seq), ('', 0,  seq))

            # test found at 0
            seq = "TATAAA"
            exactLeft, offset, trimmed_seq =  trim_proximal_primer(["TAT"], seq)
            self.assertEqual((exactLeft, offset, trimmed_seq), ('TAT', 0,  "AAA"))

            # test found at 1
            seq = "ATATAAA"
            exactLeft, offset, trimmed_seq =  trim_proximal_primer(["TAT"], seq)
            self.assertEqual((exactLeft, offset, trimmed_seq), ('TAT', 1,  "AAA"))

        def test_trim_distal(self):
            # test not found at all
            seq = "AATAT"
            exactRight, snipped, trimmed_seq =  trim_distal_primer(["GTAT"], seq)
            self.assertEqual((exactRight, snipped, trimmed_seq), ('', '',  seq))

            # test found at end exactly
            seq = "TATAAA"
            exactRight, snipped, trimmed_seq =  trim_distal_primer(["AAA"], seq)
            self.assertEqual((exactRight, snipped, trimmed_seq), ('AAA', 'AAA',  "TAT"))

            # test found at 1 back from end
            seq = "ATATAAAG"
            exactRight, snipped, trimmed_seq =  trim_distal_primer(["AAA"], seq)
            self.assertEqual((exactRight, snipped, trimmed_seq), ('AAA', 'AAAG',  "ATAT"))
            
            # test found shrinking distal find CACACA rather than CACACACA
            seq = "ATATCACACA"
            exactRight, snipped, trimmed_seq =  trim_distal_primer(["CACACACA"], seq)
            self.assertEqual((exactRight, snipped, trimmed_seq), ('CACACA', 'CACACA',  "ATAT"))

            # test can't find shrinking distal gets too short and fails and return default stuff
            seq = "ATATCACA"
            exactRight, snipped, trimmed_seq =  trim_distal_primer(["CACACACA"], seq)
            self.assertEqual((exactRight, snipped, trimmed_seq), ('', '',  seq))
            
        def test_trim_anchor(self):
            seq = "GGGCACACATTT"
            anchors = ["CACACA"]
            #find it in the exact spot we expect it
            exact, exactTrimmedOff, trimmed_sequence = trim_anchor_helper("Test Anchor", anchors, 0, 6, 3, seq)
            self.assertEqual((exact, exactTrimmedOff, trimmed_sequence), ("CACACA", "TTT", "GGGCACACA"))

            seq = "GGCACACATTT"
            anchors = ["CACACA"]
            #find it in an acceptable spot just before it so we give it a freedom of 1
            exact, exactTrimmedOff, trimmed_sequence = trim_anchor_helper("Test Anchor", anchors, 1, 6, 3, seq)
            self.assertEqual((exact, exactTrimmedOff, trimmed_sequence), ("CACACA", "TTT", "GGCACACA"))

            seq = "GGGGCACACATTT"
            anchors = ["CACACA"]
            #find it in an acceptable spot just AFTER it so we give it a freedom of 1
            exact, exactTrimmedOff, trimmed_sequence = trim_anchor_helper("Test Anchor", anchors, 1, 6, 3, seq)
            self.assertEqual((exact, exactTrimmedOff, trimmed_sequence), ("CACACA", "TTT", "GGGGCACACA"))
            

        def simple_replace_test(self, placeholder_base, new_x, new_y):
            result_template_strings = ["xA", "xC", "xG", "xT", "yA", "yC", "yG", "yT"]
            self.assertEqual(set(expand(placeholder_base + "N")), set(self.replace_template_xy(result_template_strings, new_x, new_y)))
            
    print "expanded: " + str(expand(revcomp("GTGAATCATCGAYTCTTTGAAC")))
    
    
    logger.setLevel(logging.DEBUG)        
    unittest.main()
    
    