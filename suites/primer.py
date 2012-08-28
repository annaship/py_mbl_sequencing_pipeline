# -*- coding: utf-8 -*-
# name
# sequence
# direction
# domain
# dna_region
import pipeline.primer_utils as putils
from pipeline.utils import revcomp
from pipeline.pipelinelogging import logger


class Primer:
    """Doc string here.."""
    Name = "Primer"
    def __init__(self,name,dir,domain,region,seq):
        self.name = name
        self.direction = dir
        self.original_seq = seq
        self.domain = domain
        self.region = region
        self.primer_names = {}
        self.expanded_seqs = putils.expand(self.original_seq)
        
	    
class PrimerSuite:
    """Doc string here.."""
    Name = "PrimerSuite"
    def __init__(self, runobj, domain, region, lane_key):
        self.domain = domain
        self.region = region
        
        # list of included primer classes
        self.primer_list ={} 
        self.primer_list['F']=[]        
        self.primer_list['R']=[]
        
        # list of primer sequences
        self.primer_seq_list={}
        self.primer_seq_list['F'] = []
        self.primer_seq_list['R'] = [] 
        # list of expanded primer sequences
        self.primer_expanded_seq_list={}
        self.primer_expanded_seq_list['F'] = []
        self.primer_expanded_seq_list['R'] = []        
        self.primer_names ={}
        self.primer_names_by_reverse_complement = {}
        # changes Bacteria to Bacterial for the name
        if self.domain[-1:] != 'l': 
            self.domain = self.domain + "l"
        self.name = self.domain + ":" + self.region

        # now they can specify to not use the mbl primers for this lane/runkey
        if runobj.samples[lane_key].use_mbl_primers == 1:
            suite = runobj.primer_suites[self.name]
            if suite != None:
                for key, value in suite.items():
                    direction = value['direction']
                    sequence  = value['sequence']
                    domain    = value['domain']
                    region    = value['region']
                    p = Primer(key,direction,domain,region,sequence)
                    self.primer_list[direction].append(p)
                    self.primer_seq_list[direction].append(sequence)
                    self.primer_expanded_seq_list[direction] = self.primer_expanded_seq_list[direction] + p.expanded_seqs
                    for eseq in p.expanded_seqs:
                        self.primer_names[eseq] = key
                        # we will need this for Reverse reads
                        self.primer_names_by_reverse_complement[revcomp(eseq)] = key
                
        # they may have given some forward/reverse primers on a per lane/runkey basis
        direction = 'F'
        for idx, sequence in enumerate(runobj.samples[lane_key].forward_primers):
            key = "cust_" + lane_key + "_" + direction + "_" + str(idx)
            p = Primer(key,dir,'custom','custom',sequence)
            self.primer_list[direction].append(p)
            self.primer_seq_list[direction].append(sequence)
            self.primer_expanded_seq_list[direction] = self.primer_expanded_seq_list[direction] + p.expanded_seqs
            for eseq in p.expanded_seqs:
                self.primer_names[eseq] = key
                # we will need this for Reverse reads
                self.primer_names_by_reverse_complement[revcomp(eseq)] = key
                
        direction = 'R'
        for idx, sequence in enumerate(runobj.samples[lane_key].reverse_primers):
            key = "cust_" + lane_key + "_" + direction + "_" + str(idx)
            p = Primer(key,dir,'custom','custom',sequence)
            self.primer_list[direction].append(p)
            self.primer_seq_list[direction].append(sequence)
            self.primer_expanded_seq_list[direction] = self.primer_expanded_seq_list[direction] + p.expanded_seqs
            for eseq in p.expanded_seqs:
                self.primer_names[eseq] = key
                # we will need this for Reverse reads
                self.primer_names_by_reverse_complement[revcomp(eseq)] = key
            
            
            
            