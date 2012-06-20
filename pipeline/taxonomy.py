

#########################################
#
# scriptname: Taxonomy.py
#
# Author: Susan Huse, shuse@mbl.edu
#  adapted for python by Andrew Voorhis June 2012
#
# Date: May 2008
#
# Copyright (C) 2008 Marine Biological Laborotory, Woods Hole, MA
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# For a copy of the GNU General Public License, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# or visit http://www.gnu.org/copyleft/gpl.html
#
# Keywords : remove the space before the colon and list keywords separated by a space
# 
########################################

#
# 
class Taxonomy:
    """
    Create taxonomic objects
    Return classes or full text of a taxonomy object
    
    """
    def __init__(self, tax_string):
        self.tax_string = tax_string
        if not self.tax_string:
            self.tax_string = 'Unknown'
        temp = self.tax_string.strip().split(';')
        self.data = []
        for name in temp:
            # Trim leading and trailing spaces
            self.data.append(name.strip())
            
        # Remove trailing NAs and replace internal blanks with "Unassigned"
        assigned = 0
        for i in range(7,-1,-1):
            if i in self.data and i == len(self.data)-1 and ( self.data[i] == "NA" or not self.data[i] ):
                self.data.pop()
            if not assigned and i in self.data:
                assigned = 1
            if i not in self.data and assigned == 1:
                self.data[i] = "Unassigned"

        
        self.outstring = ';'.join(self.data)
        
    def taxstring(self):
        """
        # Return the object as a ";" delimited string
        """
        #return self.outstring        
        return ';'.join(self.data)
    
    def depth(self):
        """
        # Return the depth of an object - last rank with valid taxonomy
        """
        depth = "NA"
        ranks = ("domain", "phylum", "class", "order", "family", "genus", "species", "strain" )
        for i in range(0,len(self.data)):
            if i < len(self.data) and self.data[i] != "NA" and self.data[i] != "" and self.data[i] != "Unassigned":
                depth = ranks[i]
                
        return depth
        
    def domain(self):
        """
        # Return the domain of an object
        """
        return self.data[0]
        
    def phylum(self):
        """
        # Return the phylum of an object
        """
        return self.data[1]
        
    def class_name(self):
        """
        # Return the class of an object
        """
        return self.data[2]
        
    def order(self):
        """
        # Return the order of an object
        """
        return self.data[3]
        
    def family(self):
        """
        # Return the family of an object
        """
        return self.data[4]
        
    def genus(self):
        """
        # Return the genus of an object
        """
        return self.data[5]
        
    def strain(self):
        """
        # Return the strain of an object
        """
        return self.data[6]
        
    
        
# consensus is not in the Class Taxonomy    
def consensus(taxObjects, majority):
    """
    Calculate consensus of an array of taxonomic objects
    # For an array of tax objects and a majority required, calculate a consensus taxonomy
    # Return the consensus tax object, as well as stats on the agreement
    """
    #my $self = shift;


    # Correct for percentages 1-100
    if majority <= 1:
        majority = majority * 100

    # Set up variables to store the results
    newTax = []; # consensus taxon
    rankCounts = [] # number of different taxa for each rank
    maxPcts = []    # percentage of most popular taxon for each rank
    naPcts = []     # percentage of each rank that has no taxonomy assigned
    conVote=0
    taxCount = len(taxObjects)
    minRankIndex = -1
    minRank = "NA"
    ranks = ("domain", "phylum", "class", "order", "family", "genus", "species", "strain" )

    # 
    # Calculate the Consensus
    #

    # Flesh out the taxonomies so they all have indices to 7
    
    for t in taxObjects:
        for i in range(0,7):
            if len(t.data) - 1 < i:
                # If no value for that depth, add it
                t.data.append("NA")
        
    done = False
    # For each taxonomic rank
    for i in range(0,7):
    
        # Initializes hashes with the counts of each tax assignment
        tallies = {} # for each tax value -- how many objects have this taxonomy
        rankCnt=0 # How many different taxa values are there for that rank
        maxCnt=0 # what was the size of the most common taxon
        naCnt=0 # how many are unassigned 
        topPct=0 # used to determine if we are done with the taxonomy or not

        # Step through the taxonomies and count them
        for t in taxObjects:
            if t.data[i] in tallies:
                tallies[t.data[i]] += 1
            else:
                tallies[t.data[i]] = 1

        # For each unique tax assignment
        for k in tallies.iterkeys():
            if k != "NA":
                rankCnt += 1
                minRankIndex = i 
                if tallies[k] > maxCnt:
                    maxCnt = tallies[k]
            else:
                naCnt = tallies[k];

            #print 'in newtax1',k,tallies[k],taxCount
            vote = int(( float(tallies[k]) / float(taxCount) ) * 100)
            #print 'in newtax2',k,vote,majority
            if k != "NA" and vote > topPct:
                topPct = vote
            
            if not done and vote >= majority:
                
                newTax.append(k) 
                if k != "NA":
                    conVote = vote 
            
        if topPct < majority:
            done = True

        #if ($#newTax < $i) {push (@newTax, "NA");}
        rankCounts.append(rankCnt)
        if taxCount > 0:
            maxPcts.append(int(100 * ( float(maxCnt) / float(taxCount) ) + 0.5))
            naPcts.append(int(100 * ( float(naCnt) / float(taxCount) ) + 0.5))


    taxReturn=[]
    # 0=taxObj, 1=winning vote, 2=minrank, 3=rankCounts, 4=maxPcts, 5=naPcts;
    if len(newTax) == 0:
        taxReturn.append( Taxonomy("Unknown") ) #If no consensus at all, call it Unknown
    else:
        taxReturn.append( Taxonomy(';'.join(newTax)) ) # taxonomy object for consensus
    
    #if (! $taxReturn[0]) {$taxReturn[0] = "NA";}
    if not taxReturn[0]:
        taxReturn[0] = "Unknown" # 20081126 - empty tax should be 'Unknown'
    #if taxReturn[-1] == "Unassigned":
    #    taxReturn.pop() # If resolved to an Unassigned rank, remove it.
    taxReturn.append(str(int(conVote + 0.5))) # winning majority
    if minRankIndex >= 0:
        minRank = ranks[minRankIndex]

    taxReturn.append(minRank) # lowest rank with valid assignment
    taxReturn.append(';'.join([str(i) for i in rankCounts])) # number of different taxa at each rank
    taxReturn.append(';'.join([str(i) for i in maxPcts])) # percentage of the most popular taxon (!= "NA")
    taxReturn.append(';'.join([str(i) for i in naPcts])) # percentage that are unassigned ("NA")
    return taxReturn

