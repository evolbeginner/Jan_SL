#! /bin/env python

import re
import sys
import os
import getopt

##############################################################
def show_help():
    print "python " + os.path.basename(__file__) + " <-i|--gff> [Options]"
    os.sys.exit()
    
##############################################################
features=[]
GeneIDs={}
gene_locus=None
gene_loci={}
attributes_tmp=[]

try:
    opts,args = getopt.getopt(
        sys.argv[1:],
        "hi:",
        ["gff=","feature=","help"],
    )
except getopt.GetoptError:
    print "illegal error!"
    show_help()

for op, value in opts:
    if op == '-i' or op == 'gff':
        gff = value
    elif op == '--feature':
        feature = value
    elif op == '-h' or op == '--help':
        show_help()

if not features:
    features = ['CDS']

##############################################################
for line in open(gff, 'r'):
    if line.startswith('#'):
        continue
    attributes_tmp = []
    line = line.rstrip('\n\r')
    line_array = line.split("\t")
    others, attribute = line_array[0:8], line_array[8]
    if not re.search("=GeneID",attribute):
        continue
    m = re.findall("([^;=]+)=[^;=]+", attribute)

    if others[2] == 'gene':
        m = re.search('gene=([^;=]+)',attribute)
        gene_locus = m.group(1)
        if not gene_locus in gene_loci:
            gene_loci[gene_locus] = {}
        # gene_loci[gene_locus] +=1
       
#######################################
    if not others[2] in features:
        continue

    for i in re.findall("[^;=]+=[^;=]+", attribute):
        m = re.search('([^;=]+)=([^;=]+)',i)
        if m.group(1) == 'Dbxref':
            n = re.search(r'(GeneID):([^.,;]+)',m.group(2))
            GeneID = ''.join([n.group(1),n.group(2)])
            x = re.search("Genbank:(.+)",m.group(2))
            protein_id = x.group(1)
            if not GeneID in GeneIDs:
                GeneIDs[GeneID] = {}
            GeneIDs[GeneID][protein_id] = ''
            gene_loci[gene_locus][GeneID] = ''

            attribute = 'Parent=' + '.'.join([gene_locus, str(max(len(GeneIDs[GeneID]),len(gene_loci[gene_locus])))])
            break
    print "\t".join(others) + "\t" + attribute


