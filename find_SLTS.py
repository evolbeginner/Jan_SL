#! /usr/local/bin python

import getopt
import re
import sys
import os
import sqlite3

sys.path.append(os.path.join(os.path.dirname(__file__),'lib'))
import gff_class
from selfbasics import make_dir_force


####################################################################
def read_gff(ingff,types_included=['gene'],target_regexps=['ID']):
    gene_info={}
    gff_objs = {}
    with open(ingff,'r') as f:
        for line_no,line in enumerate(f.readlines()):
            if line.startswith('#'):
                continue
            gff_item_obj = gff_class.gff_items_class(line,types_included,target_regexps)
            gff_item = gff_item_obj.parse_gff_line()
            if gff_item:
                for key,value in gff_item.iteritems():
                    value.chr = re.sub('\|.+',"",value.chr)
                gff_objs.update(gff_item)
    f.close()
    return(gff_objs)
    
    
def generate_indb(gff_objs):
    indb = "./genes_gff.sqlite3"
    os.system("rm -rf " + indb)
    cx = sqlite3.connect(indb)
    cu = cx.cursor()
    create_table_sqlite3 = 'create table genes_gff_info (id character primary key,chr character,strand character,start int,stop int)'
    cu.execute(create_table_sqlite3)
    for gene,obj in gff_objs.iteritems():
        chr = obj.chr
        strand = obj.strand
        start  = obj.start
        stop   = obj.stop
        gene   = add_quotes(gene)
        chr = add_quotes(chr)
        strand = add_quotes(strand)
        execute_sqlite3 = "insert into genes_gff_info values(" + \
                            ','.join([gene,chr,strand,start,stop]) + ')'
        cu.execute(execute_sqlite3)
    cx.commit()
    cu.close()
    cx.close()
    return(indb)


def search_db(db, gff_file, upstream): # SL_type: with_SL or without_SL
    cx = sqlite3.connect(db)
    cu = cx.cursor()
    genes_list={}
    with open(gff_file,'r') as f:
        for line_no,line in enumerate(f.readlines()):
            if line.startswith('#'):
                continue
            line = line.rstrip('\n\r')
            gff_item_obj = gff_class.gff_items_class(line,['primer'],['ID'])
            obj = gff_item_obj.parse_gff_line()
            (gene,) = obj.keys()
            chr     = add_quotes(obj[gene].chr)
            start   = int(obj[gene].start)
            stop    = int(obj[gene].stop)
            if obj[gene].strand == "-":
                # start and stop will be changed
                start, stop = stop, start
                if upstream:
                    start += upstream
            else:
                if upstream:
                    start -= upstream
            start = str(start)
            stop  = str(stop)

            execute_sqlite3s = []
            execute_sqlite3 = "select id from genes_gff_info where start BETWEEN " + start + " AND " + stop + \
                                " and " + "chr= " + chr
            execute_sqlite3s.append(execute_sqlite3)
            execute_sqlite3 = "select id from genes_gff_info where stop BETWEEN " + start + " AND " + stop + \
                                " and " + "chr= " + chr
            execute_sqlite3s.append(execute_sqlite3)
            execute_sqlite3 = "select id from genes_gff_info where start < " + start + " and stop > " + stop + \
                                " and " + "chr= " + chr
            execute_sqlite3s.append(execute_sqlite3)

            for execute_sqlite3 in execute_sqlite3s:
                get_sqlite3_searching_result(execute_sqlite3,cu)
                fetch_all_result_list = cu.fetchall()
                if fetch_all_result_list:
                    genes_list.update({fetch_all_result_list[0][0]:1})
                    break
            #if line_no+1 >= 500:   break
    cu.close()
    cx.close()
    return(genes_list)


def get_sqlite3_searching_result(execute_sqlite3,cu):
    cu.execute(execute_sqlite3) 
    return (cu)


def add_quotes(item):
    new_item = item.join(["'"]*2)
    return(new_item)


def read_experimental_SLTS(indir,types_included,target_regexps=['ID']):
    files={}
    gff_SLTS_objs={}
    for file in os.listdir(indir):
        if re.search('noSL',file):
            file['without_SL']=file
        else:
            file['with_SL']=file
    for type_of_file,file in files.iteritems():
        gff_SLTS_objs[type_of_file]={}
        with open(file,'r') as f:
            if line.startwith('#'):
                continue
            gff_item_obj = gff_class.gff_items_class(line,types_included=types_included,target_regexps=target_regexps)
            gff_item = gff_item_obj.parse_gff_line()
            if gff_item:
                gff_SLTS_objs.update(gff_item)
        f.close()
    return(gff_SLTS_objs)


def show_help():
    print "\nUsage:"
    print "python " + os.path.basename(__file__)
    ["indb=","outdb=","outdir=","ingff=","types=","no_isoform=","force","experimental_SLTS_dir="],
    print '''Required arguments:
    --ingff1    gff file used to build seqlite3 database
    --ingff2    gff file used as query 
    --indb      input database (sqlite3)
Optional arguments:
    --outdir    default: ./
    --types     types
                defalt: gene
    --upstream  number of bps upstream of a gene
    --force     if outdir has already existed, remove it!
    '''
    sys.exit()


####################################################################
indb=None
outdb=None
outdir="./"
ingff1s=[]
ingff2s=[]
types_included=['gene']
experimental_SLTS_dir="/home/sswang/huche/help/Jan/SL/data/genomic_data/Smin/original_data/real_SL_containing_genes_Smin/"
force=False
genes_list={}
upstream=None

try:
    opts, args = getopt.getopt(
        sys.argv[1:],
        "hl:i:f:o:",
        ["indb=","outdb=","outdir=","ingff=","types=","no_isoform=","force","experimental_SLTS_dir=","upstream="],
    )
except getopt.GetoptError:
    print "illegal error!"
    show_help()
for op, value in opts:
    if op == '--indb':
        indb=os.path.expanduser(value)
    elif op == '--outdb':
        outdb=os.path.expanduser(value)
    elif op == '--outdir':
        outdir=os.path.expanduser(value)
    elif op == '--ingff1':
        for i in value.split(','):
            ingff1s.append(os.path.expanduser(value))
    elif op == '--ingff2':
        for i in value.split(','):
            ingff2s.append(os.path.expanduser(value))
    elif op == "--types":
        for i in value.split(','):
            types_included.append(i)
    elif op == "experimental_SLTS_dir":
        experimental_SLTS_dir=os.path.expanduser(value)
    elif op == "--no_isoform":
        pass
    elif op == "--upstream":
        upstream = int(value)
    elif op == "--force":
        force=True


####################################################################
make_dir_force(outdir,force)

if not indb:
    if ingff1s:
        for ingff1 in ingff1s:
            gene_gff_objs = read_gff(ingff1,types_included)
            indb = generate_indb(gene_gff_objs)
    else:
        print "Either indb or ingff has to be given!"
        show_help()

# read_experimental_SLTS(experimental_SLTS_dir,['primer'],['ID'])
# for SL analyses
if not ingff2s:
    for file in os.listdir(experimental_SLTS_dir):
        ingff2s.append(os.path.join(experimental_SLTS_dir,file))

for index,ingff2 in enumerate(ingff2s):
    print ingff2
    basename = os.path.basename(ingff2)
    genes_list[basename] = search_db(indb, ingff2, upstream)

for key,value in genes_list.iteritems():
    print key
    try:
        f = open(os.path.join(outdir,key+'.sqlite3_out'),'w')
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        print "next!"
        break
    for v in value.keys():
        f.write(v+'\n')
    f.close()


