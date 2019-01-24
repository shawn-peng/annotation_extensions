#!/usr/bin/python

import re
from collections import defaultdict
import Bio.UniProt.GOA as goa
from Bio.UniProt.GOA import gafiterator, record_has
import numpy as np
import csv

#from gene import *



#%% Read homolog id missing
homologs = open("homologs/human_mouse.rpt")
homologs = csv.DictReader(homologs, delimiter='\t')

#no_swiss_genes = open("no_swiss_genes.txt", "w")
no_swiss_genes = []
for row in homologs:
    if (row["SWISS_PROT IDs"] == ""):
#        no_swiss_genes.write(row["EntrezGene ID"]+"\n")
        no_swiss_genes.append(row["EntrezGene ID"])
#        print(row)
        continue

#%% id mapping
query_params = {
'from':'P_ENTREZGENEID',
'to':'ACC',
'format':'tab',
'query':''
}

import urllib,urllib3

idmap_url = 'https://www.uniprot.org/uploadlists/'

query_results = {}

entrez_to_uniprot = {}

query_params['query'] = ' '.join(no_swiss_genes)

data = urllib.parse.urlencode(query_params).encode("utf-8")
request = urllib.request.Request(idmap_url, data)
contact = "" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib.request.urlopen(request)
page = response.read(200000).decode()
#print(page)
query_results = page

for line in page.split("\n"):
    if len(line) == 0:
        continue
    row = line.split("\t")
    for f in row[0].split(","):
        t = row[1]
        entrez_to_uniprot[f] = t

#%% Parse homolog file

homologid = 0
mouse_ids = []
human_ids = []

mouse_homologs = defaultdict(set)
human_homologs = defaultdict(set)

homologs = open("homologs/human_mouse.rpt")
homologs = csv.DictReader(homologs, delimiter='\t')
homolog_rel_file = open("homolog_rel.txt", "w")
for row in homologs:
    newid = row["HomoloGene ID"]
    if newid != homologid:
        for mouse_prot in mouse_ids:
            for human_prot in human_ids:
                mouse_homologs[human_prot].add(mouse_prot)
                human_homologs[mouse_prot].add(human_prot)
                homolog_rel_file.write("homolog(%s,%s).\n"%(human_prot,mouse_prot))
        homologid = newid
        mouse_ids = []
        human_ids = []

    geneid = row["SWISS_PROT IDs"]
    if geneid == "":
        entrezid = row["EntrezGene ID"]
        if entrezid not in entrez_to_uniprot:
            print(entrezid+" not found")
            continue
        else:
            geneid = entrez_to_uniprot[entrezid]
    
    if row["Common Organism Name"] == "mouse, laboratory":
        mouse_ids.append(geneid)
    elif row["Common Organism Name"] == "human":
        human_ids.append(geneid)
    else:
        print("unexpected organism: "+row["Common Organism Name"])


#%% Generate annotation relations
import urllib
import gzip
import shutil
import os.path


species_list = ['human', 'mouse']

for species in species_list:
    gzfile = 'goa_%s.gaf.gz'%species
    if not os.path.isfile(gzfile):
        urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/%s/goa_%s.gaf.gz'
                           %(species.upper(), species), gzfile)
    
    gaffile = 'goa_%s.gaf'%species
    if not os.path.isfile(gaffile):
        f_in = gzip.open(gzfile, 'rb')
        f_out = open(gaffile, 'wb')
        shutil.copyfileobj(f_in, f_out)
    
    #handle = "goa_human.gaf"
    #handle = "goa_mouse.gaf"
    #handle = "goa_rat.gaf"
    #handle = "goa_arabidopsis.gaf"
    #handle = "goa_zebrafish.gaf"
    #handle = "goa_chicken.gaf"
    #handle = "goa_cow.gaf"
    #handle = "goa_dog.gaf"
    #handle = "goa_pig.gaf"
    #handle = "goa_fly.gaf"
    #handle = "goa_worm.gaf"
    #handle = "goa_yeast.gaf"
    #handle = open(handle)
    
    genes = defaultdict(Gene)
    
    reg_targets = defaultdict(int)
    
    annext_cnt = 0
    
    handle = open(gaffile)
    for rec in gafiterator(handle):
        
        if rec["Annotation_Extension"] and rec["DB"] == "UniProtKB":
            # print rec["Annotation_Extension"]
            protid = rec["DB_Object_ID"]
            
            genes[protid].add_annotation(rec)
            annext_cnt += 1
    
    
    print("%d annotations has extensions" % annext_cnt)
    
    dbs = {"UniProtKB": []}

#%% Genes has_regulation_target
for gene_id, gene in genes.items():
#    print("gene", gene_id)
    dbs["UniProtKB"].append(gene_id)
    for ann in gene.annotations.values():
#        if "has_regulation_target" in ann.extensions:
#            print(ann.term)
        for rel, exts in ann.extensions.items():
            for ext_target in exts:
#                print(ext.target)
                target = ext_target.split(":")
                db = target[0]
                prot = target[1]
                
                if db not in dbs:
                    dbs[db] = []
                
                dbs[db].append(prot)


#%% id mapping
query_params = {
'from':'',
'to':'GENENAME',
'format':'tab',
'query':''
}

db_name_map = {
        'UniProtKB': 'ACC+ID',
}

query_results = {}

gene_id_map = {}

for db, ids in dbs.items():
    if db not in db_name_map:
        continue
    query_params['from'] = db_name_map[db]
    query_params['query'] = ' '.join(ids)
    
    data = urllib.parse.urlencode(query_params).encode("utf-8")
    request = urllib.request.Request(idmap_url, data)
    contact = "" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(request)
    page = response.read(200000).decode()
    print(page)
    query_results[db] = page
    
    for line in page.split("\n"):
        if len(line) == 0:
            continue
        f, t = line.split("\t")
        gene_id_map[f] = t

#%% Read gene info
prot_fullnames = {}
gene_info = open("uniprot_gene_info.tsv")
import csv
csvreader = csv.DictReader(gene_info, delimiter="\t")
for row in csvreader:
#    prot_fullnames[row["Entry"]] = row["Protein names"]
    prot_fullname = row["Protein names"]
    for genename in row["Gene names"].split(" "):
        prot_fullnames[genename] = prot_fullname

#%% Write gene full names to file
#output = open("output.txt", "w")
#for gene_id, gene in genes.items():
#    if gene_id not in gene_id_map:
#        print("prot "+gene_id+" has no gene name")
#        continue
#    output.write(prot_fullnames[gene_id_map[gene_id]])
#    output.write("\n")
#    output.write("Annotations:\n")
#    for term, ann in gene.annotations.items():
#        if "has_regulation_target" not in ann.extensions:
#            continue
#        output.write("\t"+term+"\n")
#        output.write("\tExtensions:\n")
#        for ext_rel, exts in ann.extensions.items():
#            if ext_rel != "has_regulation_target":
#                continue
#            output.write("\t"+ext_rel+"\n")
#            for ext in exts:
#                target = ext.target.split(":")
#                db = target[0]
#                prot = target[1]
#                if (db != "UniProtKB"):
#                    continue
#                if (prot not in gene_id_map):
#                    print("target "+prot+" has no gene name")
#                    output.write("\t\t"+prot+"\n")
#                    continue
#                genename = gene_id_map[prot]
##                output.write("\t"+ext.target+"\n")
#                if genename in prot_fullnames:
#                    output.write("\t\t"+prot+"\t"+prot_fullnames[genename]+"\n")
#                elif genename.upper() in prot_fullnames:
#                    output.write("\t\t"+prot+"\t"+prot_fullnames[genename.upper()]+"\n")
#        output.write("\n")
#    output.write("\n\n")
    
#%% Read gene-gene identity matrix
#import h5py
#import numpy as np
#filepath = 'gene_idn/Matrix_H_H_final.mat'
#arrays = {}
#f = h5py.File(filepath)
#for k, v in f.items():
#    arrays[k] = np.array(v)
#geneidn = arrays['M']

#%% Read gene identity matrix index
geneidn_index = {}
#from Bio import SeqIO
#i = 0
#regex = re.compile("sp\|(\w*)\|\w*")
#for record in SeqIO.parse("gene_idn/sprot_human.fasta", "fasta"):
#    prot = regex.match(record.id).group(1)
#    if prot is None:
#        print(record.id)
#    geneidn_index[prot] = i
#    i += 1
    
#%% Read gene ontology obo file
goterm_fullnames = {}
goobo = open("ontology_terms/go.obo")
state = 'header'
termid = None
termname = None
termnamespace = None
for line in goobo:
    line = line.rstrip()
    if state == 'header' or state == 'next':
        if line == '[Term]':
            state = 'term'
    elif state == 'term':
        if len(line) == 0:
            state = 'next'
            goterm_fullnames[termid] = "(%s) %s" % (termnamespace, termname)
        else:
            field, val = line.split(": ", 1)
            if field == "id":
                termid = val
            elif field == "name":
                termname = val
            elif field == "namespace":
                termnamespace = val


#%% Annotation with extension to same term
output = open("terms.txt", "w")
annotated_genes = {}
for gene_id, gene in genes.items():
    for term, ann in gene.annotations.items():
        if term not in annotated_genes:
            annotated_genes[term] = {}
        annotated_genes[term][gene_id] = ann

extcnt = 0
relcnts = []
for term, term_genes in annotated_genes.items():
#    if len(term_genes) <= 1:
#        continue
    if term in goterm_fullnames:
        output.write(term+" "+goterm_fullnames[term]+":\n")
    else:
        output.write(term+":\n")
    for gene_id, ann in term_genes.items():
        if gene_id not in gene_id_map:
            output.write("\t"+gene_id+":\n")
        else:
            genename = gene_id_map[gene_id]
            if genename in prot_fullnames:
                output.write("\t"+gene_id+"\t"+prot_fullnames[genename]+"\n")
            elif genename.upper() in prot_fullnames:
                output.write("\t"+gene_id+"\t"+prot_fullnames[genename.upper()]+"\n")
            else:
                output.write("\t"+gene_id+"\n")
        
        for rel, exts in ann.extensions.items():
            output.write("\t"+rel+"(\n")
            relcnts.append(len(exts))
            if len(exts) == 0:
                print(term, gene_id, rel)
#            for ext in exts:
            for ext_target in exts:
                extcnt += 1
                target = ext_target.split(":")
                db = target[0]
                objid = target[1]
                if (db == "UniProtKB"):
                    prot = objid
                    if (prot not in gene_id_map):
                        print("target "+prot+" has no gene name")
                        output.write("\t\t"+prot+"\n")
                        continue
                    genename = gene_id_map[prot]
    #                output.write("\t"+ext_target+"\n")
                    if prot in geneidn_index and gene_id in geneidn_index:
                        i1 = geneidn_index[prot]
                        i2 = geneidn_index[gene_id]
                        if i1 < i2:
                            i1, i2 = i2, i1
                        idn = geneidn[i1, i2]
                    else:
                        print(prot, gene_id)
                        idn = 0
                    if genename in prot_fullnames:
                        output.write("\t\t"+prot+"\t"+str(idn)+"\t"+prot_fullnames[genename]+"\n")
                    elif genename.upper() in prot_fullnames:
                        output.write("\t\t"+prot+"\t"+str(idn)+"\t"+prot_fullnames[genename.upper()]+"\n")
                    else:
                        output.write("\t\t"+prot+"\t"+str(idn)+"\n")
                elif (db == "GO"):
                    termid = ext_target
                    if (ext_target in goterm_fullnames):
                        output.write("\t\t"+termid+"\t"+goterm_fullnames[termid]+"\n")
                    else:
                        output.write("\t\t"+termid+"\n")
                else:
                    output.write("\t\t"+ext_target+"\n")
                    continue
            output.write("\t)\n")
        output.write("\n")
    output.write("\n\n")
output.close()

print("number of extensions: %d" % (extcnt))

import matplotlib
import matplotlib.pyplot as plt

relcnts = np.array(relcnts, dtype=np.int32)
cnts, n, aa = plt.hist(relcnts, range(1, relcnts.max()+2))


