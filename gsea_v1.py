### use python gsea_v1.py ${SPECIE} ${PATH_TO_SET_GENES}
import os,sys
import numpy as np
import scipy.stats as stats 
#import rpy2.robjects as robjects
#from rpy2.robjects import r
#from rpy2.robjects.packages import importr

### hg analysis
try :
    specie = sys.argv[1]
except IndexError:
    print "please provide a specie, arabidopsis or maize"
    sys.exit()

### step 1. make association file gene_goterm
desc_go = {}
desc = open("annotation/go_desc_"+specie+".txt", 'r')
desc_lines = desc.readlines()
for desc_line in desc_lines :
    desc_go[desc_line.split("\t")[0]] = desc_line.split("\t")[1].strip("\n")
pop_list = []
for gene in open("annotation/genes_"+specie+".txt", 'r') :
    pop_list.append(gene.strip("\n"))
pop_list = list(set(pop_list))
pop_count = len(pop_list)
### seed set 
set_list = []
for set_gene in open(sys.argv[2], 'r') :
    set_list.append(set_gene.strip(" \n"))
set_list = list(set(set_list))
set_count = len(set_list)

### gene: go dict
assoc_dict = dict()
assoc = open("annotation/assoc_"+specie+".txt", 'r') 
assoc_lines = assoc.readlines() 
for line in assoc_lines :
    for go in line.split("\t")[1].split(";") :
        if len(go) > 1 :
            assoc_dict.setdefault(line.split("\t")[0], []).append(go.strip("\n"))

set_go_list = list()
for gene in assoc_dict.keys() :
    if gene in set_list and assoc_dict[gene] not in set_go_list : 
        set_go_list.extend(assoc_dict[gene])
set_go_list = list(set(set_go_list))

### go: gene dict
set_go_dict = dict()
for go_term in set_go_list :
    for gene in assoc_dict.keys() :
        if go_term in assoc_dict[gene] :
            if gene not in set_go_dict.setdefault(go_term, []) :
                set_go_dict.setdefault(go_term, []).append(gene)
outfile = open(sys.argv[2].split(".")[0]+"_hg.out", 'w')
outfile.write("GO term\tDescription\tOverlap Count\tGO count\tPopulation\tSet Count\tP-value\tAdj. P-value\n")
for go_term in set_go_dict.keys() :
    go_count = 0
    for gene_name in set_list :
        if gene_name in set_go_dict[go_term] :
            go_count += 1
    m = len(set_go_dict[go_term]) ### total no. of genes in this go term 
    n = pop_count ### total no. of genes
    k = set_count ### total no. of set genes
    x = go_count  ### no. of set genes with this go term
    if go_count > 1 :
        print go_term, go_count-1, len(set_go_dict[go_term]), n-(len(set_go_dict[go_term])), k
        pval= stats.hypergeom.sf(x,n,m,k)
        adjpval = pval*x
        gocount = len(set_go_dict[go_term])
        pop = n-len(set_go_dict[go_term])
        outfile.write(go_term)
        outfile.write("\t")
        outfile.write(desc_go[go_term])
        outfile.write("\t")
        outfile.write(str(x))
        outfile.write("\t")
        outfile.write(str(m))
        outfile.write("\t")
        outfile.write(str(n))
        outfile.write("\t")
        outfile.write(str(k))
        outfile.write("\t")
        outfile.write(str(float(str(pval))))
        outfile.write("\t")
        outfile.write(str(float(str(adjpval))))
        outfile.write("\n")
outfile.close()


