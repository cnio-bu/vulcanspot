#!/usr/bin/python

import os
import re
import subprocess
import glob
import pdb

#GScore
tumorportal = {}
cgc = {}
driver = {}
es = {}
os = {}

TUMORPORTAL = open("./src/PanDrugs/databases_Feb_2018_illumina/TumorPortal.csv","r")
for line in TUMORPORTAL:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    sco = 0
    if line_a[2] == "Near significance": sco = 0.025
    if line_a[2] == "Significantly mutated": sco = 0.05
    if line_a[2] == "Highly significantly mutated": sco = 0.1

    if line_a[0] not in tumorportal.keys():
        tumorportal[line_a[0]] = sco
    else:
        if sco > tumorportal[line_a[0]]:
            tumorportal[line_a[0]] = sco
TUMORPORTAL.close()

CGC = open("./src/PanDrugs/databases_Nov_2016/CGC.tsv","r")
for line in CGC:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    cgc[line_a[0]] = 0.1
CGC.close()

DRIVER = open("./src/PanDrugs/databases_Feb_2018_illumina/srep02650-s3.csv","r")
for line in DRIVER:
    line = line.rstrip()
    line_a = line.split(",")
    sco = 0
    if line_a[7] == "High Confidence Driver": sco = 0.1
    if line_a[7] == "Candidate driver": sco = 0.05
    driver[line_a[0]] = sco
DRIVER.close()

ES = open("./src/PanDrugs/databases_Feb_2018_illumina/gene_essentiality_score.tsv","r")
for line in ES:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    if line_a[1] != "max_min_score":
        es[line_a[0]] = 0.4 * float(line_a[1])
ES.close()

OS = open("./src/PanDrugs/PanDrugsFiles/oncoscape_all_matrix_highscore.tsv","r")
for line in OS:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    if line_a[0] != "SYMBOL":
        os[line_a[0]] = 0.3 * float(line_a[36])/4
OS.close()

filei = open("./data/CCL/genes.txt","r")
fileo = open("./data/CCL/gscores.tsv","w")
gene_list = []
for line in filei:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    geneA = line_a[0]

    if line_a[0] == "GeneA":
        fileo.write("\t".join(["Gene","GScore"])+"\n")
    else:
        if geneA not in gene_list:
            gene_list.append(geneA)
            gscore = 0
            if geneA in tumorportal.keys(): gscore += tumorportal[geneA]
            if geneA in cgc.keys(): gscore += cgc[geneA]
            if geneA in driver.keys(): gscore += driver[geneA]
            if geneA in es.keys(): gscore += es[geneA]
            if geneA in os.keys(): gscore += os[geneA]

            fileo.write("\t".join([geneA,str(round(gscore,4))])+"\n")

filei.close()
fileo.close()
