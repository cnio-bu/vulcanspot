#!/usr/bin/python

import os
import subprocess
import pdb
import multiprocessing
import re

gene_drug = {}

filei = open ("./src/PanDrugs/PanDrugsFiles/Pandrugs_Feb2018.tsv", "r")
for line in filei:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    if line_a[0] != "gene_symbol":
        gene = line_a[0]
        if gene+":D" in gene_drug.keys():
            gene_drug[gene+":D"].append(line_a)
        else:
            gene_drug[gene+":D"] = []
            gene_drug[gene+":D"].append(line_a)

        if line_a[15] != "":
            genes_ind = line_a[15].split("|")
            for gen in genes_ind:
                if gen+":I" in gene_drug.keys():
                    gene_drug[gen+":I"].append(line_a)
                else:
                    gene_drug[gen+":I"] = []
                    gene_drug[gen+":I"].append(line_a)
filei.close()

filei = open ("./data/DrugPrescription/knowledgebased_pandrugs_input.tsv","r")
fileo = open ("./data/DrugPrescription/knowledgebased_pandrugs.tsv", "w")
for line in filei:
    line = line.rstrip("\n")
    line_a = line.split("\t")
    if len(line_a) < 9: line = line+"\t"
    if line_a[7] == "ONC" or line_a[7] == "UNCLASSIFIED":
        if line_a[2]+":D" in gene_drug.keys() or line_a[2]+":I" in gene_drug.keys():
            if line_a[2]+":D" in gene_drug.keys():
                write = 0
                for association in gene_drug[line_a[2]+":D"]:
                    if association[6] not in ["Withdrawn","Undefined"]:
                        fileo.write(line+"\t"+"\t".join([association[0],association[1],association[4],association[6],association[7],association[8],association[9],association[10],association[11],association[12],association[13],association[14],association[15],'DIRECT','',association[16]])+"\n")
                        write = 1
                if write == 0: fileo.write(line+"\n")
            if line_a[2]+":I" in gene_drug.keys():
                write = 0
                for association in gene_drug[line_a[2]+":I"]:
                    if association[6] not in ["Withdrawn","Undefined"] and association[12] == "target" and association[13] != "resistance":
                        fileo.write(line+"\t"+"\t".join([association[0],association[1],association[4],association[6],association[7],association[8],association[9],association[10],association[11],association[12],association[13],association[14],association[15],'INDIRECT',line_a[2],association[16]])+"\n")
                        write = 1
                if write == 0: fileo.write(line+"\n")
        else: fileo.write(line+"\n")
    else:
        if line_a[2]+":I" in gene_drug.keys():
            write = 0
            for association in gene_drug[line_a[2]+":I"]:
                if association[6] not in ["Withdrawn","Undefined"] and association[12] == "target" and association[13] != "resistance":
                    fileo.write(line+"\t"+"\t".join([association[0],association[1],association[4],association[6],association[7],association[8],association[9],association[10],association[11],association[12],association[13],association[14],association[15],'INDIRECT',line_a[2],association[16]])+"\n")
                    write = 1
            if write == 0: fileo.write(line+"\n")
        else: fileo.write(line+"\n")
filei.close()
fileo.close()
