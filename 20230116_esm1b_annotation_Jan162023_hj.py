#!/usr/bin/env python
# coding: utf-8

import sys
import re
import os
import pandas as pd


#if len(sys.argv) != 4:
#    print("please, check the input arguments")
#else:
#    vt = sys.argv[1]
#    esm1b_dir = sys.argv[2]
#    idmapping_file = sys.argv[3]


class AnnotESM1b:
    def __init__(self, vt, esm1b_dir, idmapping_file):
        self.vt = vt
        self.esm1b_dir = esm1b_dir
        self.idmapping_file = idmapping_file
        self.df = pd.DataFrame()
        self.variant = list()
        #self.aachange = list()
        #self.idmapping = list()
        #self.uniprot = list()
        self.header = list()
        #self.my_dict = dict()
        self.ft = open("test_results.txt", 'w')

    def idmappingread(self):
        self.df = pd.read_table(self.idmapping_file, sep = "\t", header = None)
        self.df.columns = ["HUMAN", "refseq_NP", "UniRef50", "uniprot", "refseq_NM", "GENE"]
        
    ### refseq_aachange_dict ###
    # make my_dict hash which the key is refseq and the value is p.[REF][POS][ALT]
    def refseq_aachange_hash(self, aachange) :
        my_dict = {}
        for ele in aachange :
            key = ele.split(":")[1]
            value = ele.split(":")[4]
            my_dict[key] = value 
        return(my_dict)
    
    ### refseq_uniprot_hash ###
    # convert my_dict to nested dictionary, my_dict3
    # my_dict3 includes refseq as key, and uniprotid, pos, ref and alt as nested value. 
    def refseq_uniprot_hash(self, my_dict):
        my_dict3 = dict()
        for key, value in my_dict.items() :
            my_dict2 = dict()
            if (key in list(self.df.refseq_NM)) :
                uniprot = list(self.df[self.df.refseq_NM == key].uniprot)[0]
                #value = value.replace("p.", "")
                if (re.search(r'p.[A-Z]\d*[A-Z]', value) != None) : 
                    value = value.replace("p.", "")
                    aachange = re.findall(r'[A-Z]\d*[A-Z]', value)[0]
                    my_dict2["uniprot"] = uniprot
                    my_dict2["pos"] = re.sub('[A-Z]', "", aachange)
                    my_dict2["ref"] = list(aachange)[0]
                    my_dict2["alt"] = list(aachange)[-1]
                    my_dict3[key] = my_dict2
                else : pass
            else : pass
        return(my_dict3)

    ### add_esm1b ###
    def add_esm1b(self, my_dict3):
        temp = list() # temp list for esm1b column to be added in variant.table
        for k,v in my_dict3.items():
            esm1b = v["uniprot"] # uniprot id for searching db
            #print(k)
            #print(v["pos"])
            ref = v["ref"]; pos = int(v["pos"]); alt = v["alt"]
            print(esm1b)
            try:
                esm1b_path = os.path.join(self.esm1b_dir, esm1b+"_LLR.table.txt")
                esm1b_df = pd.read_table(esm1b_path, header = None)
            except FileNotFoundError: # if there are not matched LLR.txt file, pass the follow code
                                        # and then back to for loop 
                continue
        
            esm1b_path = os.path.join(self.esm1b_dir, esm1b+"_LLR.table.txt") # path for each uniprot esm1b
            #print("success")
            esm1b_df = pd.read_table(esm1b_path, header = None) # read esm1b db
            esm1b_df.columns = ["position", "ref", "alt", "LLR"] # define col names
            esm1b_df["position"] = esm1b_df["position"].astype(int) # check dtype
            esm1b_df["LLR"] = esm1b_df["LLR"].astype(float) # check dtype

    
            # find matched aachange
            results = esm1b_df[(esm1b_df.position == pos) & (esm1b_df.ref == ref) & (esm1b_df.alt == alt)]
            if not results.empty:
                pos = str(pos)
                LLR = str(list(results["LLR"])[0])
                temp.append(k+":"+esm1b+":"+ref+pos+alt+":"+LLR)
                #print(temp)
            else : pass  
    
        if temp: # matching well
            temp = ';'.join(temp)
            self.variant.append(temp) # add temp for esm1b column
            #print(self.variant)
        else: # thera are not any LLR evidences
            self.variant.append("None") # add "." for esm1b column

    def writeheader(self):
        self.header = '\t'.join(self.header)
        #print(self.header)
        self.ft.write(self.header+"\n")
        
    def writetable(self):
        self.variant = '\t'.join(self.variant)
        #print(self.variant)
        self.ft.write(self.variant+"\n")
        
    def main_func(self):
        # read idmapping file
        self.idmappingread()
        
        # open variant table
        f = open(self.vt, "rt")
        # header
        line = f.readline().strip()
        #print(line)
        self.header = line.split("\t")
        self.header.append("esm1b") # add new col on header
        #print(self.header)
        
        # get index of AA.change column
        idx = list(filter(lambda x:self.header[x] == "AAChange.refGeneWithVer", range(len(self.header))))
        #print(idx)
        idx = idx[-1] # idx is the index of "AAChange.refGeneWithVer"
        self.writeheader() # write header on output file
        
        # read a line following header
        count = 0
        line = f.readline().strip()
        while line != '': 
            count += 1
            #print(count)
            self.variant = line.split("\t")
            
            #print(self.variant)
            aachange = self.variant[idx].split(",")
            if (aachange[0] == "."): # some variants do not have any aachange information
                self.variant.append("None")
                self.writetable() # write the variant line to output file
            else : 
                my_dict = self.refseq_aachange_hash(aachange)
                my_dict3 = self.refseq_uniprot_hash(my_dict)
                if my_dict3:
                    self.add_esm1b(my_dict3)
                    for i in range(len(self.variant)): # convert dtype
                        self.variant[i] = str(self.variant[i])
                    self.writetable() # write the variant line to output file
                else : 
                    self.variant.append("None")
                    self.writetable()
            line = f.readline().strip()
        
        self.ft.close()
            

my_obj = AnnotESM1b("thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_0.001_metaSVM_plDiff8_variant.table.txt", "/Volumes/hjdrive/db/esm1b/esm1b_converted_hj", "/Volumes/hjdrive/db/esm1b/20230115_idmapping_table_esm1b_match.txt")

my_obj.main_func()



