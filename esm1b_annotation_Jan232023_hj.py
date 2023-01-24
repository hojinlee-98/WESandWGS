import sys
import re
import os
import pandas as pd

if len(sys.argv) != 4:
    print("please, check the input arguments")
else:
    vt = sys.argv[1]
    esm1b_dir = sys.argv[2]
    idmapping_file = sys.argv[3]

class AnnotESM1b:
    def __init__(self, vt, esm1b_dir, idmapping_file):
        self.vt = vt
        self.esm1b_dir = esm1b_dir
        self.idmapping_file = idmapping_file
        self.variant = list()
        self.header = list()
        self.df_match = pd.DataFrame()
        self.my_dict3 = dict()
        self.ft = open("test_results.txt", 'w')

    def idmappingread(self):
        self.df = pd.read_table(self.idmapping_file, sep = "\t", header = None)
        self.df.columns = ["uniprot", "NCBI","refseq_NM.version", "refseq_NM"]
        
    ### refseq_aachange_hashing ###
    # make my_dict hash which the key is refseq and the value is p.[REF][POS][ALT]
    def refseq_aachange_hashing(self, aachange):
        my_dict = {}
        for ele in aachange :
            key = ele.split(":")[1] # Refseq NM
            value = ele.split(":")[4] # p.AAchange
            my_dict[key] = value 
        return(my_dict)
    
    ### idmapping_explorer ###
    # find RefSeq in df (idmapping) and return the uniprot IDs.
    # if there is at least one matched ID and AAchange express missense,
    # assign the new dataframe matched with uniprot IDs and Refseq ID to self.df_match
    def idmapping_explorer(self, my_dict):
        self.my_dict3 = dict()
        for key, value in my_dict.items(): # the key means RefSeq and the value means p.AAchange
            key_edit = re.sub("\.[0-9]$", "", key)
            uniprot = list(self.df[self.df.refseq_NM == key_edit].uniprot)
            if ((uniprot != None) & (re.search(r'p.[A-Z]\d*[A-Z]$', value) != None)) :
                # if there is at least one matched ID and AAchange express missense
                self.df_match = self.df[(self.df.uniprot.isin(uniprot)) & (self.df.refseq_NM == key_edit)]
                self.refseq_uniprot_hashing(key, value)
        
    ### refseq_uniprot_hashing ###
    # convert my_dict to nested dictionary, my_dict3
    # called by idmapping_explorer()
    # my_dict3 includes refseq as key, and uniprotid, pos, ref and alt as nested value. 
    def refseq_uniprot_hashing(self, key, value):
        value = value.replace("p.", "")
        for i in range(self.df_match.shape[0]):
            my_dict2 = dict()
            aachange = re.findall(r'[A-Z]\d*[A-Z]$', value)[0]
            my_dict2["uniprot"] = self.df_match["uniprot"].iloc[i]
            my_dict2["refseq_aachange"] = key
            my_dict2["refseq_idmapping"] = self.df_match["refseq_NM.version"].iloc[i]
            my_dict2["pos"] = re.sub('[A-Z]', "", aachange)
            my_dict2["ref"] = list(aachange)[0]
            my_dict2["alt"] = list(aachange)[-1]
            my_dict3_key = key+"_"+self.df_match["uniprot"].iloc[i]
            self.my_dict3[my_dict3_key] = my_dict2

    ### add_esm1b ###
    def add_esm1b(self):
        temp = list() # temp list for esm1b column to be added in variant.table
        for k,v in self.my_dict3.items():
            esm1b = v["uniprot"] # assign uniprot id for searching db to esm1b
            refseq_aachange = v["refseq_aachange"]
            refseq_idmapping = v["refseq_idmapping"]
            ref = v["ref"]; pos = int(v["pos"]); alt = v["alt"]
            try:
                esm1b_path = os.path.join(self.esm1b_dir, esm1b+"_LLR.table.txt")
                esm1b_df = pd.read_table(esm1b_path, header = None)
            except FileNotFoundError: # if there are not matched LLR.txt file, pass the follow code
                                        # and then back to for loop 
                continue

            esm1b_df.columns = ["position", "ref", "alt", "LLR"] # define col names
            esm1b_df["position"] = esm1b_df["position"].astype(int) # check dtype
            esm1b_df["LLR"] = esm1b_df["LLR"].astype(float) # check dtype
            
            # find matched aachange
            results = esm1b_df[(esm1b_df.position == pos) & (esm1b_df.ref == ref) & (esm1b_df.alt == alt)]
            if not results.empty:
                pos = str(pos)
                LLR = str(list(results["LLR"])[0])
                temp.append(refseq_aachange+":"+refseq_idmapping+":"+esm1b+":"+ref+pos+alt+":"+LLR)
            else : pass  
    
        if temp: # matching well
            temp = ';'.join(temp)
            self.variant.append(temp) # add temp for esm1b column
            self.most_deleterious_LLR()
        else: # there are not any LLR evidences
            self.variant.append(".") # add "." for esm1b column
            self.most_deleterious_LLR()
            
    ### most_deleterious_LLR ###
    # called by add_esm1b()
    # present minimum value of the LLR scores
    def most_deleterious_LLR(self):
        esm1b = self.variant[-1]
        if esm1b != "." :
            esm1b_split = esm1b.split(";")
            LLR_list = list()
            for ele in esm1b_split:
                LLR = ele.split(":")[-1]
                #LLR = float(re.findall(r'-[0-9]*.[0-9]*$', ele)[0])
                LLR_list.append(LLR)
            LLR_min = min(LLR_list)
            self.variant.append(LLR_min)
        else : self.variant.append(".")
            
    def writeheader(self):
        self.header = '\t'.join(self.header)
        self.ft.write(self.header+"\n")
        
    def writetable(self):
        self.variant = '\t'.join(self.variant)
        self.ft.write(self.variant+"\n")
        
    def main_func(self):
        # read idmapping file
        self.idmappingread()
        
        # open variant table
        f = open(self.vt, "rt")
        # header
        line = f.readline().strip()
        self.header = line.split("\t")
        self.header.append("esm1b") # add new col on header
        self.header.append("min_esm1b") # add min LLR on header

        # get index of AA.change column
        idx = list(filter(lambda x:self.header[x] == "AAChange.refGeneWithVer", range(len(self.header))))
        idx = idx[-1] # idx is the index of "AAChange.refGeneWithVer"
        self.writeheader() # write header on output file
        
        # read a line following header
        count = 0
        line = f.readline().strip()
        while line != '': 
            self.variant = line.split("\t")
            aachange = self.variant[idx].split(",")
            if (aachange[0] == "."): # some variants do not have any aachange information
                self.variant.append(".")
                self.writetable() # write the variant line to output file
            else : 
                my_dict = self.refseq_aachange_hashing(aachange)
                self.idmapping_explorer(my_dict)
                if self.my_dict3:
                    self.add_esm1b()
                    for i in range(len(self.variant)): # convert dtype
                        self.variant[i] = str(self.variant[i])
                    self.writetable() # write the variant line to output file
                else : 
                    self.variant.append(".")
                    self.writetable()
            line = f.readline().strip()
        
        self.ft.close()
            

my_obj = AnnotESM1b(vt, esm1b_dir, idmapping_file)
my_obj.main_func()
