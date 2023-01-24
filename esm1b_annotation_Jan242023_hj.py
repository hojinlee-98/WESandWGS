import sys
import re
import os
import pandas as pd
import gzip

### usage ###
# python esm1b_annotation_Jan242023_hj.py [variant table] [converted esm1b directory] [idmapping file]


if len(sys.argv) != 4:
    print("please, check the arguments")
else:
    vcf_gz = sys.argv[1]
    esm1b_dir = sys.argv[2]
    idmapping_file = sys.argv[3]

class AnnotESM1b:
    def __init__(self, vcf, esm1b_dir, idmapping_file):
        self.f = None
        self.vcf = vcf
        self.line = None
        self.esm1b_dir = esm1b_dir
        self.idmapping_file = idmapping_file
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
        for ele in aachange:
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
            if (bool(uniprot) & (re.search(r'p.[A-Z]\d*[A-Z]$', value) != None)) :
                # if there is at least one matched ID and AAchange express missense
                # (fs/non-fs INDEL do not run the follow contexts, then their my_dict3 express None)
                self.df_match = self.df[(self.df.uniprot.isin(uniprot)) & (self.df.refseq_NM == key_edit)] # matched idmapping table
                self.refseq_uniprot_hashing(key, value)

    ### refseq_uniprot_hashing ###
    # convert my_dict to nested dictionary, my_dict3
    # called by idmapping_explorer()
    # my_dict3 includes refseq as key, and uniprotid, pos, ref and alt as nested value. 
    def refseq_uniprot_hashing(self, key, value):
        value = value.replace("p.", "") # p.[A-Z][0-9]*[A-Z] to [A-Z][0-9]*[A-Z]
        for i in range(self.df_match.shape[0]): # iter nrow (.shape[0])
            my_dict2 = dict() # init my_dict2
            aachange = re.findall(r'[A-Z]\d*[A-Z]$', value)[0] # find the pattern in value
            my_dict2["uniprot"] = self.df_match["uniprot"].iloc[i]
            my_dict2["refseq_aachange"] = key
            my_dict2["refseq_idmapping"] = self.df_match["refseq_NM.version"].iloc[i]
            my_dict2["pos"] = re.sub('[A-Z]', "", aachange) # remove UPPER from aachange, then pos only remained
            my_dict2["ref"] = list(aachange)[0]
            my_dict2["alt"] = list(aachange)[-1]
            my_dict3_key = key+"_"+self.df_match["uniprot"].iloc[i] # ex) my_dict3_key means NM_123.1_ABSX-2
            self.my_dict3[my_dict3_key] = my_dict2

    ### add_esm1b ###
    # which variants are handled by this function?
    #      only missense variants are used, some of which may not match with ESM1b converted DB,
    #      in which case "esm1b=." and "min_esm1b=." are added as elements to the self.line
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

            # find matched pos,ref,alt
            results = esm1b_df[(esm1b_df.position == pos) & (esm1b_df.ref == ref) & (esm1b_df.alt == alt)]
            if not results.empty:
                pos = str(pos)
                LLR = str(list(results["LLR"])[0])
                temp.append(refseq_aachange+":"+refseq_idmapping+":"+esm1b+":"+ref+pos+alt+":"+LLR)
            else : pass

        if temp: # matching well
            temp = ','.join(temp)
            self.line.append("esm1b="+temp) # add temp for esm1b column
            self.most_deleterious_LLR()
        elif not temp: # there are not any LLR evidences (1. not matched LLR.table, 2. not matched pos,ref,alt)
            self.line.append("esm1b=.") # add "." for esm1b column
            self.most_deleterious_LLR()

    ### most_deleterious_LLR ###
    # called by add_esm1b()
    # present minimum value of the LLR scores
    def most_deleterious_LLR(self):
        esm1b = self.line[-1]
        if esm1b != "esm1b=." :
            LLR_list = list(map(lambda x:float(x.split(":")[-1]), esm1b.split(","))) # split esm1b field and extract LLRs as list
            LLR_min = min(LLR_list)
            self.line.append("min_esm1b="+str(LLR_min))
        elif esm1b == "esm1b=." :
            self.line.append("min_esm1b=.")

    def writetable(self):
        self.line = list(map(str, self.line)) # change dtype 
        self.line = ";".join(self.line)
        self.ft.write(self.line+"\n")

    def header_editer(self):
        while list(self.line)[0] == "#" :
            if re.search("ALLELE_END", self.line) == None :
                self.ft.write(self.line+"\n")
            elif re.search("ALLELE_END", self.line) != None :
                self.ft.write("##INFO=<ID=esm1b,Number=.,Type=String,Description=\"esm1b annotation provided by hj\'s script\">"+"\n")
                self.ft.write("##INFO=<ID=min_esm1b,Number=.,Type=String,Description=\"min_esm1b annotation provided by hj\'s script\">"+"\n")
                self.ft.write(self.line+"\n")
            self.line = self.f.readline().strip()

    def main_func(self):
        # read idmapping file
        self.idmappingread()

        # open variant table
        self.f = gzip.open(self.vcf, "rt")
        self.line = self.f.readline().strip()
        self.header_editer()

        while self.line != "":
            self.line = self.line.split(";")
            aachange_idx = [i for i in range(len(self.line)) if "AAChange.refGeneWithVer" in self.line[i]][0] # assign AAChange indexes
            aachange = self.line[aachange_idx].split("=")[-1]
            missing_aachange = [".", "UNKNOWN"]
            if aachange in missing_aachange: # if the variant only has dot in the AAchange columns
                # CASE1) there are not any AAChange annotation
                self.line.append("esm1b=.")
                self.line.append("min_esm1b=.")
                self.writetable()
            elif aachange != "." : # if the variant has the AAchange value 
                aachange = aachange.split(",") # split aachange using ",", now aachange is list object
                my_dict = self.refseq_aachange_hashing(aachange)
                self.idmapping_explorer(my_dict)
                if self.my_dict3: # if self.my_dict3 is not None
                    self.add_esm1b()
                    # the following two CASES are handled in self.add_esm1b()
                    # CASE2-1) no problem in annotation ESM1b 
                    # CASE2-2) besides there is AAChange value,
                    #          at least one of RefSeqNM is matched with idmapping table
                    # CASE2-3) there is AAChange value, but there are not any LLR evidences
                    self.writetable()
                elif not self.my_dict3: # if self.my_dict3 is None
                    # CASE3-1) if the variant confers fs/non-fs INDEL consequence,
                    #          ESM1b can not be used to explain the deleteriousness
                    # CASE3-2) if the variant do not matched with idmapping table. 
                    self.line.append("esm1b=.")
                    self.line.append("min_esm1b=.")
                    self.writetable()
            self.line = self.f.readline().strip()
        self.ft.close()

my_obj = AnnotESM1b(vcf_gz, esm1b_dir, idmapping_file)
my_obj.main_func()
