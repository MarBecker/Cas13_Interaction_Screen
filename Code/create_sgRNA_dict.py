import pandas as pd

#Importing the file with sgRNA sequences and corresponding gene names
guide_anno = pd.read_csv("../scree-Nraw_data/guide_anno.tsv", sep="\t")

#Removing the "_guideXX" from the "control_guideXX" names
guide_anno["spot1"] = guide_anno["spot1"].str.split("_").str[0]
guide_anno["spot2"] = guide_anno["spot2"].str.split("_").str[0]

#Creating a column with gene names, which are the concatenation of the two genes targeted
guide_anno["Gene"] = guide_anno["spot1"] + guide_anno["spot2"]

#Sort dataframe by gene to capture all combinations of two given genes in the subsequent step
guide_anno.sort_values("Gene", inplace = True)

#print(guide_anno.head(170))

#Adding consecuitve numbering to the same gene names, to provide unique names for the sgRNA arrays targeting the same two genes
genes = guide_anno["Gene"].tolist()
guides = []

for index, gene in enumerate(genes):
    if gene != genes[index-1]:
        sg = gene + "_1"
        guides.append(sg)
    elif gene == genes[index-1]:
        num = int(guides[index-1].split("_")[1])
        sg = gene + "_" + str(num+1)
        guides.append(sg)

#Adding the sgRNA names to the initial dataframe
guide_anno["sgRNA"] = guides

#print(guide_anno.head(10))

#Generating a dictionary from the raw sequences and the newly generated sgRNA names
#(Raw sequences are found in the read-count table and need to be replaced)
d = dict(zip(guide_anno["Unnamed: 0"], guide_anno["sgRNA"]))
print(list(d.items())[:10])

#Save the sgRNA dictionary to be able to track names to sequences manually
# with open("sgRNA_Dictionary.txt", "w") as file:
#     for key, value in d.items():
#         file.write(str(key) + " >>> " + str(value) + "\n")

#Load readcount table as dataframe

raw_counts = pd.read_csv("../screen_raw_data/all_raw.tsv", sep = "\t")


#Create column with sgRNA names based on the sequence column
raw_counts["sgRNA"] = raw_counts["Unnamed: 0"].map(d)

#Create column with Gene names based on sgRNA names
raw_counts["gene"] = raw_counts["sgRNA"].str.split("_").str[0]

#Remove column with raw sequences
raw_counts.drop("Unnamed: 0", axis = 1, inplace = True)

#Reorder columns to have sgRNA and gene as first 2 columns
cols = raw_counts.columns.tolist()
cols = cols[-2:] + cols[:-2]
raw_counts = raw_counts[cols]
raw_counts.set_index("sgRNA", inplace = True)

#print(raw_counts.columns)
#print(raw_counts.head())

#Export read count table
raw_counts.to_csv("../screen_raw_data/raw_counts_edit.txt", sep = "\t")

# raw_HEK = raw_counts[["sgRNA", "gene", "HEKTP1.REP1", "HEKTP1.REP2", "HEKTP3.REP1", "HEKTP3.REP2"]]
# #print(raw_HEK.head())
# raw_HEK.set_index("sgRNA", inplace = True)
# raw_HEK.to_csv("HEK_read_counts.txt", sep = "\t")