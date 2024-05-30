import sys, pandas

file_location = sys.argv[1]

sample_name = file_location.rsplit("/", 1)[1]

header_taxonomy = [
    "Percent",
    "RawNumber1", 
    "RawNumber2",
    "Rank", 
    "NCBI_Taxon_Id",
    "Name"
]

taxonomy = pandas.read_csv(file_location+"-kraken.txt", sep="\t", header=None, names=header_taxonomy, index_col=None)

taxonomy["SampleId"] = sample_name

taxonomy[["SampleId", "NCBI_Taxon_Id", "Percent"]].to_csv("taxonomy-assignment/"+sample_name+"-taxonomy.csv.gz", header=False, index=False, sep="\t", compression="gzip")