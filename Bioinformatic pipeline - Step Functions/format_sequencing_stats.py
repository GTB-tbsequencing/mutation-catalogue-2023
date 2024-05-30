import sys, pandas

file_location = sys.argv[1]

sample_name = file_location.rsplit("/", 1)[1]

header_regions_qc = [
    "Median", 
    "Coverage10x",
    "Coverage15x",
    "Coverage20x",
    "Coverage30x"
]

stats = pandas.read_csv(file_location+"_qcstats.txt", sep=",", header=None, names=header_regions_qc)

samtools_stats = []

with open(file_location+"-samtools-stats.txt", "r") as samtools_file:
    for line in [x for x in samtools_file.readlines() if x.startswith("SN")]:
        try:
            samtools_stats.append(int(line.split(":")[1].split("#")[0].strip()))
        except ValueError:
            samtools_stats.append(float(line.split(":")[1].split("#")[0].strip()))

rec = pandas.DataFrame.from_records([[sample_name] + list(list(stats[header_regions_qc].itertuples(index=False, name=None))[0]) + samtools_stats])

rec[[2, 3, 4, 5]] = rec[[2, 3, 4, 5]].round(4)

rec.to_csv("global-stats/"+sample_name+"-stats.csv.gz", sep="\t", header=False, index=False, compression="gzip")

header_regions_bed = [
    "Chromosome", 
    "Start",
    "End",
    "Locus",
    "MeanDepth"
]

mean_cov = pandas.read_csv(file_location + ".regions.bed.gz", sep="\t", header=None, names=header_regions_bed, compression="gzip", index_col=3)

depth_locus = pandas.read_csv(file_location + ".thresholds.bed.gz", sep="\t", header=0, compression="gzip").rename(columns={"region": "Locus"}).set_index("Locus")

for i in ["10X", "15X", "20X", "30X"]:
    depth_locus[i] = depth_locus[i]/(depth_locus["end"]-depth_locus["start"])

mean_cov = mean_cov.join(depth_locus).reset_index()

mean_cov["Sample"] = sample_name

mean_cov[["Sample", "Locus", "MeanDepth", "10X", "15X", "20X", "30X"]].to_csv("locus-stats/"+sample_name+"-locusseqstats.csv.gz", sep="\t", header=False, index=False, compression="gzip")