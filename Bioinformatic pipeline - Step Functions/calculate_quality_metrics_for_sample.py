import pandas, sys

threshold_file = sys.argv[1]
regions_file = sys.argv[2]

thresholds = pandas.read_csv(threshold_file, sep="\t", compression="gzip", header=0)

regions = pandas.read_csv(regions_file, sep="\t", compression="gzip", header=None, names=["Chrom", "Start", "End", "Region", "AverageCov"], index_col=False).sort_values("AverageCov").reset_index()

regions["Length"] = regions["End"]-regions["Start"]

regions["CumulativeSum"] = regions["Length"].cumsum()/regions["Length"].sum()
median = regions["AverageCov"].iloc[regions[regions["CumulativeSum"]>0.5].index[0]-1]

open(sys.argv[3], "w").write(str(median) + "," + ",".join([str(x) for x in thresholds[["10X", "15X", "20X", "30X"]].sum()/sum(thresholds["end"]-thresholds["start"])]))