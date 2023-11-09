import argparse, pandas, numpy

def main(arguments):

    annotated_file = arguments.file_name

    df = pandas.read_csv(annotated_file, sep="\t", header=None, names=["VariantId", "Annotation"], index_col=[0], dtype={"Id":int, "Annotation":str}, na_filter=False, low_memory=False)

    if df.empty:
        open("annotation_normalized.tsv", 'a').close()
    else:
        annotation = pandas.DataFrame(df.Annotation.str.split(";", expand=True)[0].str.split(",").tolist(), index=df.index).stack().str.split("|", expand=True)
        if annotation.shape[1]==19:
            annotation.columns = [
                "Allele",
                "Annotation",
                "Effect",
                "GeneName",
                "GeneId",
                "FeatureType",
                "FeatureId",
                "Biotype",
                "Rank",
                "HGVS.c",
                "HGVS.p",
                "cDNA_len",
                "CDS_len",
                "Protein_len",
                "Distance",
                "Warning1",
                "Warning2",
                "Warning3",
                "Warning4"
                ]
        #Sometimes the resulting data will lack the columns at the end
        elif annotation.shape[1]==16:
            annotation.columns = [
                "Allele",
                "Annotation",
                "Effect",
                "GeneName",
                "GeneId",
                "FeatureType",
                "FeatureId",
                "Biotype",
                "Rank",
                "HGVS.c",
                "HGVS.p",
                "cDNA_len",
                "CDS_len",
                "Protein_len",
                "Distance",
                "Warning1"
                ]
        else:
            raise ValueError

    annotation["Allele"] = annotation["Allele"].str.replace("ANN=", "")
    annotation["GeneId"] = annotation["GeneId"].str.replace("gene-", "")
    annotation["CDS_len"] = annotation["CDS_len"].str.split(pat="/",expand=True)[0]
    annotation["Protein_len"] = annotation["Protein_len"].str.split(pat="/",expand=True)[0]

    # SnpEff mishandles MNV that could be related to "Stop_retained_variant"
    # We have to fix that in various steps

    # Fix stop retained variant that are in fact early stop codon
    regexp_missense_ter = r'^p\.[A-Z][a-z]{2}([0-9]+)Ter$'
    annotation.loc[(annotation["HGVS.p"].str.match(regexp_missense_ter)) & (annotation["Annotation"]=="stop_retained_variant"), "Annotation"] = "stop_gained"
    annotation.loc[annotation["HGVS.p"].str.match(regexp_missense_ter) & (annotation["Annotation"]=="stop_gained"), "HGVS.p"] = annotation.loc[annotation["HGVS.p"].str.match(regexp_missense_ter) & (annotation["Annotation"]=="stop_gained"), "HGVS.p"].str.replace(r"Ter$", "*", regex=True) 

    # Fix stop retained variant that are missense
    regexp_missense = r'^p\.[A-Z][a-z]{2}([0-9]+)[A-Z][a-z]{2}$'
    annotation.loc[annotation["HGVS.p"].str.match(regexp_missense) & (annotation["Annotation"]=="stop_retained_variant"), "Annotation"] = "missense_variant"




    #Normalization of the table happening here:
    #Flattening of all annotations, nucleotidic (eg c.609G>A) and proteic (eg p.Ser315Thr)
    annotation = pandas.DataFrame(annotation.reset_index().set_index(["VariantId", "Annotation", "GeneId", "Distance", "CDS_len", "Protein_len"]).replace(r'^\s*$', numpy.nan, regex=True).rename(columns={"HGVS.c":"nucleotidic", "HGVS.p":"proteic"})[["nucleotidic", "proteic"]].stack().dropna()).reset_index().rename(columns={0:"Value"})

    annotation.loc[(annotation["level_6"]=="nucleotidic"), "Dist"] = annotation.loc[(annotation["level_6"]=="nucleotidic"), "CDS_len"]

    annotation.loc[(annotation["Annotation"]=="upstream_gene_variant") | (annotation["Annotation"]=="downstream_gene_variant"), "Dist"] = annotation["Distance"]

    annotation.loc[(annotation["level_6"]=="proteic"), "Dist"] = annotation["Protein_len"]

    # annotation = annotation.replace(r'^\s*$', numpy.nan, regex=True)

    # annotation["Dist"] = annotation["Dist"].astype(int)

    #Large deletions are annotated as feature ablations
    #All affected locus are in the GeneId column, separated with "&"
    feature_ablations = annotation[(annotation["Annotation"]=="feature_ablation")|(annotation["Annotation"]=="transcript_ablation")].reset_index().set_index(['index', 'VariantId', 'Annotation', 'Distance', 'CDS_len', 'Protein_len', 'level_6', "Value", 'Dist'])

    try:
        feature_ablations = pandas.DataFrame(feature_ablations.GeneId.str.split("&").tolist(), index=feature_ablations.index).stack().reset_index().drop(["level_9", "index"], axis=1).rename(columns={0:"GeneId"})
        feature_ablations["GeneId"]=feature_ablations["GeneId"].str.replace("gene-", "")
        annotation = pandas.concat([annotation.drop(annotation[(annotation["Annotation"]=="intergenic_region") | (annotation["Annotation"]=="feature_ablation") | (annotation["Annotation"]=="gene_fusion") | (annotation["Annotation"]=="bidirectional_gene_fusion") | (annotation["Annotation"]=="chromosome_number_variation") | (annotation["Annotation"]=="transcript_ablation")].index), feature_ablations])
    except (IndexError, AttributeError) as error:
        annotation = annotation.drop(annotation[(annotation["Annotation"]=="intergenic_region") | (annotation["Annotation"]=="feature_ablation") | (annotation["Annotation"]=="gene_fusion") | (annotation["Annotation"]=="bidirectional_gene_fusion") | (annotation["Annotation"]=="chromosome_number_variation") | (annotation["Annotation"]=="transcript_ablation")].index)
        
    #Corrections of a few SnpEff systematic errors
    
    #When long variants only lead to synonymous changes, sometimes the HGVS protein notation is wrong (eg p.465)
    #We replace it with p.= which is valid HGVS
    annotation["Value"] = annotation["Value"].str.replace(r"^p\.[0-9]+$", "p.=", regex=True)
    
    #Sometimes consecutive missense variants are lumped into one HGVS annotation, which is incorrect (eg p.AspGly95LeuSer)
    #We split those to replace it with two distinct entries : p.Asp95Leu and p.Gly96Ser
    #We do not keep intermediate synonymous annotation relative to the protein, i.e p.AspGlySer95LeuGlyGln becomes p.Asp95Leu and p.Ser97Gln
    #We also deal with possible stop codon mutations
    regexp_multiple_missense = r'^p\.([A-Z][a-z]{2}){2,}[0-9]+([A-Z][a-z]{2}){2,}$'
    
    regexp_stop_gained = r'^p\.[A-Z][a-z]{2}([0-9]+)\*$'
    regexp_ter_missense = r'^p\.Ter([0-9]+)[A-Z][a-z]{2}$'

    regexp_multiple_missense_start_lost = r'^p\.(Met|Val|Leu)(?:[A-Z][a-z]{2})+1\?$'
    annotation.loc[annotation["Value"].str.match(regexp_multiple_missense_start_lost) & (annotation["Annotation"]=="start_lost"), "Value"] = "p." + annotation["Value"].str.extract(regexp_multiple_missense_start_lost, expand=False) + "1?"


    multiple_missense_mutants = annotation.loc[annotation["Value"].str.match(regexp_multiple_missense)].copy()
    if not multiple_missense_mutants.empty:
        annotation = annotation[~annotation["Value"].str.match(regexp_multiple_missense)]
        multiple_missense_mutants["Order"] = (multiple_missense_mutants["Value"].str.split(r'[0-9]+', expand=True)[1].str.len()/3).astype(int).apply(lambda x:range(x))
        multiple_missense_mutants[["Ref", "Alt"]] = multiple_missense_mutants["Value"].str.split(r'[0-9]+', expand=True).apply(lambda x: x.str.replace(r"^p\.", "", regex=True).str.findall(r"[A-z][a-z]{2}"))
        multiple_missense_mutants = multiple_missense_mutants.explode(["Ref", "Alt", "Order"])
        multiple_missense_mutants = multiple_missense_mutants[multiple_missense_mutants["Ref"]!=multiple_missense_mutants["Alt"]]
        #Sometimes Dist and the number in the HGVS nomenclature are not equal, when the first change is missense...
        #So let's use the HGVS nomenclature value instead 
        multiple_missense_mutants["Dist"] = multiple_missense_mutants["Value"].str.extract("([0-9]+)").astype(int)
        multiple_missense_mutants["Dist"] = multiple_missense_mutants["Dist"]  +  multiple_missense_mutants["Order"].astype(int)
        multiple_missense_mutants["Value"] = "p." + multiple_missense_mutants["Ref"] + multiple_missense_mutants["Dist"].astype(str) +  multiple_missense_mutants["Alt"].str.replace("Ter", "*", regex=False)

        multiple_missense_mutants.loc[multiple_missense_mutants["Value"].str.match(regexp_stop_gained), "Annotation"] = "stop_gained"
        multiple_missense_mutants.loc[multiple_missense_mutants["Value"].str.match(regexp_missense), "Annotation"] = "missense_variant"
        multiple_missense_mutants.loc[multiple_missense_mutants["Value"].str.match(regexp_ter_missense), "Annotation"] = "stop_lost"

        annotation = pandas.concat([annotation[["VariantId", "Annotation", "GeneId", "Dist", "level_6", "Value"]], multiple_missense_mutants[["VariantId", "Annotation", "GeneId", "Dist", "level_6", "Value"]]])

    #For contiguous synonymous + missense mutation, the position of the first codon is reported, instead of the position of the missense codon.
    #We replace those instances with the correct values
    
    annotation.loc[annotation["Value"].str.match(regexp_missense), "Dist"] = annotation["Value"].str.extract(regexp_missense, expand=False).dropna()

    annotation[["VariantId", "Annotation", "GeneId", "Dist", "level_6", "Value"]].to_csv("annotation_normalized.tsv", header=False, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file_name", help="File name of SnpEff output")
    args = parser.parse_args()

    main(args)