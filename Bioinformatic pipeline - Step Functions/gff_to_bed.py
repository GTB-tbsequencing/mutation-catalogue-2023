import argparse
from BCBio import GFF


#This script simply takes a GFF file and converts it
#into a bed file with 4th column as locus tag

def main(arguments):

    coordinates = []

    with open(arguments.input_gff_file) as gff_handle:
        for seq in GFF.parse(gff_handle):
            for features in seq.features:
                if features.type=="gene":
                    start = str(features.location.start)
                    end = str(features.location.end)
                    name = features.qualifiers["locus_tag"][0]               
                    coordinates.append((seq.id, start, end, name))

    with open(arguments.output_bed_file, "w") as bed_handle:
        bed_handle.write("\n".join(["\t".join(line) for line in coordinates]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_gff_file")
    parser.add_argument("--output_bed_file")
    args = parser.parse_args()

    main(args)