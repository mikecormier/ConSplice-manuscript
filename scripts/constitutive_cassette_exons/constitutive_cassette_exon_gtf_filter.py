from __future__ import print_function

import argparse
import gzip
import io
import os
import sys
from collections import defaultdict

from consplice.constraint.utils import get_alternative_gene_symbols
from interlap import InterLap

# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
# ------------------------------------------------------------------------------------------------------------------------------------------------------


def arguments():

    p = argparse.ArgumentParser(
        description="Parse a GTF file and extract all constitutive and cassette exons using HEXEvent info files"
    )

    req = p.add_argument_group("Required Arguments")

    req.add_argument(
        "--gtf",
        metavar="GRCh38 gtf file",
        required=True,
        help="(Required) The file path to a GRCh38 gtf file with exon features to extract. NOTE: Both this file and the HEXEvent files need to be in the b38 of the human genome",
    )

    req.add_argument(
        "--constitutive",
        metavar="HEXEvent Constitutive Exons",
        required=True,
        help="(Required) The file path to the GRCh38 constitutive exons downloaded from HEXEvent. NOTE: must be in b38 of the human genome",
    )

    req.add_argument(
        "--cassette",
        metavar="HEXEvent Cassette Exons",
        required=True,
        help="(Required) The file path to the GRCh38 cassette exons downloaded from HEXEvent. NOTE: must be in b38 of the human genome",
    )

    req.add_argument(
        "--alt-gene-symbol",
        metavar="Alternative Gene Symbol File",
        required=True,
        help="(Required) A file that contains mappings between a canonical gene symbol to alternative gene symbols. NOTE: The file needs to have a header!. This script is set up to use the HGNC protein-coding gene mapping file: ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt, however, you are welcome to try others. There is no guarantee it will work with other mapping files.",
    )

    req.add_argument(
        "--out-file",
        metavar="Output File",
        required="True",
        help="(Required) The file path and/or name of the output gtf file to create. This file will contain only constitutive and cassette exons features from the gtf file",
    )

    p.add_argument(
        "--coordinate-based",
        metavar="HEXEvent coordinate based file",
        choices=["0-based", "1-based"],
        default="0-based",
        help="(Optional) Which coordinate base the HEXEvent files are using. Direct download files are usually 0-based files. Choices = ['0-based',1-based']. Default = '0-based'",
    )

    return p.parse_args()


# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Functions
# ------------------------------------------------------------------------------------------------------------------------------------------------------


def create_hexevent_interlap(file_path, zero_based):
    """
    create_hexevent_interlap
    ========================
    This method is used to create interlap objects of HEXEvent features. (Constitutive or cassette exons). It will iterate through
     each feature in the file and add it to the chromosome specific interlap object in the dictionary. It uses the names of the columns
     from the HEXEvent files. This function will not work if the column names have been changed. Removes any chr prefix if it exists.
     The start position is set to 1-based if it is designated as a 0-based start position.

    Parameters:
    -----------
    1) file_path: (str) Path to the HEXEvent txt file (of constitutive or cassette exon features).
    2) zero_based: (bool) True if zero-based start position else False


    Returns:
    ++++++++
    1) (dict) A dict of interlap objects. Each interlap object contains of tuple of:
                (1-based start pos, end-pos, {strand, constitLevel, inclLevel, genename})
    """

    interlap_dict = defaultdict(InterLap)

    try:
        fh = (
            gzip.open(file_path, "rt", encoding="utf-8")
            if file_path.endswith(".gz")
            else io.open(file_path, "rt", encoding="utf-8")
        )
    except IOError as e:
        print("\n!!ERROR!! unable to open file: '{}'".format(file_path))
        print(str(e))
        sys.exit(1)

    header = fh.readline().strip().split("\t")

    for line in fh:

        line_dict = dict(zip(header, line.strip().split("\t")))

        ## Skip any lines that are missing data
        if "start" not in line_dict:
            continue

        start = int(line_dict["start"]) + 1 if zero_based else int(line_dict["start"])
        end = int(line_dict["end"])

        ## create interlap object
        interlap_dict[line_dict["chromo"].replace("chr", "")].add(
            (
                start,
                end,
                {
                    "strand": line_dict["strand"],
                    "constitLevel": line_dict["constitLevel"],
                    "inclLevel": line_dict["inclLevel"],
                    "genename": line_dict["genename"],
                },
            )
        )

    fh.close()

    return interlap_dict


def overlapping_exon(
    feature_chrom,
    feature_start,
    feature_end,
    feature_strand,
    feature_gene,
    alt_gene_dict,
    interlap_object,
):
    """
    overlapping_exon
    ================
    This method is used to identify exons from a gtf file that overlap a feature from an interlap object. This method assumes
     that the start positions of the interlap object are 1-based, and that the chromosome name does not have a 'chr' prefix.
     It also will only identify if an overlap if the exon has the same start, end, strand, and gene name.

    This method expects that the interlap object is based on the HEXEvent data file and that it contains a dictionary with the 'strand',
    'constitLevel', 'inclLevel', 'strand' feature values.

    Parameters:
    -----------
    1) feature_chrom:   (str)  The gtf chromosome name
    2) feature_start:   (int)  The start position of the gtf feature
    3) feature_end:     (int)  The end position of the gtf feature
    4) feature_strand:  (str)  The strand of the gtf feature
    5) feature_gene:    (str)  The name of the gene for the gtf feature
    6) alt_gene_dict:   (dict) A dictionary of alt gene symbols
    7) interlap_object: (dict) A dictionary of interlap object features to identify overlaps from

    Returns:
    ++++++++
    1) (bool) True or False whether a matching overlap feature was found
    2) (str)  The constitLevel value from the overlapping feature or an empty string
    3) (str)  The inclLevel value from the overlapping feature or an empty string
    """

    for exon in interlap_object[feature_chrom.replace("chr", "")].find(
        (int(feature_start), int(feature_end))
    ):

        same_start = int(exon[0]) == int(feature_start)
        same_end = int(exon[1]) == int(feature_end)
        same_strand = exon[2]["strand"] == feature_strand
        same_gene = any(
            True if x == feature_gene or x in alt_gene_dict[feature_gene] else False
            for x in exon[2]["genename"].strip().split(",")
        )

        return (
            True if same_start and same_end and same_strand and same_gene else False,
            exon[2]["constitLevel"],
            exon[2]["inclLevel"],
        )

    return (False, "", "")


# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Main
# ------------------------------------------------------------------------------------------------------------------------------------------------------


def main():

    args = arguments()

    print(
        (
            "\nInput Arguments:"
            "\n================"
            "\n - gtf:               {}"
            "\n - constitutive:      {}"
            "\n - cassette:          {}"
            "\n - alt-gene-symbol    {}"
            "\n - out-file:          {}"
            "\n - coordinate-based:  {}"
        ).format(
            args.gtf,
            args.constitutive,
            args.cassette,
            args.alt_gene_symbol,
            args.out_file,
            args.coordinate_based,
        )
    )

    print("\nProcessing")
    print("==========")

    print("\n\tCreating Constitutive and Cassette exon interlap objects")
    constitutive_interlap = create_hexevent_interlap(
        args.constitutive, True if args.coordinate_based == "0-based" else False
    )
    cassette_interlap = create_hexevent_interlap(
        args.cassette, True if args.coordinate_based == "0-based" else False
    )

    print("\n\tGetting alternative gene symbols")
    alt_symbol_dict = get_alternative_gene_symbols(args.alt_gene_symbol)

    kept_exons = []

    print("\n\tOpening gtf file")
    try:
        gtf_fh = (
            gzip.open(args.gtf, "rt", encoding="utf-8")
            if args.gtf.endswith(".gz")
            else io.open(args.gtf, "rt", encoding="utf-8")
        )
    except IOError as e:
        print("\n!!ERROR!! unable to open the gtf file: '{}'".format(args.gtf))
        print(str(e))
        sys.exit(1)

    header = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    print("\n\tIdentifying constitutive and cassette exons in gtf file")
    for line in gtf_fh:

        if line[0] == "#":
            continue

        line_dict = dict(zip(header, line.strip().split("\t")))
        line_dict.update(
            {
                x.strip()
                .replace('"', "")
                .split(" ")[0]: x.strip()
                .replace('"', "")
                .split(" ")[1]
                for x in line_dict["attribute"].strip().split(";")[:-1]
            }
        )

        ## only look at exon features
        if line_dict["feature"] != "exon":
            continue

        gene_name = (
            line_dict["gene_name"]
            if "gene_name" in line_dict
            else line_dict["gene_symbol"]
        )

        exon_added = False

        ## Check for constitutive exon
        overlapped, constitLevel, inclLevel = overlapping_exon(
            line_dict["chrom"],
            line_dict["start"],
            line_dict["end"],
            line_dict["strand"],
            gene_name,
            alt_symbol_dict,
            constitutive_interlap,
        )
        if overlapped and not exon_added:
            new_line = (
                line.strip()
                + " exon_type {}; constitLevel {}; inclLevel {};".format(
                    "constitutive", constitLevel, inclLevel
                )
            )
            kept_exons.append(new_line)
            exon_added = True

        ## Check for cassette exon
        overlapped, constitLevel, inclLevel = overlapping_exon(
            line_dict["chrom"],
            line_dict["start"],
            line_dict["end"],
            line_dict["strand"],
            gene_name,
            alt_symbol_dict,
            cassette_interlap,
        )

        if overlapped and not exon_added:
            new_line = (
                line.strip()
                + " exon_type {}; constitLevel {}; inclLevel {};".format(
                    "cassette", constitLevel, inclLevel
                )
            )
            kept_exons.append(new_line)
            exon_added = True

        elif overlapped and exon_added:
            print("!!ERROR!! both constitutive and cassette exon added")
            print(line_dict["start"], line_dict["end"])
            print(
                list(
                    cassette_interlap[line_dict["chrom"].replace("chr", "")].find(
                        (int(line_dict["start"]), int(line_dict["end"]))
                    )
                )
            )

    gtf_fh.close()

    print(
        "\n\tCreating new gtf file containing only constitutive and cassette exons: '{}'".format(
            args.out_file
        )
    )

    with open(args.out_file, "w") as out_fh:

        out_fh.write("#" + "\t".join(header) + "\n" + "\n".join(kept_exons))

    print("\nDONE\n")


if __name__ == "__main__":
    sys.exit(main() or 0)
