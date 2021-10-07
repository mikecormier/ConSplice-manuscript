from __future__ import print_function

import argparse
import gzip
import io
import os
import sys
from collections import defaultdict

import numpy as np
from consplice.constraint.utils import get_alternative_gene_symbols
from interlap import InterLap

# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Global Parameters
# ------------------------------------------------------------------------------------------------------------------------------------------------------

INTERVAL_SIZE = 25

# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
# ------------------------------------------------------------------------------------------------------------------------------------------------------


def arguments():

    p = argparse.ArgumentParser(
        description="Score regions around exon features using the ConSplice regional percentile score."
    )

    req = p.add_argument_group("Required Arguments")

    req.add_argument(
        "--gtf",
        metavar="GRCh38 gtf file",
        required=True,
        help="(Required) The file path to a GRCh38 gtf file with exon features to extract.",
    )

    req.add_argument(
        "--exon-gtf",
        metavar="HEXEvent Exons GTF",
        required=True,
        help="(Required) File path to the constitutive and cassette exon gtf file created from combining a GTF file with the HEXEvent features. (Created from the 'constitutive_cassette_exon_gtf_filter.py' script)",
    )

    req.add_argument(
        "--constraint-scores",
        metavar="Regional ConSplice score file",
        required=True,
        help="(Required) The file path to the regional ConSplice score file. These scores will be added to the expanded range of constitutive and cassette exons",
    )

    req.add_argument(
        "--score-col",
        metavar="Constraint Score Column",
        required=True,
        help="(Required) The name of the column in the regional ConSplice score file to use for scoring regions from the exon-gtf",
    )

    req.add_argument(
        "--region-size",
        metavar="ConSplice region size",
        required=True,
        help="(Required) The ConSplice region size used to create regional scores",
    )

    req.add_argument(
        "--out-file",
        metavar="Output File",
        required="True",
        help="(Required) The file path and/or name of the output gtf file to create. This file will contain only constitutive and cassette exons features from the gtf file",
    )

    p.add_argument(
        "--extension-type",
        metavar="Exon extension type",
        choices=["neighboring-exon", "n-bp"],
        default="neighboring-exon",
        help="(Optional) Which extension technique to use. The type defines how far in each direction to get constraint scores from a target exon. Choices = ['neighboring-exon','n-bp']. Default = 'neighboring-exon'. 'neighboring-exon' means to get scores for the scores for all regions between the pre and proceeding exon. n-bp = the number of base pairs on each side of the target exon to use",
    )

    p.add_argument(
        "--flanking-bp",
        metavar="Flanking basepairs",
        default=1000,
        type=int,
        help="(Optional) The number of flanking basepairs to get scores for from the start and end of the target exon. Must be positive. If the neighbor-exon extension type is used then this parameter will be used for first and last exons only. Default = 1000",
    )

    p.add_argument(
        "--score-chrom-col",
        metavar="Score chrom column",
        default="chrom",
        help="(Optional) The name of the chromosome column in the constraint score file. Default = 'chrom'",
    )

    p.add_argument(
        "--score-start-col",
        metavar="Score start pos column",
        default="region_start",
        help="(Optional) The name of the start position column in the constraint score file. Default = 'region_start'",
    )

    p.add_argument(
        "--score-end-col",
        metavar="Score end pos column",
        default="region_end",
        help="(Optional) The name of the end position column in the constraint score file. Default = 'region_start'",
    )

    return p.parse_args()


# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Functions
# ------------------------------------------------------------------------------------------------------------------------------------------------------


def extract_gtf_exon_info(file_path):
    """
    extract_gtf_exon_info
    =====================
    Method to extract all exon info from a gtf file. All exon features will be extracted into a 2d list. Exons will be sorted by start
     and end position. All exons for a gene will be combined into a single 2d list.

    Parameters:
    -----------
    1) file_path: (str) The file path to the gtf file to extract exon info from

    Return:
    ++++++
    1) (dict) A dictionary of 2d lists. keys = gene names, values = 2d list of exons.
                2d list indices: 0 = start pos, 1 = end pos

    """

    exons_by_gene = defaultdict(list)

    try:
        gtf_fh = (
            gzip.open(file_path, "rt", encoding="utf-8")
            if file_path.endswith(".gz")
            else io.open(file_path, "rt", encoding="utf-8")
        )
    except IOError as e:
        print("\n!!ERROR!! unable to open the gtf file: '{}'".format(file_path))
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

        ## Add exon feature to dict
        exons_by_gene[gene_name].append(
            [
                int(line_dict["start"]),
                int(line_dict["end"]),
            ]
        )
    gtf_fh.close()

    ## sort list by genomic position
    for gene in exons_by_gene:
        ## sort by exon start position and then by exon end position.
        exons_by_gene[gene] = sorted(exons_by_gene[gene], key=lambda x: (x[0], x[1]))

    return exons_by_gene


def extract_constraint_score(file_path, chrom_col, start_col, end_col, score_col):
    """
    extract_constraint_score
    ========================
    This method is used to create interlap objects for the splicing constraint scores for fast interval tree lookup.

    Parameters:
    -----------
    1) file_path: (str) The file path to the splicing constraint regional scores file
    2) chrom_col: (str) The name of the chromosome column in the constraint file
    3) start_col: (str) The name of the start column in the constraint file
    4) end_col:   (str) The name of the end column in the constraint file
    5) score_col: (str) The name of the score column in the constraint file

    Return:
    +++++++
    1) (dict) A dictionary of interlap objects. Keys = chromosome name with any 'chr' prefix removed.
                value = the interlap object for constraint regions. (Index 0 = region start, 1 = region end, 2 = region score, 3 = gene name)
    """

    score_dict = defaultdict(InterLap)

    try:
        fh = (
            gzip.open(file_path, "rt", encoding="utf-8")
            if file_path.endswith(".gz")
            else io.open(file_path, "rt", encoding="utf-8")
        )
    except IOError as e:
        print("\n!!ERROR!! unable to open the gtf file: '{}'".format(file_path))
        print(str(e))
        sys.exit(1)

    header = fh.readline().strip().replace("#", "").split("\t")

    ## Add each score to an interlap object
    for line in fh:

        line_dict = dict(zip(header, line.strip().split("\t")))

        score_dict[line_dict[chrom_col].replace("chr", "")].add(
            (
                int(line_dict[start_col]),
                int(line_dict[end_col]),
                float(line_dict[score_col].strip()),
                line_dict["gene_symbol"],
            )
        )

    fh.close()

    return score_dict


def get_exon_expansion(target_exon, exon_list, expansion_type, flanking_bp):
    """
    get_exon_expansion
    ==================
    This method is used to get the left and right expansion of the target exon based on an expansion type.
     The expanded region will include the entire target exon with flanking space on both the 5' and 3'
     end of the target exon. This function does not take into account positive or negative strand, but
     rather provides the same flanking on either side of an exon.


     If the exon is either the first or last exon of a gene and the expansion_type == 'neighboring-exon', then:
        - if first exon, the start expansion will use the flanking_bp parameter and the end expansion will use the next exon
        - if last exon, the start expansion will use the previous exon and the end expansion will use the flanking_bp parameter

    Parameters:
    -----------
    1) target_exon:   (dict) A dictionary with gtf information about the target/principle exon
    2) exon_list:     (list) A 2d list of all sorted exons for the gene where the target exon is found
    3) epansion_type: (str)  Which type of expansion method to use. (Neighboring exons or number of bp)
    4) flanking_bp:   (int)  The number of base pairs to add to the start and end of the target exon.
                              (Only used if expansion_type = 'n-bp')

    Returns:
    ++++++++
    1) (int) The expansion start position (Strand ignored. Meaning, this may either be the start of the expansion region for + strand or the end of the region for - strand)
    2) (int) The expansion end position   (strand ignored. Meaning, this may either be the end of the expansion region for + strand for the start of the region for - strand)
    """

    target_start = int(target_exon["start"])
    target_end = int(target_exon["end"])

    exon_index = -1

    ## Iterate through the list of exons for a gene and identify the index for the target exon
    for i, exon in enumerate(exon_list):

        if exon[0] == target_start and exon[1] == target_end:
            exon_index = i
            break

    ## Check that the exon was found
    if exon_index == -1:
        print("\n!!ERROR!! Unable to identify the correct exon for the gene")
        sys.exit(1)

    range_start = -1
    range_end = -1

    ## Get the expansion start and end
    if expansion_type == "neighboring-exon":

        range_start = (
            exon_list[exon_index - 1][0]
            if exon_index != 0
            else exon_list[exon_index][0] - int(flanking_bp)
        )
        range_end = (
            exon_list[exon_index + 1][1]
            if exon_index != (len(exon_list) - 1)
            else exon_list[exon_index][1] + int(flanking_bp)
        )

    elif expansion_type == "n-bp":

        range_start = exon_list[exon_index][0] - int(flanking_bp)
        range_end = exon_list[exon_index][1] + int(flanking_bp)

    ## Check for bad expansion type
    else:
        print(
            "\n!!ERROR!! Incorrect expansion type designated: '{}'".format(
                expansion_type
            )
        )
        sys.exit(1)

    ## Check that the expansion ranges is correct
    if range_start == -1 or range_end == -1:
        print("\n!!ERROR!! Unable to identify correct expansion range")
        sys.exit(1)
    elif range_start == target_start and range_end == target_end:
        print("\n!!ERROR!! Unable to identify correct expansion range")
        sys.exit(1)

    return (range_start, range_end)


def identify_flanking_score(
    interlap_scores,
    target_start_score,
    target_end_score,
    left_start,
    right_end,
    left_region_count,
    center_region_count,
    right_region_count,
    gene_name,
    region_size,
):
    """
    identify_flanking_score
    ========================
    This method is used to get the regional constraint score for 3 different set of input regions. For each region it will identify the
     region scores based on regional constraint scores within that region and sorted by increasing genomic position on the positive strand.
     3 separate score lists are created:
        1) The first represents the 5' upstream scores from the target exon. These will be index specific, meaning they are sorted by position
        2) The second represents the scores that overlap the target exon.
        3) The third represents the 3' downstream scores from the target exon. These are also index specific

    NOTE: THIS METHOD DOES NOT CORRECT FOR STRAND

    [5' list]     [target exon list]       [3' list]
    upstream ------- target exon -------- downstream

    Parameters:
    -----------
    1) interlap_score:   (InterLap) A chromosome specific interlap object that contains regional splicing constraint for the
                                   chromosome for which the gene is on.
    2) target_start_score:  (tuple) A tuple from the interlap object that represents the target exon start splicing constraint region
    3) target_end_score:    (tuple) A tuple from the interlap object that represents the target exon end splicing constraint region
    4) left_start:          (tuple) A tuple from the interlap object that represents the left-most start splicing constraint region to look at
    5) right_end:           (tuple) A tuple from the interlap object that represents the right-most end splicing constraint region to look at
    6) left_region_count:   (int)   The number of splicing constraint regions to the left of the target exon
    7) center_region_count: (int)   The number of splicing constraint regions in the target exon
    8) right_region_count:  (int)   The number of splicing constraint regions to the right of the target exon
    9) gene_name:           (str)   The name of the gene which the target exon is in
    10) region_size:        (int)   The size of the splicing constraint score regions

    Returns:
    ++++++++
    1) (numpy array) splicing constraint scores for all regions to the left of the target exon
    2) (numpy array) splicing constraint scores for all regions that overlap the target exon
    3) (numpy array) splicing constraint scores for all regions to the right of the target exon
    """

    def get_index(score_pos, end_pos, region_size, region_count):
        """
        get_index
        =========
        Get the array index of a region based on region size, total number or regions in the array,
         and the exact position of the region.

        Parameters:
        -----------
        1) score_pos: (int) The start position of the current score region
        2) end_pos:   (int) The end position for the entire region
        3) region_size: (int) The size of the score regions
        4) region_count: (int) The number of total regions for a predetermined range

        Returns:
        ++++++++
        (int) The 0 based index of where that region would exist in a list/array
        """

        index = int(
            int(region_count) - ((int(end_pos) - int(score_pos)) / int(region_size))
        )

        return index

    ## Set up empty numpy arrays with the size of each section
    left_scores = np.empty(int(left_region_count))
    left_scores[:] = np.nan

    center_scores = np.empty(int(center_region_count))
    center_scores[:] = np.nan

    right_scores = np.empty(int(right_region_count))
    right_scores[:] = np.nan

    ## Region = left start <-> center start - 1
    if left_region_count >= 1:
        for score in interlap_scores.find((left_start[0], target_start_score[0] - 1)):

            ## Check that the gene name matches:
            if score[3] != gene_name:
                continue

            ## Get the array index
            array_index = get_index(
                score[0], target_start_score[0], region_size, left_region_count
            )

            if array_index < 0:
                print("Left")
                print("Below Zero")

            ## Add score to numpy array
            left_scores[array_index] = score[2]

    ## Region = center start <-> center end
    if center_region_count >= 1:
        for score in interlap_scores.find((target_start_score[0], target_end_score[1])):

            ## Check that the gene name matches:
            if score[3] != gene_name:
                continue

            ## Get the array index
            array_index = get_index(
                score[0], target_end_score[1] + 1, region_size, center_region_count
            )

            if array_index < 0:
                print("Center")
                print("Below Zero")

            ## Add score to numpy array
            center_scores[array_index] = score[2]

    ## Region = center_end + 1 <->  right end
    if right_region_count >= 1:
        for score in interlap_scores.find((target_end_score[1] + 1, right_end[1])):

            ## Check that the gene name matches:
            if score[3] != gene_name:
                continue

            ## Get the array index
            array_index = get_index(
                score[0], right_end[1] + 1, region_size, right_region_count
            )

            if array_index < 0:
                print("right")
                print("Below Zero")

            ## Add score to numpy array
            right_scores[array_index] = score[2]

    return (left_scores, center_scores, right_scores)


def get_interval_scores(
    target_exon,
    left_region_start,
    right_region_end,
    interlap_scores,
    gene_name,
    interval_size,
):

    """
    get_interval_scores
    ===================
    This method is to get regional constraint scores from left to right by a set interval. This method uses the exon start and end positions
     to identify which genomic positions to use for regional score lookup. The intervals will extend until the region start position is reached
     and the region end position is reached.

    This method will create 5 lists.
        1) The Left to Right  oriented (5' -> 3' for + strand) list of region scores up to the exon start position separated by interval size
        2) The region score for the exon start position
        3) The Left to Right  oriented (5' -> 3' for + strand) list of region scores for 4 regions in the exon
        4) The region score for the exon end position
        5) The Left to Right  oriented (5' -> 3' for + strand) list of region scores from the exon end position to the region end position separated by interval size

    Parameters:
    -----------
    1) target_exon:        (dict) A dictionary that represents a gtf entry for this exon
    2) left_region_start:   (int) The left most genomic position of the region of interest
    3) right_region_start:  (int) The right most genomic position of the region of interest
    4) interlap_score: (InterLap) A chromosome specific interlap object that contains regional splicing constraint for the
                                   chromosome for which the gene is on.
    5) gene_name:           (str)   The name of the gene which the target exon is in
    6) interval_size:       (int) The interval size to use to separate entries by

    Returns:
    ++++++++
    1) (list) List of region scores separated by interval size left of the exon start
    2) (list) List of region scores for the exon start position (Should be n of 1)
    3) (list) List of region scores for intervals inside the exon
    4) (list) List of region scores for the exon end position (Should be n of 1)
    5) (list) List of region scores separated by interval size right of the exon end

    """

    exon_start = int(target_exon["start"])
    exon_end = int(target_exon["end"])

    ## Get left scores separated by the interval size
    temp_pos = exon_start - interval_size
    left_scores = []
    index = 0
    while temp_pos > left_region_start:

        ## get the scores for this position
        temp_scores = []
        for score in interlap_scores.find((temp_pos, temp_pos)):

            ## Check that the gene name matches:
            if score[3] != gene_name:
                continue

            temp_scores.append(score[2])

        if len(temp_scores) > 1:
            print("!!ERROR!! Too many left scores")

        ## If there is no region score, stop. don't continue if there is no region or a missing region
        if len(temp_scores) < 1:
            break

        ## Add score to list
        left_scores.append(temp_scores[0])

        ## Check the index
        assert len(left_scores) - 1 == index

        ## Update temp_pos
        temp_pos -= interval_size
        index += 1

    ## Need to flip the left list so it is correctly ordered left to right
    left_scores.reverse()

    ## Get right scores separated by the interval size
    temp_pos = exon_end + interval_size
    right_scores = []
    index = 0
    while temp_pos < right_region_end:

        ## get the scores for this position
        temp_scores = []
        for score in interlap_scores.find((temp_pos, temp_pos)):

            ## Check that the gene name matches:
            if score[3] != gene_name:
                continue

            temp_scores.append(score[2])

        if len(temp_scores) > 1:
            print("!!ERROR!! Too many right scores")

        ## If there is no region score, stop. don't continue if there is no region or a missing region
        if len(temp_scores) < 1:
            break

        right_scores.append(temp_scores[0])

        ## Check the index
        assert len(right_scores) - 1 == index

        ## Update temp_pos
        temp_pos += interval_size
        index += 1

    ## Get the exon start score
    left_exon_score = []
    for score in interlap_scores.find((exon_start, exon_start)):
        ## Check that the gene name matches:
        if score[3] != gene_name:
            continue

        left_exon_score.append(score[2])

    if len(left_exon_score) > 1:
        print("!!ERROR!! To many start scores")

    ## Get the exon end score
    right_exon_score = []
    for score in interlap_scores.find((exon_end, exon_end)):
        ## Check that the gene name matches:
        if score[3] != gene_name:
            continue

        right_exon_score.append(score[2])

    if len(right_exon_score) > 1:
        print("!!ERROR!! To many end scores")

    ## Center Scores
    exon_scores = []
    for pos in [
        exon_start + interval_size,
        exon_start + (interval_size * 2),
        exon_end - (interval_size * 2),
        exon_end - interval_size,
    ]:
        temp_scores = []
        for score in interlap_scores.find((pos, pos)):
            ## Check that the gene name matches:
            if score[3] != gene_name:
                continue

            temp_scores.append(score[2])

        if len(temp_scores) > 1:
            print("!!ERROR!! To many center scores")

        ## If there is no region score, stop. don't continue if there is no region or a missing region
        if len(temp_scores) < 1:
            break

        exon_scores.append(temp_scores[0])

    return (left_scores, left_exon_score, exon_scores, right_exon_score, right_scores)


def score_regions(
    region_start,
    region_end,
    target_exon,
    interlap_scores,
    region_size,
    gene_name,
    interval_size,
):
    """
    score_regions
    =============
    This method is used to get the regional constraint score regions between a start and end position.
     8 separate score lists are created:
        1) The first represents the 5' upstream scores from the target exon. These will be index specific, meaning they are sorted by position
        2) The second represents the scores that overlap the target exon.
        3) The third represents the 3' downstream scores from the target exon. These are also index specific
        4) The fourth represents the constraint scores at different 25x basepair intervals from the 5' side of the exon. That is index -1 == -25 bases from the 5' end of the exon, -2 == -50, -3 = -75, ..., -10 == -250. etc.
        5) The fifth represents the constraint score at the 5' end of the intron. (Relative exon 5' start position == 0)
        6) The sixth represents the splicing constraint scores in the target exon at different 25x basepair intervals. That is, index 0 == 25 bases from the 5' end of the exon, 1 = 50 bases from the 5' end, 2 = -50 from 3' end, and 3 = -25 from 3' end.
        7) The seventh represents the constraint score at the 3' end of the intron. (Relative exon 3' start position == 0)
        8) The eighth represents the constraint scores at different 25x basepair intervals from the 3' side of the exon. That is index 0 == 25 bases from the 3' end of the exon, 1 == 50, 3 = 75, ..., 10 == 250. etc.

    Lists 4-8 use a 25x basepair interval.

    All lists are sorted by position from 5' to 3'. Scores on genes on the negative strand are are re-sorted to
     give the scores from 5' to 3' so that the positive and negative strand scores can be compared.

    [5' list]     [target exon list]       [3' list]
    upstream ------- target exon -------- downstream

    regional scores are identified in using the identify_flanking_score() function, which returns a left, center, and right list of scores. (NOTE: These are not strand specific. This function updates the data to strand specific positions)
    25x interval scores are identified using the get_interval_scores() function, which returns the left, exon start, exon center, exon end, and right interval lists. (NOTE: These are not strand specific. This function updates the data to strand specific positions)

    In this method, the 5' and/or the 3' regions are adjusted to the most possible regions given the total region parameters when either one
     has missing regions. That is, if the region start or region end position is outside of available score regions, the start and/or
     end positions are adjusted to the closest score region to the start or end position.


    Parameters:
    -----------
    1) region_start:   (int) The start position of the region
    2) region_end:     (int) The end position of the region
    3) target_exon:    (dict) A dictionary with gtf information about the target/principle exon
    4) interlap_score: (dict) A by chromosome interlap object with regional splicing constraint scores
    5) region_size:    (int)  The size of the splicing constraint score regions
    6) gene_name:      (str)  The name of the gene which the target exon is in
    7) interval_size:  (int)  The size of the interval to use to get scores

    Returns:
    ++++++++
    1) (numpy array) position sorted 5' upstream scores (scores upstream of the target exon)
    2) (numpy array) position sorted center scores (scores for the target exon)
    3) (numpy array) position sorted 3' downstream scores (scores downstream of the target exon)
    4) (list)        position sorted 5' upstream scores by interval up to exon start
    5) (list)        n of 1 list score for exon start position
    6) (list)        position sorted 5' -> 3' scores within the target exon
    7) (list)        n of 1 list score for exon end position
    8) (list)        position sorted 3' downstream scores by interval from exon end position
    """

    chrom = target_exon["chrom"].replace("chr", "")

    region_size = int(region_size)

    ## Get the start and end score positions for the entire region and the target exon
    start_scores = list(
        interlap_scores[chrom].find(
            (int(target_exon["start"]), int(target_exon["start"]))
        )
    )
    end_scores = list(
        interlap_scores[chrom].find((int(target_exon["end"]), int(target_exon["end"])))
    )

    ## Check for missing scores
    if not start_scores or not end_scores:
        return (["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"])

    target_start_score = []
    ## Get the gene specific score info
    for score in start_scores:
        if score[3] == gene_name:
            target_start_score = score
            break

    target_end_score = []
    for score in end_scores:
        if score[3] == gene_name:
            target_end_score = score
            break

    ## Check for missing scores
    if not target_end_score or not target_end_score:
        return (["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"])

    ## Get the left and right scores
    temp_start_pos = int(region_start)
    left_start = target_start_score

    ## Iterate over each possible region upstream of the target start position to find the
    ### most possible regions given the range parameters
    while temp_start_pos < int(target_exon["start"]):

        left_scores = list(
            interlap_scores[chrom].find((int(temp_start_pos), int(temp_start_pos)))
        )

        if left_scores:
            for score in left_scores:
                if score[3] == gene_name:
                    left_start = score
            break
        else:
            temp_start_pos += region_size

    ## Get the left and right scores
    temp_end_pos = int(region_end)
    right_end = target_end_score

    ## Iterate over each possible region downstream of the target end position to find the
    ### most possible regions given the range parameters
    while temp_end_pos > int(target_exon["end"]):

        right_scores = list(
            interlap_scores[chrom].find((int(temp_end_pos), int(temp_end_pos)))
        )

        if right_scores:
            for score in right_scores:
                if score[3] == gene_name:
                    right_end = score
            break
        else:
            temp_end_pos -= region_size

    ## Get the number of regions for each section
    left_region_count = (target_start_score[0] - left_start[0]) / region_size

    center_region_count = (
        (target_end_score[1] + 1) - target_start_score[0]
    ) / region_size

    right_region_count = (right_end[1] - target_end_score[1]) / region_size

    ## Get the left, center, and right region scores
    left_scores, center_scores, right_scores = identify_flanking_score(
        interlap_scores[chrom],
        target_start_score,
        target_end_score,
        left_start,
        right_end,
        left_region_count,
        center_region_count,
        right_region_count,
        gene_name,
        region_size,
    )

    ## Get scores by interval
    (
        left_interval_scores,
        left_exon_score,
        exon_interval_scores,
        right_exon_score,
        right_interval_scores,
    ) = get_interval_scores(
        target_exon,
        region_start,
        region_end,
        interlap_scores[chrom],
        gene_name,
        interval_size,
    )

    ## Fix negative strand regions. (Flip 5' and 3')
    if target_exon["strand"] == "-":

        ## flip all arrays so they are oriented 5' to 3' from left to right
        ## Change the "right_score" to the left score position and the "left_score" to the right score position
        ##
        ##     5'                                                           3'

        right_interval_scores.reverse()
        right_exon_score.reverse()
        exon_interval_scores.reverse()
        left_exon_score.reverse()
        left_interval_scores.reverse()

        return (
            np.flip(right_scores) if right_scores.size > 0 else ["NA"],
            np.flip(center_scores) if center_scores.size > 0 else ["NA"],
            np.flip(left_scores) if left_scores.size > 0 else ["NA"],
            right_interval_scores if len(right_interval_scores) > 0 else ["NA"],
            right_exon_score if len(right_exon_score) > 0 else ["NA"],
            exon_interval_scores if len(exon_interval_scores) > 0 else ["NA"],
            left_exon_score if len(left_exon_score) > 0 else ["NA"],
            left_interval_scores if len(left_interval_scores) > 0 else ["NA"],
        )

    else:
        ## Return the left to right 5' to 3' positive strand scores
        ##
        ###    5'                                               3'
        return (
            left_scores if left_scores.size > 0 else ["NA"],
            center_scores if center_scores.size > 0 else ["NA"],
            right_scores if right_scores.size > 0 else ["NA"],
            left_interval_scores if len(left_interval_scores) > 0 else ["NA"],
            left_exon_score if len(left_exon_score) > 0 else ["NA"],
            exon_interval_scores if len(exon_interval_scores) > 0 else ["NA"],
            right_exon_score if len(right_exon_score) > 0 else ["NA"],
            right_interval_scores if len(right_interval_scores) > 0 else ["NA"],
        )


# ------------------------------------------------------------------------------------------------------------------------------------------------------
## Main
# ------------------------------------------------------------------------------------------------------------------------------------------------------


def main():

    args = arguments()

    assert (
        args.flanking_bp > 0
    ), "\n!!ERROR!! The flanking basepair input parameter must be a positive value\n"

    print(
        (
            "\nInput Arguments:"
            "\n================"
            "\n - gtf:               {}"
            "\n - exon-gtf:          {}"
            "\n - constraint-scores: {}"
            "\n - score-col:         {}"
            "\n - region-size:       {}"
            "\n - out-file:          {}"
            "\n - extension-type:    {}"
            "\n - flanking-bp:       {}bp ({})"
            "\n - score-chrom-col:   {}"
            "\n - score-start-col:   {}"
            "\n - score-end-col:     {}"
        ).format(
            args.gtf,
            args.exon_gtf,
            args.constraint_scores,
            args.score_col,
            args.region_size,
            args.out_file,
            args.extension_type,
            args.flanking_bp,
            "Used for all exons"
            if args.extension_type == "n-bp"
            else "Used only for first and last exons",
            args.score_chrom_col,
            args.score_start_col,
            args.score_end_col,
        )
    )

    print("\nProcessing")
    print("==========")

    print("\n\tExtracting exon info from the gtf file")
    exons_by_gene = extract_gtf_exon_info(args.gtf)

    print("\n\tCreating an interval tree for fast regional splicing constraint lookup")
    scores_interlap = extract_constraint_score(
        args.constraint_scores,
        args.score_chrom_col,
        args.score_start_col,
        args.score_end_col,
        args.score_col,
    )

    print("\n\tParsing constitutive and cassette exon gtf file.")
    print("\n\tAdding constraint scores to the expanded range")
    try:
        gtf_fh = (
            gzip.open(args.exon_gtf, "rt", encoding="utf-8")
            if args.exon_gtf.endswith(".gz")
            else io.open(args.exon_gtf, "rt", encoding="utf-8")
        )
    except IOError as e:
        print("\n!!ERROR!! unable to open the gtf file: '{}'".format(args.exon_gtf))
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

    new_lines = []

    for j, line in enumerate(gtf_fh):

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

        gene_name = (
            line_dict["gene_name"]
            if "gene_name" in line_dict
            else line_dict["gene_symbol"]
        )

        ## Get expansion positions
        expansion_start, expansion_end = get_exon_expansion(
            line_dict, exons_by_gene[gene_name], args.extension_type, args.flanking_bp
        )

        ## Get regional scores for exon
        (
            upstream_list,
            center_list,
            downstream_list,
            upstream_interval_list,
            exon_start_list,
            exon_interval_list,
            exon_end_list,
            downstream_interval_list,
        ) = score_regions(
            expansion_start,
            expansion_end,
            line_dict,
            scores_interlap,
            args.region_size,
            gene_name,
            INTERVAL_SIZE,
        )

        ## create new output line
        new_line = [
            line_dict["chrom"],
            line_dict["start"],
            line_dict["end"],
            line_dict["strand"],
            line_dict["gene_name"],
            line_dict["gene_id"],
            line_dict["transcript_id"],
            line_dict["gene_type"],
            line_dict["exon_type"],
            line_dict["constitLevel"],
            line_dict["inclLevel"],
            expansion_start,
            expansion_end,
            ",".join(map(str, upstream_list)),
            ",".join(map(str, center_list)),
            ",".join(map(str, downstream_list)),
            ",".join(map(str, upstream_interval_list)),
            ",".join(map(str, exon_start_list)),
            ",".join(map(str, exon_interval_list)),
            ",".join(map(str, exon_end_list)),
            ",".join(map(str, downstream_interval_list)),
        ]

        new_lines.append("\t".join(map(str, new_line)))

    gtf_fh.close()

    ## Create output file
    print("\n\tCreating output file: '{}'".format(args.out_file))

    output_header = [
        "#chrom",
        "start",
        "end",
        "strand",
        "gene_name",
        "gene_id",
        "transcript_id",
        "gene_type",
        "exon_type",
        "constitLevel",
        "inclLevel",
        "expansion_region_start",
        "expansion_region_end",
        "Five_Prime_Scores",
        "Center_Scores",
        "Three_Prime_Scores",
        "Five_Prime_Interval_Scores",
        "Five_Prime_Exon_start",
        "Exon_Interval_Scores",
        "Three_Prime_Exon_End",
        "Three_Prime_Interval_Scores",
    ]

    with open(args.out_file, "w") as out:
        out.write("## Expansion Type: {}\n".format(args.extension_type))
        out.write("## Flanking bp: {}\n".format(args.flanking_bp))
        out.write("## Scores from: {}\n".format(args.constraint_scores))
        out.write("## Score region size: {}\n".format(args.region_size))
        out.write("## Interval Size: {}\n".format(INTERVAL_SIZE))
        out.write("## All scores are oriented 5' to 3'\n")
        out.write("\t".join(output_header) + "\n")
        out.write("\n".join(new_lines))

    print("\nDONE\n")


if __name__ == "__main__":
    sys.exit(main() or 0)
