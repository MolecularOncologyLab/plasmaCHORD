import numpy as np
import pandas as pd 
import os 
import re
import subprocess
import matplotlib.pyplot as plt 
import pickle
import pyreadr
#from rpy2.robjects import r, pandas2ri



import argparse

def main():
    parser = argparse.ArgumentParser(description="Process an input variable.")
    parser.add_argument("input_var", type=int, help="An integer variable passed from the command line")
    args = parser.parse_args()
    
    #print("Received variable:", args.input_var)
    # Do additional processing here
    return args

if __name__ == "__main__":
    args = main()  # Capture the returned args

# Now you can use args outside of main()
#print("Outside main, input_var is:", args.input_var)


fake_id = pd.read_excel('BR36-PtsID-FakeIDS.xlsx')
fake_id_dic = {}
for index, row in fake_id.iterrows():
    fake_id_dic[row['Fake ID']] = row['ID']
df = pd.read_excel('BR36_mol_response_09202022.xlsx')
positions_df = pd.read_excel('Anagnostou_NatureMEdicine2023ctDNAresponse_supptables.xlsx',
                            sheet_name='ST6', skiprows=1)
positions_df = positions_df[positions_df['Cellular Origin']!='Germline']
column_dic = {'C1':'plasma.pgdx.id.C1', 'C2':'plasma.pgdx.id.C2', 'C3':'plasma.pgdx.id.C3'}
path_bam = '/dcs05/scharpf/data/jrwhite/az-download-2025-02/'



def SNV_search(start_pos, md_tag, variant_pos, ref_allele, alt_allele):
    """
    Determines if a read contains a specific variant based on the MD tag and other parameters.

    Parameters:
    - start_pos: int, starting position of the read in the reference genome.
    - md_tag: str, MD tag from the BAM file indicating mismatches.
    - variant_pos: int, position of the variant in the reference genome.
    - ref_allele: str, reference allele at the variant position.
    - alt_allele: str, alternative allele at the variant position.

    Returns:
    - True if the read contains the variant, False otherwise.
    """
    current_pos = start_pos
    mismatches = []  # List to store positions of mismatches

    # Process MD tag
    md_parts = md_tag.replace('MD:Z:', '').split('^')  # Remove MD:Z: and split deletions
    for part in md_parts:
        num = ''
        for char in part:
            if char.isdigit():
                num += char
            else:
                if num:
                    current_pos += int(num)
                    num = ''
                # For each mismatch, check if it matches the variant position and ref allele
                if current_pos == variant_pos and char == ref_allele:
                    return True  # Found the variant
                current_pos += 1  # Move to next position
        if num:  # Handle trailing numbers indicating matches
            current_pos += int(num)

    return False  # Variant not found in this read


import re

def has_deletion(start_pos, md_tag, cigar, variant_pos, deleted_alleles):
    
    """
    Determines if a read contains a specific variant based on the MD tag, CIGAR string, and other parameters.

    Parameters:
    - start_pos: int, starting position of the read in the reference genome.
    - md_tag: str, MD tag from the BAM file indicating mismatches and deletions.
    - cigar_string: str, CIGAR string indicating the alignment of the read to the reference.
    - variant_pos: int, position of the variant in the reference genome.
    - ref_alleles: str, reference bases at the variant position (for SNVs) or deleted (for deletions).

    Returns:
    - True if the read contains the variant, False otherwise.
    """
    # Check CIGAR string for deletion
    cigar_deletion = False
    current_pos = start_pos
    for length, op in re.findall(r'(\d+)([MID])', cigar):
        length = int(length)
        if op == 'M':
            current_pos += length
        elif op == 'D':
            if current_pos <= variant_pos < current_pos + length and length == len(deleted_alleles):
                cigar_deletion = True
                break
            current_pos += length
    
    # Check MD tag for deletion
    md_deletion = False
    md_pattern = re.compile(r'(\d+)|\^([ACGTN]+)')
    md_pos = start_pos
    for match in md_pattern.finditer(md_tag):
        if match.group(1):  # Matching bases
            md_pos += int(match.group(1))
        elif match.group(2):  # Deletion
            del_seq = match.group(2)
            if md_pos == variant_pos and del_seq == deleted_alleles:
                md_deletion = True
                break
    
    return cigar_deletion or md_deletion


def INS_search(start_pos, md_tag, cigar_string, variant_pos, insertion_seq):
    """
    Determines if a read contains a specific insertion based on the MD tag, CIGAR string, and other parameters.

    Parameters:
    - start_pos: int, starting position of the read in the reference genome.
    - md_tag: str, MD tag from the BAM file indicating mismatches and matches.
    - cigar_string: str, CIGAR string indicating the alignment of the read to the reference, including insertions.
    - variant_pos: int, position of the insertion in the reference genome.
    - insertion_seq: str, sequence of the insertion.

    Returns:
    - True if the read contains the insertion, False otherwise.
    """
    current_pos = start_pos
    insertion_detected_via_cigar = False

    # Parse CIGAR string
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    
    # Check for insertion using CIGAR string
    for count, op in cigar_ops:
        count = int(count)
        if op == 'I':  # Insertion to the reference
            # Check if the insertion is at the variant position
            if current_pos == variant_pos:
                insertion_detected_via_cigar = True
                break  # Found the insertion at the expected position
            # No position change in reference for insertions, but adjust for analysis context
        elif op in 'MDN=X':  # For match, deletion, skipped region, sequence match/mismatch
            current_pos += count  # Adjust current position

    # Use MD tag to support alignment verification, if necessary
    # For insertion detection, MD tag analysis might be less directly relevant
    # but could confirm the sequence integrity around the insertion point

    return insertion_detected_via_cigar
import re

def substitution_search(start_pos, md_tag, variant_pos, ref_allele, alt_allele):
    """
    Determines if a read contains a specific substitution (multi-nucleotide variant)
    based on the MD tag and other parameters.

    For example, if the variant is a substitution of ref "AA" to alt "CC",
    the MD tag is expected to show a block of mismatches corresponding to the
    reference bases (here, "AA") at the appropriate positions.

    Parameters:
    - start_pos: int, starting reference position of the read.
    - md_tag: str, the MD tag from the BAM file (e.g. "MD:Z:5AA10").
    - variant_pos: int, the starting reference position of the variant.
    - ref_allele: str, the reference allele sequence for the variant (e.g. "AA").
    - alt_allele: str, the alternate allele sequence (e.g. "CC").
                 (Note: This function uses the MD tag to check for the variant,
                 so alt_allele is provided for consistency but is not used here.)

    Returns:
    - True if the read contains the substitution (i.e. there is a block of mismatches
      at positions variant_pos to variant_pos+len(ref_allele)-1 that matches ref_allele),
      False otherwise.
    """
    # Remove any "MD:Z:" prefix
    md = md_tag.replace("MD:Z:", "")
    
    # We'll scan the MD tag and record the reference positions that are mismatched.
    # According to the MD tag format, digits indicate a number of matching bases,
    # letters indicate a mismatch (the letter is the reference base), and deletion
    # segments are prefixed by '^'. In this function we ignore deletion segments.
    pos = start_pos  # current reference position
    mismatch_dict = {}  # maps reference positions (int) -> reference base (str) reported as mismatch
    
    # Compile a regex that will match a number (match length), a single letter (mismatch),
    # or a deletion (e.g. '^AC').
    pattern = re.compile(r'(\d+)|([ACGTN])|(\^[ACGTN]+)')
    for m in pattern.finditer(md):
        if m.group(1):  # a block of matching bases
            match_length = int(m.group(1))
            pos += match_length
        elif m.group(2):  # a mismatch letter
            # Record that at this reference position, the read did not match;
            # the MD tag letter shows what the reference base is.
            mismatch_dict[pos] = m.group(2)
            pos += 1
        elif m.group(3):
            # A deletion segment (e.g. '^ACGT'): for our purposes (a substitution),
            # we ignore deletions. (Following the style in the deletion function,
            # we do not increment the reference position for deletions.)
            continue

    # Now check if the consecutive reference positions from variant_pos to
    # variant_pos+len(ref_allele)-1 are all mismatches and that their reported
    # reference bases (from the MD tag) match the expected ref_allele.
    observed_ref = ""
    for offset in range(len(ref_allele)):
        check_pos = variant_pos + offset
        if check_pos in mismatch_dict:
            observed_ref += mismatch_dict[check_pos]
        else:
            # If any base in the region is not a mismatch, then the substitution
            # is not observed in this read.
            return False

    # If the block of mismatches exactly equals the expected reference allele,
    # we conclude that the read supports the substitution event.
    if observed_ref == ref_allele:
        return True
    else:
        return False


def determine_variant_type(variant):
    """
    Determines the type of genetic variants (SNV, insertion, or deletion) based on the given list of variant strings.
    
    Parameters:
    - variants: list of str, variant strings formatted as 'chr_position-position_ref_alt'
    
    Returns:
    - List of tuples with variant description and its type (SNV, insertion, deletion).
    """

    ref, alt = variant.split('_')[-2], variant.split('_')[-1]
    if variant.split('_')[-1]=='':
        return 'deletion'
    elif (len(ref) ==1) and (len(alt) ==1):
        return 'SNV'
    #elif ref =="(null)":
    elif '__' in variant:
        return 'insertion'
    #elif alt =="(null)":
    #    return 'deletion'
    #elif len(alt) < len(ref):
    #    return 'sub-del'
    else:
        return 'substitution'
    
    
import re

def reverse_replacements(string2):
    # Replace the first '.' with ':'
    string1 = re.sub(r'\.', ':', string2, count=1)

    # Replace the second '.' with '>'
    string1 = re.sub(r'\.', '>', string1, count=1)

    # Replace '-' with '_'
    #string1 = string1.replace('-', '_')
    a = string1.split('--')[0]
    b = string1.split('--')[1].replace('-', '_')
    string1 = '__'.join([a,b])

    return string1

def get_fragment_bases(chromosome, start, end, reference_path="/Users/drabiza1/Documents/Reference/hg19.fa"):
    """
    Retrieves the bases for a specified fragment from a reference genome using samtools faidx.

    Parameters:
    - chromosome: str, the chromosome of the fragment.
    - start: int, the start position of the fragment.
    - end: int, the end position of the fragment.
    - reference_path: str, the path to the reference genome FASTA file.

    Returns:
    - str, the bases of the fragment, or an error message if the command fails.
    """
    region = f"{chromosome}:{start}-{end}"
    command = ["samtools", "faidx", reference_path, region]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        bases = ''.join(result.stdout.split('\n')[1:])  # Extract the bases from the output
        return bases
    except subprocess.CalledProcessError as e:
        return False
    
 
def calculate_fragment_coordinates(chromosome, start_pos, tlen):
    """
    Calculate the mapping coordinates (beginning and end) of a DNA fragment.

    Parameters:
    - chromosome: str, the chromosome where the fragment maps.
    - start_pos: int, the 1-based start position of the read on the chromosome.
    - tlen: int, the observed template length (TLEN) from the BAM file.

    Returns:
    - A tuple containing the chromosome, beginning, and end positions of the fragment.
    """
    if tlen > 0:
        # The read is the first in the pair and maps before its mate
        fragment_start = start_pos
        fragment_end = start_pos + tlen - 1
    else:
        # The read is the second in the pair and maps after its mate
        fragment_start = start_pos + tlen  # TLEN is negative
        fragment_end = start_pos - 1  # Adjust because the end position is inclusive

    return (chromosome.split('.fa')[0], fragment_start, fragment_end)


def fragment_info(sam_row):
    if sam_row[3] <= sam_row[7]:
        chrom, start, end = calculate_fragment_coordinates(sam_row[2],sam_row[3], sam_row[8])
    else:
        chrom, start, end = calculate_fragment_coordinates(sam_row[2],sam_row[7], abs(sam_row[8]) )
    return chrom, start, end, abs(sam_row[8])


def get_motifs(sam_row):
    if sam_row[3] <= sam_row[7]:
        chrom, start, end = calculate_fragment_coordinates(sam_row[2],sam_row[3], sam_row[8])
    else:
        chrom, start, end = calculate_fragment_coordinates(sam_row[2],sam_row[7], abs(sam_row[8]) )
    sequence = get_fragment_bases(chrom, start, end)
    return sequence[0:4].upper(), sequence[-4:].upper()
    
    
import os
#import concurrent.futures

#num_cores = os.cpu_count()  # Gets the total number of cores
#print("Number of cores:", num_cores)




#Mutation_Type = {'WBC_matched':'other', 'biopsy_matched':'mutant', 'IMPACT-BAM_matched':'mutant'}
#We are going to define all non-wild as mutant 

cols = ["id_chr_pos",  "start" , "end", "width", "mut_loc" ,"mutation", 
       "mut_read", "seq_read", "motif1" ,"motif2" ,"wild", "mutant"]
variant_reads = {}
cnt = 0 

cnter = 0 

#missing = []
still_missing = []
dir_name = 'BR36_Reads'

def Process_Sample(row):
        temp_df = []

        chrom = row['Nucleotide Position (Genomic hg19)'].split('_')[0]
        start = row['Nucleotide Position (Genomic hg19)'].split('_')[1].split('-')[0]

        output_file = row['Nucleotide Position (Genomic hg19)'].split('-')[0]
        reads_df = pd.read_csv(f"{dir_name}/{row['Patient ID']}.{row['Timepoint']}.{output_file}.sam",  sep ='\t', header = None, 
                    on_bad_lines='skip')

        vtype = determine_variant_type(row['Nucleotide Position (Genomic hg19)'])




        reads_df=reads_df[reads_df[5]!='*']
        total_reads = reads_df.shape[0]
        for i in range(9,reads_df.shape[1]):
            if sum(reads_df[i].str.contains('MD')) > 1000: 
                MD_index=int(i)


        vreads = []
        if vtype =='SNV':

            #reads_df = reads_df[reads_df[MD_index].str.contains('[ACGT]')]
            variant_pos = int(start)
            variant = row['Nucleotide Position (Genomic hg19)']
            ref_allele, alt_allele = variant.split('_')[-2], variant.split('_')[-1]
            #print(variant_pos, ref_allele)
            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                motif1, motif2 = get_motifs(row1)
                chrom, frag_start, frag_end , length = fragment_info(row1)
                if SNV_search(start_pos, md_tag, variant_pos, ref_allele, alt_allele) == True:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'mutant', 2 , alt_allele ,motif1, motif2, 
                                  ref_allele, alt_allele ])
                else:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'wild', 2 , ref_allele ,motif1, motif2, 
                                  ref_allele, alt_allele ])
            variant_reads = pd.DataFrame(temp_df, columns = cols)
            variant_reads.to_csv('Processed1/'+row['Patient ID']+'.'+row['Timepoint']+'.'+output_file+'_nofamily_galp_processed.rds.csv', 
             index = False)



        elif vtype =='substitution':

            #reads_df = reads_df[reads_df[MD_index].str.contains('[ACGT]')]
            variant_pos = int(start)
            variant = row['Nucleotide Position (Genomic hg19)']
            ref_allele, alt_allele = variant.split('_')[-2], variant.split('_')[-1]

            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                motif1, motif2 = get_motifs(row1)
                chrom, frag_start, frag_end , length = fragment_info(row1)
                if substitution_search(start_pos, md_tag, variant_pos, ref_allele, alt_allele) == True:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'mutant', 2 , alt_allele ,motif1, motif2, 
                                  ref_allele, alt_allele ])
                else:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'wild', 2 , ref_allele ,motif1, motif2, 
                                  ref_allele, alt_allele ])
            variant_reads = pd.DataFrame(temp_df, columns = cols)
            variant_reads.to_csv('Processed1/'+row['Patient ID']+'.'+row['Timepoint']+'.'+output_file+'_nofamily_galp_processed.rds.csv', 
             index = False)



        elif vtype =='insertion':

            insertion_seq = row['Nucleotide Position (Genomic hg19)'].split('_')[-1]
            variant_pos = int(start)


            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                cigar_string = row1[5]
                motif1, motif2 = get_motifs(row1)
                chrom, frag_start, frag_end , length = fragment_info(row1)
                if INS_search(start_pos, md_tag, cigar_string, variant_pos, insertion_seq) == True:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'mutant', 2 , insertion_seq ,motif1, motif2, 
                                  '-', insertion_seq ])
                else:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'wild', 2 , '-' ,motif1, motif2, 
                                  '-', insertion_seq ])
            variant_reads = pd.DataFrame(temp_df, columns = cols)
            variant_reads.to_csv('Processed1/'+row['Patient ID']+'.'+row['Timepoint']+'.'+output_file+'_nofamily_galp_processed.rds.csv', 
             index = False)

        elif vtype =='deletion':



            deleted_alleles=row['Nucleotide Position (Genomic hg19)'].split('_')[-2]
            variant_pos = int(start) #- len(deleted_alleles)


            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                cigar_string = row1[5]
                motif1, motif2 = get_motifs(row1)
                chrom, frag_start, frag_end , length = fragment_info(row1)
                if has_deletion(start_pos, md_tag, cigar_string, variant_pos, deleted_alleles)== True:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'mutant', 2 , '-' ,motif1, motif2, 
                                  deleted_alleles, '-' ])
                else:
                    temp_df.append([row['Nucleotide Position (Genomic hg19)'],frag_start, frag_end , length, 
                                    start, 'wild', 2 , deleted_alleles ,motif1, motif2, 
                                  deleted_alleles, '-' ])
            variant_reads = pd.DataFrame(temp_df, columns = cols)
            variant_reads.to_csv('Processed1/'+row['Patient ID']+'.'+row['Timepoint']+'.'+output_file+'_nofamily_galp_processed.rds.csv', 
             index = False)

#if __name__ == '__main__':
#    rows = [row for _, row in positions_df.iterrows()]
#    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
#        results = list(executor.map(Process_Sample, rows))
Process_Sample(positions_df.iloc[args.input_var])    
