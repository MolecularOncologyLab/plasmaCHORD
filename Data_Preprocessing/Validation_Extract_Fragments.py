#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('pylab', 'inline')
import pandas as pd 
import os 
import re
import subprocess
import matplotlib.pyplot as plt 
import pickle
import pyreadr


# # Import validation

# In[28]:


#import variant excel
variant_df= pd.read_csv('variant_data.csv')


# In[3]:


#FIlter to remove other type of variants 
#variant_df = variant_df[variant_df['Origin']=='Tumor']


# In[29]:


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



# In[30]:


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



# In[31]:


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


# In[32]:


def determine_variant_type(variant):
    """
    Determines the type of genetic variants (SNV, insertion, or deletion) based on the given list of variant strings.
    
    Parameters:
    - variants: list of str, variant strings formatted as 'chr_position-position_ref_alt'
    
    Returns:
    - List of tuples with variant description and its type (SNV, insertion, deletion).
    """

    
    if '__' in variant:
        variant_type = 'insertion'
    elif variant.endswith('_'):
        variant_type = 'deletion'
    else:
        variant_type = 'SNV'
    
    
    return variant_type


# In[1]:


#df = pd.read_csv(f'{dir_name}/{output_file}.sam',  sep ='\t', header = None, 
#                error_bad_lines=False, warn_bad_lines=False)


# ## Use Samtools to extract fragments at each position
# .. 
# ## Once Samtools was run, import the fragments 

# In[ ]:


dir_name = 'Validation-Reads/'


# In[ ]:


variant_reads = {}
cnt = 0 

missing = []
dir_name = 'Validation-Reads/'

#FOR NOW remove TNP and DNP 
variant_df[variant_df['Variant_Type']!='TNP']
variant_df[variant_df['Variant_Type']!='DNP']

for ind, row in variant_df.iterrows():
    processed_df = pd.DataFrame()
    
    
    try:
    

        chrom = row['Chr']
        start = row['Start_hg19']
        path = row['bamlink']
        output_file = row["index"].replace('_','-').replace('>','.').replace(':', '.')


        reads_df = pd.read_csv(f'{dir_name}/{output_file}.sam',  sep ='\t', header = None, 
                    error_bad_lines=False, warn_bad_lines=False)
        vtype = row['Variant_Type']
        reads_df=reads_df[reads_df[5]!='*']

        #Find the MD tage column 
        for i in range(9,reads_df.shape[1]):
            if sum(reads_df[i].str.contains('MD')) > 100: 
                MD_index=int(i)


        vreads = []
        temp_df = []
        if vtype =='SNP':

            reads_df = reads_df[reads_df[MD_index].str.contains('[ACGT]')]
            variant_pos = int(row['Start_hg19'])
            ref_allele = row['Base_from']
            alt_allele = row['Base_to']
            #print(variant_pos, ref_allele)
            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                if SNV_search(start_pos, md_tag, variant_pos, ref_allele, alt_allele) == True:
                    #vreads.append(row1[0])
                    temp_df.append(row['index'], )
            variant_reads[ind] = reads_df[reads_df[0].isin(vreads)]

        elif vtype =='INS':
            reads_df = reads_df[reads_df[5].str.contains('[ISD]')]
            variant_pos = int(row['End_hg19'])
            #IF INSERTIONS are problematic change 
            insertion_seq = row['vars'].split('>')[-1]

            #print(variant_pos, ref_allele)
            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                cigar_string = row1[5]
                if INS_search(start_pos, md_tag, cigar_string, variant_pos, insertion_seq) == True:
                    vreads.append(row1[0])
            variant_reads[ind] = reads_df[reads_df[0].isin(vreads)]

        elif vtype =='DEL':
            reads_df = reads_df[reads_df[5].str.contains('[ID]')]

             #IF DELETIONS are problematic change 
            deleted_alleles=row['Base_from']
            #deleted_alleles = row['vars'].split('_')[-1].split('>')[0]
            variant_pos = int(row['End_hg19']) - len(deleted_alleles)
            #variant_pos = int(row['Start_hg19'])
            #variant_pos = int(row['vars'].split(':')[-1].split('_')[0])+2





            #print(variant_pos, ref_allele)
            for ind1, row1 in reads_df.iterrows():
                start_pos = row1[3]
                md_tag = row1[MD_index]
                cigar_string = row1[5]
                if has_deletion(start_pos, md_tag, cigar_string, variant_pos, deleted_alleles)== True:
                    vreads.append(row1[0])      
            variant_reads[ind] = reads_df[reads_df[0].isin(vreads)]

    except:
        missing.append(output_file)
            #print(SNV_search(start_pos, md_tag, variant_pos, ref_allele, alt_allele))
            #if 'C' in md_tag:
             #   print(start_pos, md_tag)
    cnt+=1
    #print(len(vreads),row['Variant_Type'] , row['MAF'])
    #if len(vreads)==0:
    #        print(row['Variant_Type'])
    


# In[136]:


for i in variant_reads.keys():
    print(variant_reads[i].shape[0], variant_df.loc[i]['distinct.mutant.reads'])


# In[143]:


with open('variant_reads.pkl', 'wb') as f:
    pickle.dump(variant_reads, f)


# In[ ]:


with open('variant_reads.pkl', 'rb') as f:
    variant_reads = pickle.load(f)


# In[ ]:




