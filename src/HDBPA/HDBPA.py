import argparse
import pandas as pd
import numpy as np
import re
import pysam
from pysam import FastaFile
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm, collections
from tqdm import tqdm

EOL = '\n'
    
def HDBPA(args):
    if args.mode not in ['NT', 'AA']:
        print('Unrecognized mode.')
        sys.exit(1)

    pop_size_1 = reform_aligned_reads(args.SAM1, args.r, args.mode)
    pop_size_2 = reform_aligned_reads(args.SAM2, args.f, args.mode)  
    
    df_phylo = phylogeny_all(args.SAM1.split('.sam')[0] + '_reformed.csv', 
                  args.SAM2.split('.sam')[0] + '_reformed.csv')
    
    df_phylo.to_csv(os.path.join('/'.join(args.SAM1.split('/')[:-1]), 'HDBPA_output.csv'))
    
    output_plot = os.path.join('/'.join(args.SAM1.split('/')[:-1]), 'HDBPA_plot.png')
    plot(fc_2(df_phylo, pop_size_1, pop_size_2), output_plot)
    
def get_fasta_tag(file):
    with open(file, 'r') as f:
        for line in f:
            if '>' in line:
                seq_tag = line.rstrip().rsplit('>')[1]
                    
    return seq_tag
    
def reform_aligned_reads(input_samfile, ref, mode):
    df_reforemed = sequence_in_backbone(CCIGAR_data(input_samfile), ref)

    # store reformed aligned sequences
    output_fasta = input_samfile.split('.sam')[0] + '_reformed.fa'
    write_file_bulk(df_reforemed, output_fasta)
    if mode == 'AA':
        output_csv = input_samfile.split('.sam')[0] + '_reformed.csv'
        mutation_tracker_csv(output_fasta, output_csv, ref)
    return len(df_reforemed)

def CCIGAR_data(input_sam_file):
    
    QNAME_list, POS_list,MAPQ_list,  CIGAR_list, SEQ_list = [], [], [], [], []
    with open(input_sam_file, 'r') as f:
        for line in f:
            if '@' not in line:
                QNAME_list.append(line.rstrip().rsplit('\t')[0])
                POS_list.append(line.rstrip().rsplit('\t')[3])
                MAPQ_list.append(line.rstrip().rsplit('\t')[4])
                CIGAR_list.append(line.rstrip().rsplit('\t')[5])
                SEQ_list.append(line.rstrip().rsplit('\t')[9])
    
    df = pd.DataFrame(
        {'QNAME': QNAME_list, 
         'POS': POS_list, 
         'MAPQ': MAPQ_list,
         'CIGAR': CIGAR_list, 
         'SEQ': SEQ_list}
    )
    
    # remove duplicate sequence, keep the one with a higher MAPping Quality
    df.sort_values('MAPQ', ascending=False).drop_duplicates('QNAME').sort_index()

    return df

def seqlen_CIGAR(CIGAR):
    seqlen = 0
    alignment_match = re.findall(r'(?P<M_count>[\d]+)M', CIGAR)
    ref_del = re.findall(r'(?P<D_count>[\d]+)D', CIGAR)
    
    for item in alignment_match:
        seqlen += int(item)
    for item in ref_del:
        seqlen += int(item)
    
    return seqlen

def seq_CIGAR(CIGAR, POS, ref_sequence, SEQ_old):

    # get the reference sequence
    sequences_object = pysam.FastaFile(ref_sequence)
    # get the seqID
    seq_tag = get_fasta_tag(ref_sequence)
    
    # read the CIGAR string
    # store the int and characters separately in a data frame
    match_x = re.findall(r'[\d]+', CIGAR)
    match_y = re.findall(r'[\D]+', CIGAR)
    
    data = {'x': match_x, 'y': match_y}
    df = pd.DataFrame(data)

    start_pos_ref = POS - 1
    start_pos_seq = 0
    seq = ''
    
    # built the query sequence into the reference sequence
    # remove insertions in the query sequence
    # fill the deletions in the query sequence using the corresponding sequence in the ref
    for i in range(len(df)):
        if df.iloc[i]['y'] == 'M': # alignment match
            seq += SEQ_old[start_pos_seq: start_pos_seq + int(df.iloc[i]['x'])]
            start_pos_seq += int(df.iloc[i]['x'])
            start_pos_ref += int(df.iloc[i]['x'])
        elif df.iloc[i]['y'] == 'D': # del to the ref
            seq += sequences_object.fetch(seq_tag, start_pos_ref, 
                                          start_pos_ref + int(df.iloc[i]['x']))
            start_pos_ref += int(df.iloc[i]['x'])
        elif df.iloc[i]['y'] == 'I': # ins to the ref
            start_pos_seq += int(df.iloc[i]['x'])
        elif df.iloc[i]['y'] == 'S': # soft clipping in the query sequence
            start_pos_seq += int(df.iloc[i]['x'])
    
    return seq

def sequence_in_backbone(df, ref_sequence):  
    
    # get the reference sequence and its seqID and length
    sequences_object = pysam.FastaFile(ref_sequence)
    seq_tag = get_fasta_tag(ref_sequence)
    ref_len = len(sequences_object.fetch(seq_tag))
    
    seq_list = []
    seq_len_list = []
    
    for i in range(len(df)):
        seq = ''
        POS = int(df.iloc[i]['POS']) - 1
        CIGAR = df.iloc[i]['CIGAR']
        SEQ_old = df.iloc[i]['SEQ']
        
        seqlen = seqlen_CIGAR(CIGAR)
        flanking_1 = sequences_object.fetch(seq_tag, 0, POS)
        flanking_2 = sequences_object.fetch(seq_tag, POS + seqlen, ref_len)
        
        seq = flanking_1 + seq_CIGAR(CIGAR, POS, ref_sequence, SEQ_old) + flanking_2
        seq_list.append(seq)
        seq_len_list.append(len(seq))
        if len(seq) != ref_len:
            print('ERROR!')
            sys.exit()
        
    df['Output_seq'] = seq_list
    df['Output_seq_length'] = seq_len_list
    
    return df

def write_file_bulk(df, output_file_path):
    with open(output_file_path, 'w') as fo:
        for i in range(len(df)):
            fo.write(EOL.join([
                f">{df.iloc[i]['QNAME']}",
                f"{df.iloc[i]['Output_seq']}{EOL}"])
                    )

def mutations_identifier(input_file_path, ref):
    # the data should be stored as ReadID-(Mutations)
    
    line_count = 0
    seq = ''
    
    ReadID_list = []
    mutations = []
    mutation_list = []
    
    with open(input_file_path, 'r') as f:
        for line in f:
            line_count += 1
            if line_count % 2 == 0:
                seq = line.rstrip()
                ReadID_list.append(line_count//2)
                mutations = translator(seq, ref) # only non-synonymous mutations here
                mutation_list.append(mutations)
                    
    data = {'Read_id':ReadID_list, 'Mutations':mutation_list}
    df = pd.DataFrame(data)
    return df

def mutation_tracker_csv(input_file_path, output_file_path, ref):

    df = mutations_identifier(input_file_path, ref)
    df.to_csv(output_file_path)
    
def codon_lib_generator():
    Base1_list = ['T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T',
                  'T', 'T', 'T', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                  'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                  'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G',
                  'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G']
    Base2_list = ['T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A',
                  'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C',
                  'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T',
                  'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G',
                  'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A',
                  'G', 'G', 'G', 'G']
    Base3_list = ['T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G']
    AA_list = ['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C',
               'C', '*', 'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H',
               'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T',
               'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 'V', 'V', 'V', 'V',
               'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G']
    codon_lib = {}
    
    for i in range(len(Base1_list)):
        codon = Base1_list[i] + Base2_list[i] + Base3_list[i]
        codon_lib[codon] = AA_list[i]
    return codon_lib
    
def ORF_start_definer_HIV(position):
    # when p6 overlaps with PR, here not consider p6.
    # assume p6 ends at where PR starts
    
    start_position = 0
    region = ''
    
    if position < 2253:
        # to the end of p1
        if position <= 1185:
            start_position = 790
            region = 'MA'
        elif position <= 1878:
            start_position = 1186
            region = 'CA'
        elif position <= 1920:
            start_position = 1879
            region = 'p2'
        elif position <= 2085:
            start_position = 1921
            region = 'NC'
        elif position <= 2133:
            start_position = 2086  
            region = 'p1'
        else:
            start_position = 2134
            region = 'p6'
            
    elif position > 2253:
        if position <= 2549:
            start_position = 2253
            region = 'PR'
        elif position <= 3869:
            start_position = 2550
            region = 'RT'
        elif position <= 4229:
            start_position = 3870
            region = 'RNaseH'
        elif position <= 5093:
            start_position = 4230
            region = 'IN'
            
    return start_position, region

def translator(nuc_seq, ref_sequence):
    # customized for HIV HXB2 sequence.
    # only capture non-synonymous amino acid mutations here.
    
    codon_lib = codon_lib_generator()
    
    # get ref sequence and its sedID
    sequences_object = pysam.FastaFile(ref_sequence)
    seq_tag = get_fasta_tag(ref_sequence)
    aa_mutation_list = []
    nuc_mutation_list = []
    start_position, i = 0, 0
    
    for i in range(789, 2252, 3):
        codon_ref = ''
        codon_seq = ''
        aa_ref = ''
        aa_seq = ''
        aa_mutation = ''
        nuc_mutation = ''
        pos_list = [i+1, i+2, i+3]
        
        start_position, region = ORF_start_definer_HIV(pos_list[0])
        nuc_pos = (pos_list[0]-start_position)//3 + 1
        
        codon_ref = sequences_object.fetch(seq_tag, i, i+3)
        codon_seq = nuc_seq[i:i+3]
        aa_ref = codon_lib[codon_ref]
        aa_seq = codon_lib[codon_seq]
        
        if aa_seq != aa_ref:
            aa_mutation = region + '_' + aa_ref + str(nuc_pos) + aa_seq
            
            for i, j in zip(list(codon_ref), find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += str(i)
            for i, j in zip(pos_list, find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += '_'
                    nuc_mutation += str(i)
                    nuc_mutation += '_'
            for i, j in zip(list(codon_seq), find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += str(i)
            
            aa_mutation_list.append(aa_mutation)
            nuc_mutation_list.append(nuc_mutation)
           
    for i in range(2252, 5093, 3):
        codon_ref = ''
        codon_seq = ''
        aa_ref = ''
        aa_seq = ''
        aa_mutation = ''
        nuc_mutation = ''
        pos_list = [i+1, i+2, i+3]
        
        start_position, region = ORF_start_definer_HIV(pos_list[0])
        nuc_pos = (pos_list[0]-start_position)//3 + 1
        
        codon_ref = sequences_object.fetch(seq_tag, i, i+3)
        codon_seq = nuc_seq[i:i+3]
        aa_ref = codon_lib[codon_ref]
        aa_seq = codon_lib[codon_seq]
            
        if aa_seq != aa_ref:
            aa_mutation = region + '_' + aa_ref + str(nuc_pos) + aa_seq
            
            for i, j in zip(list(codon_ref), find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += str(i)
            for i, j in zip(pos_list, find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += '_'
                    nuc_mutation += str(i)
                    nuc_mutation += '_'
            for i, j in zip(list(codon_seq), find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += str(i)
            
            aa_mutation_list.append(aa_mutation)
            nuc_mutation_list.append(nuc_mutation)
        
    return aa_mutation_list
    
def find_different_nuc(codon_ref, codon_seq):
    
    compare_list = []
    
    for i, j in zip(codon_ref, codon_seq):
        compare_list.append(i==j)
        
    return compare_list
    
def phylogeny_all(sample1, sample2, threshold = 0):
    
    parent_list = []
    progeny_list = []
    distance_list = []
    dif_list = []
    Afre_list = []
    Pfre_list = []
    
    Patterns1 = pattern_abundance_calculator(sample1)['Pattern'].tolist()
    Abundance1 = pattern_abundance_calculator(sample1)['Abundance'].tolist()
    
    Patterns2 = pattern_abundance_calculator(sample2)['Pattern'].tolist()
    Abundance2 = pattern_abundance_calculator(sample2)['Abundance'].tolist()

    for P2 in tqdm(Patterns2):
        distance = len(P2[2:-2].split("', '"))
        for P1 in Patterns1:
            if get_distance(P1[2:-2].split("', '"), P2[2:-2].split("', '"))[0] < distance:
                distance, dif = get_distance(P1[2:-2].split("', '"), P2[2:-2].split("', '"))
                ancestor = P1
                
            else:
                pass

        parent_list.append(Patterns1.index(ancestor) + 1)
        progeny_list.append(Patterns2.index(P2) + 1)
        distance_list.append(distance)
        dif_list.append(dif)
        Afre_list.append(Abundance1[Patterns1.index(ancestor)])
        Pfre_list.append(Abundance2[Patterns2.index(P2)])
        
    data = {
        'Parent_strain':parent_list,
        'Progeney_strain': progeny_list,
        'Distance': distance_list,
        'Differences': dif_list,
        'Parent_fre': Afre_list,
        'Progeny_fre': Pfre_list
    }
    df = pd.DataFrame(data)
    return df
    
def pattern_abundance_calculator(input_file_path):
    
    df = pd.read_csv(input_file_path)

    mutations_list = df['Mutations'].tolist()
    abundance_lib = {}

    for item in mutations_list:
        if item in abundance_lib:
            abundance_lib[item] += 1
        else:
            abundance_lib[item] = 1

    pattern_list = list(abundance_lib.keys())
    abundance_list = list(abundance_lib.values())

    data_pattern = {'Pattern': pattern_list,
                 'Abundance':abundance_list}

    df_pattern = pd.DataFrame(data_pattern)
    df_pattern = df_pattern.sort_values(by = 'Abundance', ascending=False)
    
    return df_pattern

def get_distance(list_a, list_b):
    dif_add = list(set(list_a) - set(list_b))
    for i in range(len(dif_add)):
        dif_add[i] = '+ ' + dif_add[i]
    dif_minus = list(set(list_b) - set(list_a))
    for i in range(len(dif_minus)):
        dif_minus[i] = '- ' + dif_minus[i]
    return len(dif_add + dif_minus), dif_add + dif_minus
    
def plot(df, output_plot):
    blues = cm.get_cmap("Blues", 10)
    df['Parent_fre'] = df['Parent_fre']/sum( df['Parent_fre'].tolist())
    df['Progeny_fre'] = df['Progeny_fre']/sum( df['Progeny_fre'].tolist())
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111, aspect='equal')

    x = [0, 0, 1, 1]
    y1, y2, y3, y4 = 0, 0, 0, 0
    for i in range(len(df)):
        y2 += df.iloc[i]['Parent_fre']
        y3 += df.iloc[i]['Progeny_fre']
        y = [y1, y2, y3, y4]
        ax.add_patch(patches.Polygon(xy=list(zip(x,y)), color=blues(i)))
        # add text: MRCA ID
        plt.text(-0.2, y1 + (y2 - y1)/2,
                 f'A{int(df.iloc[i]["Parent_strain"])}',
                 verticalalignment='center'
                )
        # add a short horizontal reference line
        lines = [[(-0.01, y1 + (y2 - y1)/2,), 
                  (0, y1 + (y2 - y1)/2,)]]
        lc = collections.LineCollection(lines, color = 'black', linewidths=1)
        ax.add_collection(lc)

        y1 = y2
        y4 = y3
    
    # add time points on x-axis
    plt.text(0, -0.1, 'Pre-treatment')
    plt.text(0.75, -0.1, 'Post-treatment')
    
    # add white margins
    ax.autoscale()
    ax.margins(0.2)
    
    plt.xticks([]),plt.yticks([])
    plt.savefig(output_plot, bbox_inches='tight')
    return

def fc_2(df_phy, pop_size_1, pop_size_2):
    Parent_strain_list = []
    Afre_list = []
    Pfre_list = []
    fold_list = []
    
    for item in list(set(df_phy['Parent_strain'].tolist())):
        df = df_phy[df_phy['Parent_strain'] == item]
        
        Parent_strain_list.append(item)
        Afre_list.append(round(df['Parent_fre'].tolist()[0]/pop_size_1,4))
        Pfre_list.append(round(sum(df['Progeny_fre'].tolist())/pop_size_2,4))

    df_fc = pd.DataFrame({
        'Parent_strain': Parent_strain_list,
        'Parent_fre': Afre_list,
        'Progeny_fre':Pfre_list
    })

    df_fc['Fold_change'] = df_fc['Progeny_fre']/df_fc['Parent_fre']

    return df_fc.sort_values(by = 'Progeny_fre', ascending = False).head(10)