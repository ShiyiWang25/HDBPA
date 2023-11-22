import argparse


def _parse_args(args, **kwargs):

    parser = argparse.ArgumentParser(
                        prog='Hamming-distance based phylogenetic analysis',
                        description=('Link sequences to their MRCAs detected'
                                     'at a previous time point.'
                                     'Reveal genetic changes in individual genomes'
                                     'as well as in subpopulations.'),
                        epilog='')
    
    parser.add_argument('--SAM1',
                        type=str,
                        help='Input aligned reads in the 1st sample.')  
    parser.add_argument('--SAM2',
                        type=str,
                        help='Input aligned reads in the 2nd sample.')  
    parser.add_argument('-r',
                        type=str,
                        help='Import the reference sequence.')     
    parser.add_argument('--mode',
                        type=str,
                        default = 'AA',
                        choices=['NT', 'AA'],
                        help='Calculate the Hamming distance on the nucleotide level or'
                        'the non-synonymous amino acid level.')
    return parser.parse_args(args)