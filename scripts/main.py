import os 

import click
import matplotlib.pyplot as plt
import numpy as np

import arbre_suff as suff


@click.command()
@click.option('-k', '--kmer_size', default=31, help='Kmer size used for the alignment')
@click.option('-f', '--reads_file', help='Relative path to the fasta file')
def main(kmer_size, reads_file):
    arbre = suff.arbre_suff(reads_file,kmer_size)
    arbre_s = suff.nettoyer_arbre(arbre,1)

        

if __name__ == '__main__':
    main()