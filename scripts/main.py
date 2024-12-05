import os 

import click
import matplotlib.pyplot as plt
import numpy as np

import arbre_suff as suff
import assemblage

@click.command()
@click.option('-k', '--kmer_size', default=31, help='Kmer size used for the alignment')
@click.option('-f', '--reads_file', help='Relative path to the fasta file')
def main(kmer_size, reads_file):
    testfile = '/home/guilhem/cours/Bioinformatique avancée/Cours Annie/TP Assemblage-20241114/assembly_project/data/test.fq'
    test_size = 3
    arbre = suff.arbre_suff(testfile,test_size)
    suff.printArbre(arbre.A)
    assemblage.assemble('AAA', arbre, test_size)
    # arbre_s = suff.nettoyer_arbre(arbre,1)
    # suff.printArbre(arbre_s)
        

if __name__ == '__main__':
    main()