import os 

import click
import matplotlib.pyplot as plt
import numpy as np

def read_file(file,kmer_size):

    isHeader = False
    result = {}

    with open(file) as f:
        for line in f:
            if(line.startswith('@') or line.startswith('>')):
                isHeader = True
            elif(line.startswith('+')):
                isHeader = False
            elif(isHeader == True):
                i=0
                while( i+kmer_size < len(line)):
                    kmer = line[i:i+kmer_size]    
                    if(kmer in result):
                        result[kmer] += 1

                    else:
                        result[kmer] = 1
                    i += 1

    return result

def plot_counts(kmer_dict):
    kmer_counts = {}
    for value in kmer_dict.values():
        if(value in kmer_counts):
            kmer_counts[value] += 1
        else:
            kmer_counts[value] = 1
    
    sorted_dict = sorted(kmer_counts.items())
    x_val =[]
    y_val = []
    for val in sorted_dict:
        x_val += [val[0]]
        y_val += [val[1]]

    plt.figure(figsize=(10, 6))
    plt.bar(x_val, y_val)
    plt.xlabel('K-mers')
    plt.ylabel('Counts')
    plt.title('K-mer Counts')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('./img/kmer_distrib.png')
    plt.show()

@click.command()
@click.option('-k', '--kmer_size', default=31, help='Kmer size used for the alignment')
@click.option('-f', '--reads_file', help='Relative path to the fasta file')
def main(kmer_size, reads_file):
    kmer_dict = read_file(reads_file,kmer_size)
    plot_counts(kmer_dict)

if __name__ == '__main__':
    main()