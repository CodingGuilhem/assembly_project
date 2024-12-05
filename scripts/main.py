import copy
import time
import click
import arbre_suff as suff

@click.command()
@click.option('-k', '--kmer_size', default=31, help='Kmer size used for the alignment')
@click.option('-f', '--reads_file', required=True, help='Path to the fastq file')
@click.option('-t', '--threshold_representation', required=True, type=int, help='Threshold to delete the underrepresented K-mer (if k < threshold we delete the K-mer)')
@click.option('-o', '--output_dir', required=True, help='Output file name (with  relative path if needed)')
def main(kmer_size, reads_file, threshold_representation, output_dir):
    start_time = time.time()
    print("Création de l'arbre des suffixes ...")
    arbre_suffixes = suff.arbre_suff(reads_file, kmer_size)
    arbre_suffixes = suff.nettoyer_arbre(arbre_suffixes, threshold_representation)

    print('Création du graph de Bruijn ...')
    arbre_suffixes_suivie = copy.deepcopy(arbre_suffixes)
    dic_de_bruijn = suff.graph_de_bruijn(arbre_suffixes, arbre_suffixes_suivie)
    dic_de_bruijn = {key: value for key, value in dic_de_bruijn.items() if not isinstance(value, int)}

    print('Assemblage ...')
    suff.make_fasta(dic_de_bruijn, output_dir)
    end_time = time.time()
    print(f"Temps pour l'assemblage : {end_time - start_time:.2f} secondes")

if __name__ == '__main__':
    main()