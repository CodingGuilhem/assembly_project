# on parcourt le fichier fasta
# on ne lit que les lignes avec les sequences des reads
# au fur et a mesure de la lecture des sequences, on recuperes les k-mers que l on place dans des dictionnaires

read1 = "ATGCTATGCT"

kmer_len = 3
kmer = {}

if len(read1) < kmer_len:
    print("Error: the length of the sequence is less than the size of the k-mer.")
else:
    i=0
    while i+kmer_len<=len(read1): # verifier que le read est au moins aussi long que la longueur de k-mer choisie
        substring = read1[i:i+kmer_len]
        # ajouter ou incrementer l occurrence du k-mer
        if substring in kmer:
            kmer[substring] += 1
        else:
            kmer[substring] = 1
        i += 1

print(kmer)

