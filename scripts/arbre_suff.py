import matplotlib.pyplot as plt
import numpy as np



def kmer_dict (file,seuil,k):
    file_path = file
    kmer_dict = {}
    with open(file_path, 'r') as fastq_file:
        line_count = 0
        for line in fastq_file:
            line_count += 1
            if (line_count - 2) % 4 == 0:
                if len(line.strip()) >= k:
                    for i in range(len(line.strip())-(k-1)):
                        kmer = line.strip()[i:i+k]
                        if kmer not in kmer_dict : kmer_dict[kmer] = 1
                        else : kmer_dict[kmer] += 1

    return {k: v for k, v in kmer_dict.items() if v > seuil}


def plot_kmerdict(kmer_dict):
    freq_dict = {}
    for freq in kmer_dict.values():
        if freq not in freq_dict:
            freq_dict[freq] = 1
        else:
            freq_dict[freq] += 1

    # Trier les fréquences par ordre croissant
    sorted_freqs = sorted(freq_dict.items())
    print(sorted_freqs)
    # Préparer les données pour le plot
    frequencies = [x[0] for x in sorted_freqs]
    counts = [x[1] for x in sorted_freqs]

    # Créer le bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(frequencies, counts)
    plt.title('Distribution des k-mers par fréquence')
    plt.xlabel('Occurence des k-mers')
    plt.ylabel('Nombre de k-mers')
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.show()



class NoeudArbre:
    def __init__(self, etiquette):
        self.etiquette = etiquette
        self.A = None
        self.C = None
        self.G = None
        self.T = None
        self.parent = None 
        self.number = 0

    def printNoeud(self):
        print(f"etiquette : {self.etiquette}")
        if self.A != None:
            print(f"nucleotide A")
        if self.C != None:
            print(f"nucleotide C")
        if self.T != None:
            print(f"nucleotide T")
        if self.G != None:
            print(f"nucleotide G")
        if self.estFeuille():
            print("feuille")
        if self.parent != None:
            print(f"parent : {self.parent.etiquette}")
        else:
            print("racine")
        print(f"numero : {self.number} ")

    def estFeuille(self):
        if(self.A == None and self.C == None and self.T == None and self.G == None):
            return True 
        else: return False


def findPath(node):
    parent = node.Parent
    path = node.etiquette 
    if parent == None:
        return parent.etiquette
    else:
        path += findPath(parent)

    return path

def printArbre(noeud):
    if noeud.estFeuille():
        noeud.printNoeud()        
    else:
        if(noeud.A != None):
            printArbre(noeud.A)
        if(noeud.T!= None):
            printArbre(noeud.T)
        if(noeud.C != None):
            printArbre(noeud.C)
        if(noeud.G != None):
            printArbre(noeud.G)
        
        noeud.printNoeud()

def arbre_suff(file_path,k):
    root = NoeudArbre("Racine")
    with open(file_path, 'r') as fastq_file:
        line_count = 0
        for line in fastq_file:
            line_count += 1
            if (line_count - 2) % 4 == 0 and len(line.strip()) >= k:
                for i in range(len(line.strip())-(k-1)):
                    noeud_courand = root
                    kmer = line.strip()[i:i+k]
                    for nuc in range(len(kmer)):
                        if getattr(noeud_courand,kmer[nuc]): noeud_courand = getattr(noeud_courand,kmer[nuc])
                        else : 
                            setattr(noeud_courand,kmer[nuc],NoeudArbre(kmer[nuc]))
                            parent = noeud_courand
                            noeud_courand = getattr(noeud_courand,kmer[nuc])
                            setattr(noeud_courand,'parent', parent)
                            setattr(noeud_courand,'etiquette', kmer[nuc])
                        if nuc == len(kmer)-1 : 
                            setattr(noeud_courand,'number', getattr(noeud_courand,'number')+1) 
                            
    return root          

def nettoyer_arbre(noeud, n):
    # Parcours récursif des enfants A, C, G, T
    for nucleotide in ['A', 'C', 'G', 'T']:
        enfant = getattr(noeud, nucleotide)
        if enfant is not None:
            nettoyer_arbre(enfant,n)
            if getattr(enfant, 'A') is None and getattr(enfant, 'C') is None and \
               getattr(enfant, 'G') is None and getattr(enfant, 'T') is None and \
               enfant.number <= n: setattr(noeud, nucleotide, None)
    return noeud
            
def compter_feuilles(noeud):
    if getattr(noeud, 'A') is None and getattr(noeud, 'C') is None and \
       getattr(noeud, 'G') is None and getattr(noeud, 'T') is None: return 1 
    
    count = 0
    for nucleotide in ['A', 'C', 'G', 'T']:
        enfant = getattr(noeud, nucleotide)
        if enfant is not None: count += compter_feuilles(enfant)
    return count



