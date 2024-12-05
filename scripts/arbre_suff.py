import copy
import sys
from collections import defaultdict, deque
import time 

sys.setrecursionlimit(20000)

class NoeudArbre:
    def __init__(self, etiquette):
        self.etiquette = etiquette
        self.A = None
        self.C = None
        self.G = None
        self.T = None
        self.parent = None 
        self.number = 0

    def chemin(self):
        chemin = ""
        noeud = self
        while noeud.parent:
            chemin = noeud.etiquette + chemin
            noeud = noeud.parent

        return chemin
    
    def est_feuille(self):
        return not (self.A or self.C or self.G or self.T)
    
    def sup(self,etiquette_noeud=''):
        noeud_courant = self


        for nucleotide in etiquette_noeud:
            if noeud_courant:
                noeud_courant = getattr(noeud_courant,nucleotide)

        if noeud_courant and noeud_courant.est_feuille() and noeud_courant.parent:
            setattr(noeud_courant.parent,noeud_courant.etiquette,None)
            noeud_courant.parent.sup()


def arbre_suff(file_path, k):
    # Initialisation de la racine de l'arbre avec une étiquette non spécifique "Racine".
    root = NoeudArbre("Racine")
    
    # Ouverture du fichier FASTQ situé au chemin spécifié par 'file_path'.
    with open(file_path, 'r') as fastq_file:
        line_count = 0  # Compteur de lignes pour naviguer dans le format FASTQ.
        
        # Parcourir chaque ligne du fichier FASTQ.
        for line in fastq_file:
            line_count += 1  # Incrémentation du compteur de lignes
            
            # Sélection des lignes contenant des séquences
            if (line_count - 2) % 4 == 0 and len(line.strip()) >= k:
                # Pour chaque position de départ possible pour un k-mer dans la ligne.
                for i in range(len(line.strip()) - (k - 1)):
                    noeud_courant = root  # Commencer à la racine pour chaque nouveau k-mer.
                    kmer = line.strip()[i:i+k]  # Extraction du k-mer à la position i.
                    
                    # Pour chaque nucléotide dans le k-mer.
                    for nuc in range(len(kmer)):
                        # Si le nucléotide existe déjà comme enfant, se déplacer vers ce noeud.
                        if getattr(noeud_courant, kmer[nuc]):
                            noeud_courant = getattr(noeud_courant, kmer[nuc])
                        else:
                            # Sinon, créer un nouveau noeud pour ce nucléotide et le lier au noeud courant.
                            setattr(noeud_courant, kmer[nuc], NoeudArbre(kmer[nuc]))
                            parent = noeud_courant
                            noeud_courant = getattr(noeud_courant, kmer[nuc])
                            setattr(noeud_courant, 'parent', parent)
                            setattr(noeud_courant, 'etiquette', kmer[nuc])
                        
                        # Si c'est le dernier nucléotide du k-mer, incrémenter le compteur de ce noeud.
                        if nuc == len(kmer) - 1:
                            setattr(noeud_courant, 'number', getattr(noeud_courant, 'number') + 1)
                            
    # Retourne la racine de l'arbre après avoir ajouté tous les k-mers du fichier.
    return root


def nettoyer_arbre(noeud: NoeudArbre, n):
    # Parcours des enfants du nœud pour chaque nucléotide possible (A, C, G, T).
    for nucleotide in ['A', 'C', 'G', 'T']:
        # Récupération de l'enfant correspondant au nucléotide courant.
        enfant = getattr(noeud, nucleotide)

        # Si l'enfant existe, procéder à un nettoyage récursif de cet enfant.
        if enfant is not None:
            nettoyer_arbre(enfant, n)
            # Après le retour de la récursion, vérifier si l'enfant est une feuille et si son compteur est inférieur ou égal à n.
            if enfant.est_feuille() and enfant.number <= n:
                # Si les conditions sont remplies, supprimer cet enfant en mettant son lien à None.
                setattr(noeud, nucleotide, None)

    # La fonction retourne le nœud actuel après avoir potentiellement modifié ses liens enfants.
    return noeud


def compter_feuilles(noeud):
    if getattr(noeud, 'A') is None and getattr(noeud, 'C') is None and \
       getattr(noeud, 'G') is None and getattr(noeud, 'T') is None: return 1 
    
    count = 0
    for nucleotide in ['A', 'C', 'G', 'T']:
        enfant = getattr(noeud, nucleotide)
        if enfant is not None: count += compter_feuilles(enfant)
    
    return count


def make_dic(arbre_suffixes: NoeudArbre, arbre_suffixes_suivie: NoeudArbre, noeud_courant, dic_de_bruijn, dic_suivie, n):
    # Obtenir le suffixe du noeud courant pour suivre le chemin suivant
    suffixe = noeud_courant[1:]
    
    # Si le noeud courant n'est pas déjà dans le dictionnaire à cet indice, l'ajouter
    if noeud_courant not in dic_de_bruijn[n]:
        dic_de_bruijn[n][noeud_courant] = []
    
    # Associer le suffixe au noeud courant dans le dictionnaire de suivi
    dic_suivie[suffixe] = [noeud_courant, n]
    
    # Commencer à suivre le chemin à partir de l'arbre des suffixes
    chemin_next = arbre_suffixes
    
    # Supprimer le noeud courant de l'arbre suivie pour éviter la répétition
    arbre_suffixes_suivie.sup(noeud_courant)
    
    # Suivre le chemin à partir du suffixe pour trouver le prochain noeud
    for nuc in suffixe:
        chemin_next = getattr(chemin_next, nuc)
        # Si le chemin est interrompu, retourner les dictionnaires mis à jour
        if chemin_next is None:
            return arbre_suffixes_suivie, dic_de_bruijn, dic_suivie

    # Explorer toutes les transitions possibles à partir du noeud courant
    for nuc in ['A', 'C', 'G', 'T']:
        if getattr(chemin_next, nuc):
            # S'assurer que l'index n est correct et mis à jour
            while isinstance(dic_de_bruijn[n], int):
                n = dic_de_bruijn[n]
            if noeud_courant not in dic_de_bruijn[n]:
                dic_de_bruijn[n][noeud_courant] = []
            # Obtenir le prochain noeud et son chemin
            noeud_next = getattr(chemin_next, nuc)
            noeud_next_chemin = noeud_next.chemin()
            dic_de_bruijn[n][noeud_courant].append(noeud_next_chemin)
            
            # Ajouter des connexions ou mettre à jour des chemins dans le graphe de De Bruijn
            if noeud_next_chemin[1:] not in dic_suivie:
                # Si le suffixe du chemin suivant (noeud_next_chemin sans le premier caractère)
                # n'est pas déjà dans le dictionnaire dic_suivie, cela signifie que nous n'avons pas encore
                # traité ce chemin. Nous devons donc le traiter en appelant récursivement make_dic.
                make_dic(arbre_suffixes, arbre_suffixes_suivie, noeud_next_chemin, dic_de_bruijn, dic_suivie, n)
            else:
                # Si le suffixe est déjà dans dic_suivie, cela signifie que nous avons déjà traité
                # ce chemin et nous devons gérer les duplications potentielles ou les connexions redondantes.
                n_prime = dic_suivie[noeud_next_chemin[1:]][1]  # Obtenir l'identifiant associé au chemin déjà traité.
                if n == n_prime:
                    # Si l'identifiant du noeud actuel est le même que celui trouvé dans dic_suivie,
                    # il n'y a pas de nouvelle branche, donc nous copions les connexions existantes pour ce chemin
                    # dans le graphe de De Bruijn et nous supprimons ce noeud de l'arbre_suffixes_suivie pour éviter des traitements redondants.
                    dic_de_bruijn[n][noeud_next_chemin] = dic_de_bruijn[n_prime][dic_suivie[noeud_next_chemin[1:]][0]]
                    arbre_suffixes_suivie.sup(noeud_next_chemin)
                else:
                    # Si les identifiants diffèrent, cela indique que les chemins se sont croisés à un point précédent
                    # et que nous devons fusionner les informations de graphe sous un même identifiant pour assurer
                    # la cohérence du graphe. Toutes les connexions et les références du noeud actuel sont transférées
                    # sous l'identifiant de l'autre noeud (n_prime).
                    dic_de_bruijn[n_prime].update(dic_de_bruijn[n])
                    for i in dic_de_bruijn[n].keys():
                        # Mettre à jour dic_suivie pour refléter le changement d'identifiant pour tous les noeuds concernés.
                        dic_suivie[i[1:]][1] = n_prime
                    dic_de_bruijn[n] = n_prime  # Mettre à jour l'identifiant actuel pour correspondre à n_prime
                    n = n_prime  # Utiliser n_prime comme le nouvel identifiant pour les opérations futures.


    return arbre_suffixes_suivie, dic_de_bruijn, dic_suivie



def graph_de_bruijn(arbre_suffixes: NoeudArbre, arbre_suffixes_suivie: NoeudArbre, dic_de_bruijn={}, dic_suivie={}, nb=1, noeud=None):
    # Continuer tant que l'arbre de suffixes suivie n'est pas entièrement réduit à une feuille
    while not arbre_suffixes_suivie.est_feuille():
        # Initialiser le dictionnaire pour ce nombre si ce n'est pas déjà fait
        if nb not in dic_de_bruijn:
            dic_de_bruijn[nb] = {}
            noeud = arbre_suffixes_suivie
        
        # Trouver le premier noeud non-feuille
        while not noeud.est_feuille():
            for nuc in ['A', 'C', 'G', 'T']:
                if getattr(noeud, nuc):
                    noeud = getattr(noeud, nuc)
                    break
        
        # Obtenir le chemin complet du noeud courant
        noeud_courant = noeud.chemin()
        dic_de_bruijn[nb][noeud_courant] = []
        arbre_suffixes_suivie.sup(noeud_courant)
        # Appliquer la fonction make_dic pour mettre à jour le graphe
        arbre_suffixes_suivie, dic_de_bruijn, dic_suivie = make_dic(arbre_suffixes, arbre_suffixes_suivie, noeud_courant, dic_de_bruijn, dic_suivie, nb)
        nb += 1
        
        # Appel récursif pour continuer la construction du graphe
        graph_de_bruijn(arbre_suffixes, arbre_suffixes_suivie, dic_de_bruijn, dic_suivie, nb)

    # Libérer l'arbre de suffixes suivie à la fin
    del arbre_suffixes_suivie
    return dic_de_bruijn




def find_eulerian_path(graph,start_node):
    stack = []
    chemin = deque()
    noeud_courant = start_node or next(iter(graph))  # Choisir un point de départ
    stack.append(noeud_courant)

    while stack:
        # Si le nœud actuel a des voisins restants
        if graph[noeud_courant]:
            stack.append(noeud_courant)  # Empiler le nœud actuel
            # Sélectionner le prochain voisin (une arête sortante) et le supprimer de la liste
            noeud_next = graph[noeud_courant].pop()
            noeud_courant = noeud_next  # Se déplacer vers le nœud suivant
        else:
            # Si aucun voisin restant, ajouter ce nœud au chemin et revenir en arrière
            chemin.appendleft(noeud_courant)
            noeud_courant = stack.pop()  # Retour à un nœud précédent

    return list(chemin)


def make_fasta(dic_de_bruijn, file_path):
    with open(file_path, "w") as file: pass
    
    for i in dic_de_bruijn.keys():
        graph =dic_de_bruijn[i]

        # Calculer les degrés entrants et sortants
        degree_sortant = defaultdict(int)
        degree_entrant = defaultdict(int)
        for node, neighbors in graph.items():
            degree_sortant[node] += len(neighbors)
            for neighbor in neighbors:
                degree_entrant[neighbor] += 1

        # Trouver les points de départ et d'arrivée
        start_node = None
        end_node = None
        for node in set(degree_sortant.keys()).union(set(degree_entrant.keys())):
            out = degree_sortant[node]
            inp = degree_entrant[node]
            if out > inp:
                start_node = node
            elif inp > out:
                end_node = node

        # Ajouter une arête artificielle pour gérer les cas où un chemin est nécessaire
        if start_node and end_node: graph[end_node].append(start_node)

        # Trouver le chemin eulérien avec l'algorithme de Hierholzer
        eulerian_path = find_eulerian_path(graph,start_node)

        # Supprimer l'arête artificielle si ajoutée
        if start_node and end_node:
            eulerian_path = eulerian_path[:-1]

        # Création du contig a partir du chemin eulérien
        start = eulerian_path[0]
        for j in eulerian_path[1:]: start+=j[-1]

        # Mise a jour du fasta de sortie
        lines_to_add = [f">{i-1}_len{len(start)}\n",f"{start}\n"]
        with open(file_path, "a") as file:  # 'a' pour append
            file.writelines(lines_to_add)


def main(file_path,k,threshold_representation,new_file_path):
    start_time = time.time()
    print("Création de l'arbre des suffixes ...")
    arbre_suffixes = arbre_suff(file_path,k) 
    arbre_suffixes = nettoyer_arbre(arbre_suffixes,threshold_representation)

    print('Création du graph de Bruijn ...')
    arbre_suffixes_suivie = copy.deepcopy(arbre_suffixes)
    dic_de_bruijn = graph_de_bruijn(arbre_suffixes,arbre_suffixes_suivie)
    dic_de_bruijn = {key: value for key, value in dic_de_bruijn.items() if not isinstance(value , int)}

    print('Assemblage ...')
    make_fasta(dic_de_bruijn,new_file_path)
    end_time = time.time()
    print(f"Temps pour l'assemblage : {end_time - start_time:.2f} secondes")
 


file_path = 'data/reads.fastq.fq'
k=17
threshold_representation = 1
new_file_path = 'data/V2.fa'
main(file_path,k,threshold_representation,new_file_path)
