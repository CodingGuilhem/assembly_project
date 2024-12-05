from typing import Node, Dict

import arbre_suff

def findBestWeight(currentNode: Node[arbre_suff.NoeudArbre],path: str):
    """
    Find the best weight (the greatest) in the currentNode tree (or subtree)
    """
    weight = []
    path = []
    # Find all the weight possible with their node
    for nucl in ['A','C','T','G']:
        nextNode = getattr(currentNode,nucl) 
        if nextNode == None:
            return nextNode.number, nextNode.findPath(nextNode)
        else:
            weight,path += findBestWeight(nextNode)
        
    print(max(weight))
    print(path,weight)


def trouveProchainKmer(string: str, racineArbreKmer, tailleK=31):
    kPossible = []
    for i in range(1, len(string)):
        kPossible.append(string[i:len(string)])  # On selectionne les kmer possibles dans l'ordre décroissant pour selectionner le plus grand possible
    print(kPossible)
    bestKmer = {'a':0}
    for kmer in kPossible:
        currentNode = racineArbreKmer
        kmerLength = len(kmer)
        pos = 0
        # Go to the last node of the Kmer if possible 
        while kmerLength != 1:
            nextNode = getattr(currentNode,kmerLength[pos])
            if(nextNode is not None):
                currentNode = nextNode
                kmerLength -= 1
            else:
                currentNode = None
                break
        
        #If we went to the last node possible of the Kmer we look for the best value in this subtree
        if currentNode != None:
            bestWeight = findBestWeight(currentNode)

        


def assemble(kmerDebut,racineArbreKmer,tailleK = 31):
    res = kmerDebut
    prochainKmer = trouveProchainKmer(res,racineArbreKmer,tailleK) 
    # while(prochainKmer != None):
    #     res = assembler(res,prochainKmer)
    #     assemble(kmerDebut= res,racineArbreKmer)
    
    # return res