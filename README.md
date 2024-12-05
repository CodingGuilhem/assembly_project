# Projet assemblage  
## Introduction

Ce projet consiste à créer un assembleur de génome de nos propre main
## Prérequis

- Python 3.10.12
- Bibliothèques Python (click, copy)

## Installation

1. Clonez le dépôt :
    ```bash
    git clone https://github.com/votre-utilisateur/assembly_project.git
    ```
2. Installez les dépendances :
    ```bash
    pip install -r pyrequirements.txt
    ```

## Utilisation

1. Préparez vos données de séquençage en format FASTQ.
2. Exécutez le script d'assemblage :
    ```bash
    python script/main.py -k kmer_size -f fastQ_file -t threshold -o output_directory
    ```

## Résultats

Les résultats de l'assemblage seront disponibles dans le répertoire spécifié. Vous y trouverez les fichiers suivants :
- `Assamblage.fasta` : les contigs assemblés


## Auteurs

- Hugo Bellavoir
- Amel Benarbia
- Guilhem Biosse

## Licence

Ce projet est sous licence MIT. Voir le fichier [LICENSE](LICENSE) pour plus de détails.