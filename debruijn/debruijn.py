#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import pytest
#import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "BEN HAMZA Rayed"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["BEN HAMZA Rayed"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "BEN HAMZA Rayed"
__email__ = "rayed1997@outlook.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):
#prend un seul argument correspondant au fichier fastq et retourne un générateur de séquences
    with open(fastq_file, "r") as file:
        for line_number, n in enumerate(file):
            if line_number % 4 == 1:
                yield n.strip()



def cut_kmer(read, kmer_size):
#retourne un générateur de k-mer   
    decal = 0
    seq_length = len(read)
    while decal + kmer_size <= seq_length:
        yield read[decal:(decal + kmer_size)]
        decal += 1

def build_kmer_dict(fastq_file, kmer_size):
#retourne un dictionnaire ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
    reads = read_fastq(fastq_file)
    kmer_count = {}
    for read in reads:
        kmers = cut_kmer(read, kmer_size)
        for kmer in kmers:
            if not kmer in kmer_count:
                kmer_count[kmer] = 1
            else:
                kmer_count[kmer] += 1
    return kmer_count


def build_graph(kmer_dict):
#Créer un arbre orienté  et pondéré représentant tous les k-mers préfixes et suffixes existant. 
    graph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        p = kmer[:-1]
        s = kmer[1:]
        graph.add_edge(p, s, weight=kmer_dict[kmer])
    return graph


def get_starting_nodes(graph):
#retourne une liste de noeuds d’entrée
    entry_node_list = []
    for node in graph.nodes():
        if not list(graph.predecessors(node)):
            entry_node_list.append(node)
    return entry_node_list

def get_sink_nodes(graph):
#retourne une liste de noeuds de sortie 
    exit_node_list = []
    for node in graph.nodes():
        if not list(graph.successors(node)):
            exit_node_list.append(node)
    return exit_node_list


def get_contigs(graph, starting_nodes, ending_nodes):
#Retourne une liste de tuple(contig, taille du contig)
    contigs = []
    for node_start in starting_nodes:
        for node_end in ending_nodes:
                path = nx.shortest_path(graph, node_start, node_end)
                contig = "".join([node[0] for node in path[:-1]] + [path[-1]])
                contigs.append((contig, len(contig)))
    return contigs

def fill(text, width=80):
#Retour chariot tous les 80 caractères
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
#Prend une liste de tuple et un nom de fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta
    entt = ">contig_Numéro{} len={}\n"
    with open(output_file, "w") as fichiersortie:
        for i, contig in enumerate(contigs_list):
            fichiersortie.write(entt.format(i, contig[1]))
            fichiersortie.write(fill(contig[0]))

def std(data):
    #Retourne l'écart type
    return statistics.stdev(data)

def path_average_weight(graph, path):
    #Retourne un poids moyen.
    total = 0
    for node_1, node_2 in zip(path[:-1], path[1:]):
            total += graph[node_1][node_2]["weight"]
    return total/(len(path)-1) if total else 0

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    #Retourne un graphe nettoyé des chemins indésirables
    for path in path_list:
        graph.remove_nodes_from(path[(not delete_entry_node):(None if delete_sink_node else -1)])
    return graph



def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    best_weight_indexes = [i for i, weight in enumerate(weight_avg_list)
                           if weight == max(weight_avg_list)]
    best_length_and_weights = [length for i, length in enumerate(path_length)
                               if i in best_weight_indexes]
    best_path_indexes = [i for i in best_weight_indexes
                         if path_length[i] == max(best_length_and_weights)]
    best_path_index = random.choice(best_path_indexes)
    graph = remove_paths(graph, path_list[:best_path_index]+path_list[(best_path_index+1):],
                         delete_entry_node, delete_sink_node)
    return graph



def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_lengths = [len(path_list) for path in path_list]
    avg_path_weights = [path_average_weight(graph, path) for path in path_list]

    return select_best_path(graph, path_list, path_lengths, avg_path_weights)

def simplify_bubbles(graph):
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    for ancestor_node in starting_nodes:
        for descendant_node in sink_nodes:
            current_entry = ancestor_node
            current_exit = descendant_node
            successors = list(graph.successors(current_entry))
            predecessors = list(graph.predecessors(current_exit))
            while len(successors) < 2 and successors:
                current_entry = successors[0]
                successors = list(graph.successors(current_entry))
            while len(predecessors) < 2 and predecessors:
                current_exit = predecessors[0]
                predecessors = list(graph.predecessors(current_exit))
            if list(nx.all_simple_paths(graph, current_entry, current_exit)):
                graph = solve_bubble(graph, current_entry, current_exit)
    return graph

def solve_entry_tips(graph, starting_nodes):
    graph = simplify_bubbles(graph)
    tips = []
    for start_node in starting_nodes:
        current_node = start_node
        path = [current_node]
        successors = list(graph.successors(current_node))
        predecessors = list(graph.predecessors(current_node))
        while len(successors) < 2 and len(predecessors) < 2 and successors:
            current_node = successors[0]
            path.append(current_node)
            successors = list(graph.successors(current_node))
            predecessors = list(graph.predecessors(current_node))
        tips.append(path)
    path_lengths = [len(path) for path in tips]
    avg_path_weights = [path_average_weight(graph, path) for path in tips]
    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_entry_node=True)

    return graph

def solve_out_tips(graph, ending_nodes):
    graph = simplify_bubbles(graph)
    tips = []
    for sink_node in ending_nodes:
        current_node = sink_node
        path = [current_node]
        successors = list(graph.successors(current_node))
        predecessors = list(graph.predecessors(current_node))
        while len(successors) < 2 and len(predecessors) < 2 and predecessors:
            current_node = predecessors[0]
            path.append(current_node)
            successors = list(graph.successors(current_node))
            predecessors = list(graph.predecessors(current_node))
        tips.append(path[::-1])

    path_lengths = [len(path) for path in tips]
    avg_path_weights = [path_average_weight(graph, path) for path in tips]
    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_sink_node=True)
    return graph


def draw_graph(graph, graphimg_file):
    """Dessin du graphe
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)





#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    print("kmer size: ", args.kmer_size)
    #print(build_kmer_dict(args.fastq_file, args.kmer_size))
    print("number of kmers : ", len(build_kmer_dict(args.fastq_file, args.kmer_size)))
    graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))
    

    starting_nodes = get_starting_nodes(graph)
    exit_nodes = get_sink_nodes(graph)
    print(starting_nodes)
    print(exit_nodes)


if __name__ == '__main__':
    main()
