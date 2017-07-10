# DTLMedian.py
# Written July 2017 by Andrew Ramirez and Eli Zupke
# Utilizes previous code to find the "median" reconciliation of a pair of gene and species trees


# -1. DATA STRUCTURE QUICK REFERENCE:
#
#
#   DTL Reconciliation graph:
#       { mapping_node: [event1, event2, ... eventn, number] ...}
#   Event:
#       ('type', child_mapping_node1, child_mapping_node2)
#
#   Mapping node:
#       ('gene_node','SPECIES_NODE')
#   or in loss or contemporary event nodes:
#       (None, None)
#
#
#   (edge) trees:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#       aka:
#       {root_edge: (root_edge[0], root_edge[1], child1_edge, child2_edge) ...}
#
#   vertex_trees:
#       {'N':('C1','C2') ...}
#

from operator import itemgetter
import pandas as pd
import matplotlib.pyplot as plt
import DTLReconGraph as RG
import NewDiameter


def preorderMappingNodeSort(preorderGeneNodeList, preorderSpeciesList, mappingNodeList):
    """Sorts a list of mapping nodes into double preorder (the gene node is more significant than the species node, and)
    both are in preorder in the final list."""

    # In order to sort the mapping nodes, we need a way to convert them into numbers. These two lookup tables allow
    # us to achieve a lexicographical ordering with gene nodes more significant than species nodes.
    geneLevelLookup = {}
    speciesLevelLookup = {}

    # By multiplying the gene node keys by the number of species nodes, we can ensure that a mapping node with a later
    # gene node always comes before one with a later species node, because a gene node paired with the last species
    # node will be one level less than the next gene node paired with the first species node.
    geneMultiplier = len(preorderSpeciesList)

    for i, node in enumerate(preorderGeneNodeList):
        geneLevelLookup[node] = i * geneMultiplier

    for i, node in enumerate(preorderSpeciesList):
        speciesLevelLookup[node] = i

    # The lambda function looks up the level of both the gene node and the species nodes and adds them together to
    # get a number to give to the sorting algorithm for that mapping node. The gene node is weighted far more heavily
    # than the species node to make sure it is always more significant.
    sortedList = sorted(mappingNodeList, key=lambda node: geneLevelLookup[node[0]] + speciesLevelLookup[node[1]])

    return sortedList


def generateScores(preorderMappingNodeList, DTLReconGraph, geneRoot):
    """
    :param preorderMappingNodeList: A list of all mapping nodes in DTLReconGraph in double preorder
    :param DTLReconGraph:           The DTL reconciliation graph that we are scoring
    :param geneRoot:                The root of the gene tree
    :return:                        0. A file structured like the DTLReconGraph, but with the lists of events replaced
                                    with dicts, where the keys are the events and the values are the scores of those
                                    events, and
                                    1. The number of MPRs in DTLReconGraph.
    """

    # Initialize the dictionary that will store mapping node and event counts (which also acts as a memoization
    # dictionary)
    counts = dict()

    # Initialize the very start count, for the first call of countMPRs
    count = 0

    # Loop over all given minimum cost reconciliation roots
    for mappingNode in preorderMappingNodeList:
        if mappingNode[0] == geneRoot:
            count += countMPRs(mappingNode, DTLReconGraph, counts)

    # Initialize the scores dict. This dict contains the frequency score of each
    scores = dict()
    for mappingNode in preorderMappingNodeList:
        scores[mappingNode] = 0.0

    # This entry is going to be thrown away, but it seems neater to just let calculateScoresOfChildren
    # add scores to an unused entry than to check to see if they are (None, None) in the first place.
    scores[(None, None)] = 0.0

    # The scored graph is like the DTLReconGraph, except instead of individual events being in a list, they are the
    # keys of a dictionary where the values are the frequency scores of those events.
    eventScores = {}

    for mappingNode in preorderMappingNodeList:

        # If we are at the root of the gene tree, then we need to initialise the score entry
        if mappingNode[0] == geneRoot:
            scores[mappingNode] = counts[mappingNode] / float(count)
        calculateScoresForChildren(mappingNode, DTLReconGraph, eventScores, scores, counts)

    return eventScores, count


def countMPRs(mappingNode, DTLReconGraph, counts):
    """
    :param mappingNode: an individual mapping node that maps a node
    for the parasite tree onto a node of the host tree, in the format
    (p, h), where p is the parasite node and h is the host node
    :param DTLReconGraph: A DTL reconciliation graph (see data structure quick reference at top of file)
    :param counts: a dictionary representing the running memo that is passed
    down recursive calls of this function. At first it is just an empty
    dictionary (see above function), but as it gets passed down calls, it collects
    keys of mapping nodes and values of MPR counts. This memo improves runtime
    of the algorithm
    :return: the number of MPRs spawned below the given mapping node in the graph
    """

    # Search the counts dictionary for previously calculated results (this is the memoization)
    if mappingNode in counts:
        return counts[mappingNode]

    # Base case, occurs if being called on a child produced by a loss or contemporary event
    if mappingNode == (None, None):
        return 1

    # Initialize a variable to keep count of the number of MPRs
    count = 0

    # Loop over all event nodes corresponding to the current mapping node
    for eventNode in DTLReconGraph[mappingNode]:

        # Save the children produced by the current event
        mappingChild1 = eventNode[1]
        mappingChild2 = eventNode[2]

        # Add the product of the counts of both children (over all children) for this event to get the parent's count
        counts[eventNode] = countMPRs(mappingChild1, DTLReconGraph, counts) * countMPRs(mappingChild2, DTLReconGraph, counts)
        count += counts[eventNode]

    # Save the result in the counts
    counts[mappingNode] = count

    return count


def calculateScoresForChildren(mappingNode, DTLReconGraph, scoredGraph, scores, counts):
    """
    This function calculates the frequency score for every mapping node that is a child of an event node that is a
    child of the given mapping node, and stores them in scoredGraph.
    :param mappingNode:     The mapping node that is the parent of the two scores we will compute
    :param DTLReconGraph:   The DTL reconciliation graph (see data structure quick reference at top of file)
    :param scoredGraph:     The scored DTL reconciliation graph (see data structure quick reference at top of file)
    :param scores:          The score for each mapping node (which will ultimately be thrown away) that this function
                            helps build up
    :param counts:          The counts generated in countMPRs
    :return:                Nothing, but scoredGraph is built up.
    """
    events = DTLReconGraph[mappingNode]

    assert scores[mappingNode] != 0
    # This multiplier is arcane magic that we all immediately forgot how it works, but it gets the job done.
    multiplier = float(scores[mappingNode]) / counts[mappingNode]
    # Iterate over every event
    for eventNode in events:

        scoredGraph[eventNode] = multiplier * counts[eventNode]

        # Save the children produced by the current event
        mappingChild1 = eventNode[1]
        mappingChild2 = eventNode[2]
        scores[mappingChild1] += scoredGraph[eventNode]
        scores[mappingChild2] += scoredGraph[eventNode]

def findMedian(DTLReconGraph, eventScores, postorderMappingNodes, MPRRoots):
    """
    :param DTLReconGraph: A dictionary representing a DTL Recon Graph.
    :param eventScores: A dictionary with event nodes as keys and values corresponding to the frequency of
    that events in MPR space for the recon graph
    :param postorderMappingNodes: A list of the mapping nodes in a possible MPR, except sorted first in
    postorder by species node and postorder by gene node
    :param MPRRoots: A list of mapping nodes that could act as roots to an MPR for the species and
    gene trees in question, output from the findBestRoots function in DTLReconGraph.py
    :return: A new dictionary which is has the same form as a DTL reconciliation graph except every
    mapping node only has one event node, along with the number of median reconciliations for the given DTL
    reconciliation graph, as well as the root of the median MPR for the given graph. Thus, this graph will
    represent a single reconciliation: the median reconciliation.
    """

    # Note that for symmetric median reconciliation, each frequency must have 0.5 subtracted from it

    # Initialize a dict that will store the running total frequency sum incurred up to the given mapping node,
    # and the event node that directly gave it that frequency sum. Keys are mapping nodes, values are tuples
    # consisting of an event node and the corresponding running total frequency sum up to that mapping node
    sum_freqs = dict()

    # Initialize the variable to store the number of median reconciliations
    n_med_recons = 0

    # Loop over all mapping nodes for the gene tree
    for map_node in postorderMappingNodes:

        # Get the events for the current mapping node and their frequencies, in a tuple in that order
        events = [(event, eventScores[event]) for event in DTLReconGraph[map_node]]

        # Easily find the event with the best frequency
        best_event = max(events, key=itemgetter(1))

        # Consider the children - and thus future event nodes - of the event
        if best_event[0][0] == 'C':  # Check for a contemporaneous event
            sum_freqs[map_node] = (best_event[0], 0.5)  # C events have freq 1, so 1 - 0.5 = 0.5
        elif best_event[0][0] == 'L':  # Losses are also special since they produce only 1 child
            sum_freqs[map_node] = (best_event[0], best_event[1] - 0.5 + sum_freqs[best_event[0][1]][1])
        else:  # The only other cases to consider are speciation, duplication, or transfer, which have 2 children
            sum_freqs[map_node] = (best_event[0], best_event[1] - 0.5 + sum_freqs[best_event[0][1]][1] +
                                   sum_freqs[best_event[0][2]][1])

    # Get all possible roots of the graph, and their running frequency scores, in a list, for later use
    possible_root_combos = [(root, sum_freqs[root][1]) for root in MPRRoots]

    # Find the best root and frequency combo
    best_root = max(possible_root_combos, key=itemgetter(1))

    # Find the number of median reconciliations
    for root in possible_root_combos:

        # best[1] is the best running total frequency found in the roots
        if root[1] == best_root[1]:
            n_med_recons += 1

    # Now extract the actual root
    # Note we convert it to a list so we can use it as input for buildDTLReconGraph from the file of the same name
    med_root = best_root[0]

    # Adjust the sum_freqs dictionary so we can use it with the buildMedianReconGraph function
    for root in sum_freqs:

        # We place the event tuples into lists so they work well with the diameter algorithm
        sum_freqs[root] = [sum_freqs[root][0]]  # Only use the event, no longer the associated frequency sum

    med_recon_graph = buildMedianReconGraph(sum_freqs, med_root)

    return med_recon_graph, n_med_recons, med_root


def buildMedianReconGraph(eventDict, root):
    """
    :param eventDict: a dictionary with mapping nodes for keys and values which are the single event that mapping
    node may have in a median reconciliation, as a tuple but each of these tuples are the single event in a list.
    :param root: the mapping node at which the median reconciliation or a subgraph of the median
    is starting at
    :return: a DTL Reconciliation Graph in the form returned in DTLReconGraph.py, except here the only
    reconciliation represented is the median - i.e., only events and mapping nodes valid in the median are
    represented
    """

    # Initialize the dict to be returned for this subgraph
    subgraph_recon_dict = dict()

    # From the get go, we need to save the current subgraph root and its event
    subgraph_recon_dict.update({root: eventDict[root]})

    # Check for a loss
    if eventDict[root][0][0] == 'L':
        subgraph_recon_dict.update(buildMedianReconGraph(eventDict, eventDict[root][0][1]))

    # Check for events that produce two children
    elif eventDict[root][0][0] in ['T', 'S', 'D']:
        subgraph_recon_dict.update(buildMedianReconGraph(eventDict, eventDict[root][0][1]))
        subgraph_recon_dict.update(buildMedianReconGraph(eventDict, eventDict[root][0][2]))

    return subgraph_recon_dict


# Test function
def t(filename='le1', D=2, T=3, L=1):

    # Get all of the important values to start our test
    species_tree, gene_tree, dtl_recon_graph, menpmn, mdenpmn, data, mpr_count, best_roots = RG.reconcile(filename, D, T, L)

    # Reformat gene tree and get info on it, as well as for the species tree in the following line
    postorder_gene_tree, gene_tree_root, gene_node_count = NewDiameter.reformat_tree(gene_tree, "pTop")
    postorder_species_tree, species_tree_root, species_node_count = NewDiameter.reformat_tree(species_tree, "hTop")

    # Get a list of the mapping nodes in preorder
    preorder_mapping_node_list = preorderMappingNodeSort(postorder_gene_tree, postorder_species_tree, dtl_recon_graph.keys())

    # Find the dictionary for frequency scores for the given mapping nodes and graph, as well as the given gene root
    scoresDict = generateScores(list(reversed(preorder_mapping_node_list)), dtl_recon_graph, gene_tree_root)
    median_reconciliation, n_meds, _ = findMedian(dtl_recon_graph, scoresDict[0], preorder_mapping_node_list, best_roots)

    # Clean up the reconciliation graph
    NewDiameter.clean_graph(dtl_recon_graph, gene_tree_root)

    # Use the diameter algorithm to find the diameter between the recon graph and its median
    diameter = NewDiameter.new_diameter_algorithm(postorder_species_tree, postorder_gene_tree, gene_tree_root, median_reconciliation, dtl_recon_graph, True, False)

    # Since this function contains a TON of information in its variables, any value can be changed to be the return
    # value to allow the user to get pretty much any info they need
    return data, menpmn, mdenpmn#diameter


# Test function to check the skewedness of the number of event nodes per mapping node
def checkSkew():

    # We want to loop through all files
    for i in range(1, 5666):

        # Start building the number of the tree of life data file
        filenum = str(i).zfill(4)

        try:
            data, datamean, datamedian = t('TreeLifeData/COG' + filenum + '.newick')
            plt.hist(data)
            plt.title('File: %s' % ('TreeLifeData/COG' + filenum + '.newick'))
            plt.ylabel('Count')
            plt.xlabel('Number of event nodes in a mapping node')
            print('mean: %s, median: %s' % (datamean, datamedian))
            plt.show()
            plt.clf()

        except IOError:
            print('File %s does not exist' % ('TreeLifeData/COG' + filenum + '.newick'))


# See how the mean number of event nodes per mapping node is distributed
def checkDistrOfMean():

    # Initialize our data list
    data = list()

    # This is just for us to check if there is ever such skew that the median is not 1
    nmediansnonzero = 0

    # We want to loop through all files
    for i in range(1, 5666):

        # Start building the number of the tree of life data file
        filenum = str(i)
        filenum = ('0' * (4 - len(filenum))) + filenum

        try:
            datapts, datamean, datamedian = t('TreeLifeData/COG' + filenum + '.newick')
            if datamedian != 1:
                nmediansnonzero += 1
            data.append(datamean)
        except IOError:
            print('File %s does not exist' % ('TreeLifeData/COG' + filenum + '.newick'))

    # Save the data to csv in case matplotlib formatting isn't good enough for whatever use necessary
    pd.DataFrame(data, columns=['menpmn']).to_csv('MeanEventNodesPerMappingNode.csv', index=False)

    print('%d data files had a median that was not 1' % nmediansnonzero)
    plt.hist(data)
    plt.title('Distribution of mean event nodes per mapping node')
    plt.ylabel('Count')
    plt.xlabel('Mean number of event nodes in a mapping node for a certain data file')
    plt.show()
