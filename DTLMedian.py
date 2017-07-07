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
#   Scored DTL reconciliation graph:
#       { mapping_node: {event1:score, event2:score, ...} ...}
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

import DTLReconGraph as RG
import NewDiameter
from operator import itemgetter


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
    geneMultiplier = len(mappingNodeList)

    for i, node in enumerate(preorderGeneNodeList):
        geneLevelLookup[node] = i * geneMultiplier

    for i, node in enumerate(preorderSpeciesList):
        speciesLevelLookup[node] = i

    # The lambda function looks up the
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

    # Initialize the memo
    memo = dict()

    # Initialize the very start count, for the first call of countMPRs
    count = 0

    # Loop over all given minimum cost reconciliation roots
    for mappingNode in preorderMappingNodeList:
        if mappingNode[0] == geneRoot:
            count += countMPRs(mappingNode, DTLReconGraph, memo)

    # Initialise the scores dict. This dict contains the frequency score of each
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
            scores[mappingNode] = memo[mappingNode] / float(count)
        calculateScoresForChildren(mappingNode, DTLReconGraph, eventScores, scores, memo)

    return eventScores, count


def countMPRs(mappingNode, DTLReconGraph, memo):
    """
    :param mappingNode: an individual mapping node that maps a node
    for the parasite tree onto a node of the host tree, in the format
    (p, h), where p is the parasite node and h is the host node
    :param DTLReconGraph: A DTL reconciliation graph (see data structure quick reference at top of file)
    :param memo: a dictionary representing the running memo that is passed
    down recursive calls of this function. At first it is just an empty
    dictionary (see above function), but as it gets passed down calls, it collects
    keys of mapping nodes and values of MPR counts. This memo improves runtime
    of the algorithm
    :return: the number of MPRs spawned below the given mapping node in the graph
    """

    # Search the memo dictionary for previously calculated results
    if mappingNode in memo:
        return memo[mappingNode]

    # Base case, occurs if being called on a child produced by a loss or contemporary event
    if mappingNode == (None, None):
        return 1

    # Initialize a variable to keep count of the number of MPRs
    count = 0

    # Loop over all event nodes corresponding to the current mapping node
    for eventNode in DTLReconGraph[mappingNode][:-1]:

        # Save the children produced by the current event
        mappingChild1 = eventNode[1]
        mappingChild2 = eventNode[2]

        # Add the product of the counts of both children (over all children) for this event to get the parent's count
        memo[eventNode] = countMPRs(mappingChild1, DTLReconGraph, memo) * countMPRs(mappingChild2, DTLReconGraph, memo)
        count += memo[eventNode]

    # Save the result in the memo
    memo[mappingNode] = count

    return count


def calculateScoresForChildren(mappingNode, DTLReconGraph, scoredGraph, scores, memo):
    """
    This function calculates the frequency score for every mapping node that is a child of an event node that is a
    child of the given mapping node, and stores them in scoredGraph.
    :param mappingNode:     The mapping node that is the parent of the two scores we will compute
    :param DTLReconGraph:   The DTL reconciliation graph (see data structure quick reference at top of file)
    :param scoredGraph:     The scored DTL reconciliation graph (see data structure quick reference at top of file)
    :param scores:          The score for each mapping node (which will ultimately be thrown away) that this function
                            helps build up
    :param memo:            The counts generated in countMPRs
    :return:                Nothing, but scoredGraph is built up.
    """
    events = DTLReconGraph[mappingNode]

    assert scores[mappingNode] != 0
    # This multiplier is arcane bullshit that we all immediately forgot how it works, but it gets the job done.
    multiplier = float(scores[mappingNode])/memo[mappingNode]

    # Iterate over every event
    for eventNode in events[:-1]:

        scoredGraph[eventNode] = multiplier * memo[eventNode]

        # Save the children produced by the current event
        mappingChild1 = eventNode[1]
        mappingChild2 = eventNode[2]
        scores[mappingChild1] += scoredGraph[eventNode]
        scores[mappingChild2] += scoredGraph[eventNode]


def t(file="le1"):
    host, paras, graph, count, _ = RG.reconcile(file, 2, 3, 1)
    preorderGeneList, geneRoot, _ = NewDiameter.reformat_tree(paras, "pTop")
    preorderSpeciesList, speciesRoot, _ = NewDiameter.reformat_tree(host, "hTop")
    preorderGeneList = reversed(preorderGeneList)  # Actually postorder before being reversed
    preorderSpeciesList = reversed(preorderSpeciesList)  # Same as above
    preorderMappingNodeList = preorderMappingNodeSort(preorderGeneList, preorderSpeciesList, graph.keys())
    return generateScores(preorderMappingNodeList, graph, geneRoot)


# TODO: fix bug
def construct_median(DTLDict):
    """"""


def findMedian(DTLDict, postorderMappingNodes, MPRRoots):
    """
    :param DTLDict: A dictionary representing a DTL Recon Graph, but in a
    different format than what would be returned in DTLReconGraph.py. This should be a nested dictionary
    (as returned in the first value by generateScores), with keys being mapping nodes, and the values being
    other dictionaries. Each of these inner dictionaries has event nodes as keys (which should have a key for
    every event node corresponding to the previously indexed mapping node within a DTL Recon Graph) and
    frequency scores as float values, in the range (0, 1].
    :param postorderMappingNodes: A list of the mapping nodes in a possible MPR, except sorted first in
    postorder by species node and postorder by gene node
    :param MPRRoots: A list of mapping nodes that could act as roots to an MPR for the species and
    gene trees in question, output from the findBestRoots function in DTLReconGraph.py
    :return: A new dictionary which is has the same form as a DTL reconciliation graph except every
    mapping node only has one event node. Thus, this graph will represent a single reconciliation: the
    median reconciliation.
    """
    print(DTLDict)
    print('')
    # Note that for symmetric median reconciliation, each frequency must have 0.5 subtracted from it

    # Initialize a dict that will store the running total frequency sum incurred up to the given mapping node,
    # and the event node that directly gave it that frequency sum. Keys are mapping nodes, values are tuples
    # consisting of an event node and the corresponding running total frequency sum up to that mapping node
    sum_freqs = dict()

    # Initialize the dict that will represent the passed dict in a form that can be used by buildDTLReconGraph from RG
    DTLReconDict = dict()

    # Initialize the variable to store the number of median reconciliations
    n_med_recons = 0

    # Loop over all mapping nodes for the gene tree
    for map_node in postorderMappingNodes:

        # Get the events for the current mapping node and their frequencies, in a tuple in that order
        events = DTLDict[map_node].items()

        # Easily find the event with the best frequency
        best = max(events, key=itemgetter(1))

        # Consider the children - and thus future event nodes - of the event
        if best[0][0] == 'C':  # Check for a contemporaneous event
            sum_freqs[map_node] = (best[0], 0.5)  # C events have freq 1, so 1 - 0.5 = 0.5
        elif best[0][0] == 'L':  # Losses are also special since they produce only 1 child
            sum_freqs[map_node] = (best[0], best[1] - 0.5 + sum_freqs[best[0][1]][1])  # Consider the child's score too
        else:  # The only other cases to consider are speciation, duplication, or transfer, which have 2 children
            sum_freqs[map_node] = (best[0], best[1] - 0.5 + sum_freqs[best[0][1]][1] + sum_freqs[best[0][2]][1])

    # Get all possible roots of the graph, and their running frequency scores, in a list, for later use
    possible_root_combos = [(root, sum_freqs[root][1]) for root in MPRRoots]

    # Find the best root and frequency combo
    best = max(possible_root_combos, key=itemgetter(1))

    # Find the number of median reconciliations
    for root in possible_root_combos:

        # Best[1] is the best running total frequency found in the roots
        if root[1] == best[1]:
            n_med_recons += 1

    # Now extract the actual root
    # Note we convert it to a list so we can use it as input for buildDTLReconGraph from the file of the same name
    med_root = [best[0]]

    # Find the best root, of possible MPR start roots, to start on that gives the median
    print(best[1])
    print('')
    # Convert the given DTLDict to a form that will play nice with DTLReconGraph
    for mapping_node in DTLDict.keys():

        # We need the garbage 0 value at the end for the dictionary to work with buildDTLReconGraph
        DTLReconDict[mapping_node] = list(DTLDict[mapping_node].keys()) + [0]

    # Build the final median reconciliation graph
    med_recon_graph = RG.buildDTLReconGraph(med_root, DTLReconDict, {})

    return med_recon_graph, n_med_recons


def s(filename='example'):
    h, p, g, _, bestRoots = RG.reconcile(filename, 2, 3, 1)
    new_ptree = NewDiameter.reformat_tree(p, 'pTop')[0]
    new_htree = NewDiameter.reformat_tree(h, 'hTop')[0]

    result = findMedian(t(filename)[0], preorderMappingNodeSort(new_ptree, new_htree, g.keys()), bestRoots)

    return result
