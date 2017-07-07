# ReconGraphProperties.py
# Written July 2017 by Andrew Ramirez and Eli Zupkes

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
import Diameter

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
    scores[(None,None)] = 0.0

    # The scored graph is like the DTLReconGraph, except instead of individual events being in a list, they are the
    # keys of a dictionary where the values are the frequency scores of those events.
    scoredGraph = {}

    for mappingNode in preorderMappingNodeList:
        scoredGraph[mappingNode] = {}

        # If we are at the root of the gene tree, then we need to initialise the score entry
        if mappingNode[0] == geneRoot:
            scores[mappingNode] = memo[mappingNode] / float(count)
        calculateScoresForChildren(mappingNode, DTLReconGraph, scoredGraph, scores, memo)
    for thing in scores:
        if len(thing) > 2:
            print "{0} : {1}".format(thing, scores[thing])

    return scoredGraph, count



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

        scoredGraph[mappingNode][eventNode] = multiplier * memo[eventNode]

        # Save the children produced by the current event
        mappingChild1 = eventNode[1]
        mappingChild2 = eventNode[2]
        scores[mappingChild1] += scoredGraph[mappingNode][eventNode]
        scores[mappingChild2] += scoredGraph[mappingNode][eventNode]



def t(file="le1"):
    host, paras, graph, count = RG.reconcile(file,1,4,1)
    preorderGeneList, geneRoot, _ = Diameter.reformat_tree(paras,"pTop")
    preorderSpeciesList, speciesRoot, _ = Diameter.reformat_tree(host,"hTop")
    preorderGeneList = reversed(preorderGeneList)
    preorderSpeciesList = reversed(preorderSpeciesList)
    preorderMappingNodeList = preorderMappingNodeSort(preorderGeneList, preorderSpeciesList,  graph.keys())
    scores = generateScores(preorderMappingNodeList, graph, geneRoot)
    for k in scores[0]:
        print "{0} : {1}".format(k, scores[0][k])
    print count