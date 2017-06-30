# DTLReconGraph.py  --  was previously named DP.py
# Ran Libeskind-Hadas, June 2015
# The basic DP algorithm for reconciling pairs of trees

# Altered and expanded by Carter Slocum and Annalise Schweickart
# Altered and expanded by Andrew Ramirez and Eli Zupke

# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".
# The final DTL reconciliation graph is represented as a dictionary, with
# keys representing mapping nodes, and the associated values corresponding
# to events that can occur at each mapping node in an MPR, as well as other small
# values/information (see docstring of DP for more info on these values).

# Edited by Annalise Schweickart and Carter Slocum, July 2015 to return
# the DTL reconciliation graph that uses frequency scoring, as well as the
# number of reconciliations of the host and parasite trees.

# Edited by Andrew Ramirez and Eli Zupke in late June 2017 to return the number of MPRs for a given data set
# with an improved time efficiency relative to previous methods, usually by a factor of greater than or
# equal to 2.

# All necessary import statements
import copy
import newickFormatReader
import Greedy

Infinity = float('inf')


def preorder(tree, rootEdgeName):
    """
    :param tree: host or parasite tree (see description above)
    :param rootEdgeName: the name associated with the root of the given tree
    :return: list of edges in the given tree in preorder (high to low edges). See
    tech report for more information on post- or pre-order.
    """

    value = tree[rootEdgeName]
    _, _, leftChildEdgeName, rightChildEdgeName = value

    # Base case
    if leftChildEdgeName is None:  # Then rightChildEdgeName == None also
        return [rootEdgeName]
    else:  # Recursive call
        return [rootEdgeName] + \
                preorder(tree, leftChildEdgeName) + \
                preorder(tree, rightChildEdgeName)


def postorder(tree, rootEdgeName):
    """ The parameters of this function are the same as that of preorder above, except it
    returns the edge list in postorder (low to high edges; see tech report for more
    info on post- or pre-order)."""

    value = tree[rootEdgeName]
    _, _, leftChildEdgeName, rightChildEdgeName = value

    # Base case
    if leftChildEdgeName is None:  # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:  # Recursive call
        return postorder(tree, leftChildEdgeName) + \
               postorder(tree, rightChildEdgeName) + \
               [rootEdgeName]


def DP(hostTree, parasiteTree, phi, D, T, L, useScores=False):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, and duplication cost (D),
        transfer cost (T), loss cost (L), and a bool indicating whether score frequency (defined
        later in docstring) should be outputted, and returns the DTL reconciliation
        graph in the form of a dictionary, as well as the total cost of the
        optimal Maximum Parsimony Reconciliation and the number of maximum
        parsimony reconciliations. Note that the DTL reconciliation graph
        is returned as a dictionary with mapping nodes for keys, and values
        corresponding to lists which include all valid event nodes (along with
        the number of MPRs in which each event occurs, the last element in each
        tuple in the list) for a given mapping node for the MPR, minimum cost associated
        with the mapping node as the last element in the list. The frequency score
        corresponds to the number of MPRs in which eachThe notation and dynamic programming
        algorithm are explained in the tech report. Cospeciation is assumed to cost 0. """

    # A, C, O, and bestSwitch are all defined in tech report. Keys are edges and values are as defined in tech report
    A = {}
    C = {}
    O = {}
    bestSwitch = {}

    # Keeps track of events, children, and scores. Keys are vertex pairs and values are event representations in a tuple
    eventsDict = {}

    # Dictionary to keep track of minimum reconciliation cost for each (vp, vh)
    minCost = {}

    # Keeps track of which vertex mappings 'gave' O its cost for the corresponding edges
    oBest = {}

    # Keeps track of switch locations. Keys are edges, values are edges to send the key edges to for transfers
    bestSwitchLocations = {}

    # Calculates and keeps track of the frequency scoring of each event
    Score = {}

    # Following logic taken from tech report, we loop over all ep and eh
    for ep in postorder(parasiteTree, "pTop"):

        # Get the parasite tree info in the format
        # (vp top, vp bottom, edge of child 1, edge of child 2)
        _, vp, ep1, ep2 = parasiteTree[ep]

        # If there's no child 1, there's no child 2 and vp is a tip
        if ep1 is None:
            vpIsATip = True
            pChild1 = None
            pChild2 = None
        else:
            vpIsATip = False

            # Save end node names for the parasite's children
            pChild1 = ep1[1]
            pChild2 = ep2[1]

        # Begin looping over host edges
        for eh in postorder(hostTree, "hTop"):

            # Similar format to that of the parasite tree above
            _, vh, eh1, eh2 = hostTree[eh]

            # Initialize entries for this iteration of ep and eh
            eventsDict[(vp, vh)] = []
            oBest[(vp, vh)] = []

            # Same logic as for the parasite tree above
            if eh1 is None:
                vhIsATip = True
                hChild1 = None
                hChild2 = None
            else:
                vhIsATip = False

                # Save end node names for the host's children
                hChild1 = eh1[1]
                hChild2 = eh2[1]
                
            # Compute A(ep, eh)

            if vhIsATip:

                # Check if the tips map to one another
                if vpIsATip and phi[vp] == vh:

                    # The cost of matching mapped tips (thus, their edges) is 0
                    A[(ep, eh)] = 0

                    # Create a contemporary event
                    Amin = [("C", (None, None), (None, None), 1.0)]

                    # Give a frequency of 1 to this event
                    Score[(vp, vh)] = 1.0
                else:

                    # Non-matched tips can't reconcile
                    Score[(vp, vh)] = Infinity
                    A[(ep, eh)] = Infinity
                    Amin = [Infinity]
            else:

                # Compute Co and create event list to add to eventsDict

                if not vpIsATip:

                    # Calculate cospeciation cost assuming the cost is 0
                    COepeh = min(C[(ep1, eh1)] + C[(ep2, eh2)],
                                 C[(ep1, eh2)] + C[(ep2, eh1)])
                    coMin = []  # List to keep track lowest cost speciation
                    if COepeh == C[(ep2, eh1)] + C[(ep1, eh2)]:
                        coMin.append(("S", (pChild2, hChild1),
                                      (pChild1, hChild2), (Score[(pChild2, hChild1)] *
                                                           Score[(pChild1, hChild2)])))
                    if COepeh == C[(ep1, eh1)] + C[(ep2, eh2)]:
                        coMin.append(("S", (pChild1, hChild1),
                                      (pChild2, hChild2), (Score[(pChild1, hChild1)] *
                                                           Score[(pChild2, hChild2)])))
                else:
                    COepeh = Infinity
                    coMin = [Infinity]
                    Score[(vp, vh)] = Infinity

                # Compute L and create event list to add to eventsDict
                LOSSepeh = L + min(C[(ep, eh1)], C[(ep, eh2)])
                lossMin = []  # List to keep track of lowest cost loss

                # Check which (or maybe both) option produces the minimum
                if LOSSepeh == L + C[(ep, eh1)]:
                    lossMin.append(("L", (vp, hChild1), (None, None),
                                    Score[(vp, hChild1)]))
                if LOSSepeh == L + C[(ep, eh2)]:
                    lossMin.append(("L", (vp, hChild2), (None, None),
                                    Score[(vp, hChild2)]))

                # Determine which event occurs for A[(ep, eh)]
                A[(ep, eh)] = min(COepeh, LOSSepeh)

                # Record event occurring for A[(ep, eh)] (as Amin) by seeing which
                # event(s) produces least cost
                if COepeh < LOSSepeh:
                    Amin = coMin
                elif LOSSepeh < COepeh: 
                    Amin = lossMin
                else:
                    Amin = lossMin + coMin

            # Compute C(ep, eh)
            # First, compute D
            if not vpIsATip:

                # Calculate the cost of a duplication event
                DUPepeh = D + C[(ep1, eh)] + C[(ep2, eh)]

                # List to keep track of lowest cost duplication event
                dupList = ("D", (pChild1, vh), (pChild2, vh),
                           (Score[(pChild1, vh)] * Score[(pChild2, vh)]))
            else:
                DUPepeh = Infinity
                dupList = [Infinity]

            # Next, Compute T and create event list to add
            # to eventsDict using bestSwitchLocations
            if not vpIsATip:
                switchList = []  # List to keep track of lowest cost switch

                # Calculate the cost of a switch/transfer event
                SWITCHepeh = T + min(C[(ep1, eh)] + bestSwitch[(ep2, eh)],
                                     C[(ep2, eh)] + bestSwitch[(ep1, eh)])

                # If ep2 switching has the lowest cost or equal to the other
                if (C[(ep1, eh)] + bestSwitch[(ep2, eh)]) <= (C[(ep2, eh)] +
                                                              bestSwitch[(ep1, eh)]):

                    # Search for the optimal switch location by searching through the best switch
                    # locations for the given child and vh pair
                    for location in bestSwitchLocations[(pChild2, vh)]:

                        # Proposed new landing site
                        currentLoc = location[1]

                        # Proposed switch to a leaf, which is an impossible event
                        if currentLoc is None:
                            Score[(pChild1, currentLoc)] = Infinity
                            Score[(pChild2, currentLoc)] = Infinity

                        # Append the proposed event to the list of possible switches
                        switchList.append(("T", (pChild1, vh), (pChild2,
                                           currentLoc), (Score[(pChild1, vh)] *
                                                         Score[(pChild2, currentLoc)])))

                # If ep1 switching has the lowest cost or equal to the other
                elif (C[(ep2, eh)] + bestSwitch[(ep1, eh)]) <= (C[(ep1, eh)] +
                                                                bestSwitch[(ep2, eh)]):

                    # Search for the optimal switch location by searching through the best switch
                    # locations for the given child and vh pair
                    for location in bestSwitchLocations[(pChild1, vh)]:

                        # Proposed new landing site
                        currentLoc = location[1]

                        # Proposed switch to a leaf, which is an impossible event
                        if currentLoc is None:
                            Score[(pChild1, currentLoc)] = Infinity
                            Score[(pChild2, currentLoc)] = Infinity

                        # Append the proposed event to the list of possible switches
                        switchList.append(("T", (pChild2, vh),
                                           (pChild1, currentLoc), (Score[(pChild2, vh)] *
                                                                   Score[(pChild1, currentLoc)])))
            else:  # vp is a tip
                SWITCHepeh = Infinity
                switchList = [Infinity]

            # Compute C[(ep, eh)] and add the event or events with that cost
            # to the dictionary eventsDict
            C[(ep, eh)] = min(A[(ep, eh)], DUPepeh, SWITCHepeh)

            # Add the minimum costs for the current edges to the minCost dict
            minCost[(vp, vh)] = C[(ep, eh)]

            # Find which events produce the optimal cost for these edges and add to eventDict
            if C[(ep, eh)] == DUPepeh:
                eventsDict[(vp, vh)].append(dupList)
            if C[(ep, eh)] == SWITCHepeh:
                eventsDict[(vp, vh)].extend(switchList)
            if C[(ep, eh)] == A[(ep, eh)]:
                eventsDict[(vp, vh)].extend(Amin)

            # Calculate O for eh's children

            # Scan through all of the keys (e.g. (a, A)) recorded so far
            # Note an event is stored in eventsDict in the form
            # [event (str), pair1 (tuple, str), pair2 (tuple, str), score (float)]
            for key in eventsDict:
                mapScore = 0  # Initialize frequency scoring for each event

                # Search the dict for events related to the current key
                for event in eventsDict[key]:

                    # This filters actual events, since events are stored as lists
                    if isinstance(event, tuple):

                        # Frequency scores are the last element, and that's what this extracts
                        mapScore += event[-1]

                # Accumulate that and set it as the score for the current key
                Score[key] = mapScore

            # Remove all 'impossible' events from the options
            if minCost[(vp, vh)] == Infinity:
                del minCost[(vp, vh)]
                del eventsDict[(vp, vh)]

            # Compute oBest[(vp, vh)], the source of O(ep, eh)
            if vhIsATip: 
                O[(ep, eh)] = C[(ep, eh)]
                oBest[(vp, vh)] = [(vp, vh)]              
            else:

                # Compute O(ep, eh) if vh is not a tip
                O[(ep, eh)] = min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)])

                # oMin helps us easily find which value (between C, O for child 1, and O for child 2) produces
                # O for this edge. Knowing what its indices represent, we search through to see which produce O
                oMin = [i for i, e in enumerate([C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)]])
                        if e == O[(ep, eh)]]

                # Corresponds to C
                if 0 in oMin:
                    oBest[(vp, vh)].append((vp, vh))

                # Corresponds to the O table for each child
                if 1 in oMin:
                    oBest[(vp, vh)].extend(oBest[(vp, hChild1)])
                if 2 in oMin:
                    oBest[(vp, vh)].extend(oBest[(vp, hChild2)])

        # Compute bestSwitch values
        bestSwitch[(ep, "hTop")] = Infinity
        bestSwitchLocations[(vp, hostTree["hTop"][1])] = [(None, None)]
        for eh in preorder(hostTree, "hTop"):

            # Redefine the host information for this new loop
            _, vh, eh1, eh2 = hostTree[eh]

            # Is vh a tip?
            if eh1 is None:  # Then eh2 == None too and vh is a tip!
                vhIsATip = True
                hChild1 = None
                hChild2 = None
            else:
                vhIsATip = False
                hChild1 = hostTree[eh][2][1]
                hChild2 = hostTree[eh][3][1]

            # Find best cost for a switch to occur (bestSwitch)
            # and the location to which the edge switches (bestSwitchLocations)   
            if not vhIsATip:

                # Initialize lists for switch locations
                bestSwitchLocations[(vp, hChild1)] = []
                bestSwitchLocations[(vp, hChild2)] = []

                # Compute the switch costs
                bestSwitch[(ep, eh1)] = min(bestSwitch[(ep, eh)], O[(ep, eh2)])
                bestSwitch[(ep, eh2)] = min(bestSwitch[(ep, eh)], O[(ep, eh1)])

                # Add best switch locations for child 1
                if bestSwitch[(ep, eh1)] == bestSwitch[(ep, eh)] and \
                   bestSwitchLocations[(vp, vh)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild1)].extend(
                        bestSwitchLocations[(vp, vh)])
                if bestSwitch[(ep, eh1)] == O[(ep, eh2)] and \
                   oBest[(vp, hChild2)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild1)].extend(
                        oBest[(vp, hChild2)])

                # Add best switch locations for child 2
                if bestSwitch[(ep, eh2)] == bestSwitch[(ep, eh)] and \
                   bestSwitchLocations[(vp, vh)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild2)].extend(
                        bestSwitchLocations[(vp, vh)])
                if bestSwitch[(ep, eh2)] == O[(ep, eh1)] and \
                   oBest[(vp, hChild1)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild2)].extend(
                        oBest[(vp, hChild1)])

    # Append the cost of each event to the corresponding eventsDict entry
    for key in eventsDict:
        eventsDict[key].append(minCost[key])

    # Create the list of minimum cost mapping nodes involving root of parasite tree
    treeMin = findBestRoots(parasiteTree, minCost)

    # Build the reconciliation graph as a dictionary, with keys as mapping nodes and values as event nodes
    DTLReconGraph = buildDTLReconGraph(treeMin, eventsDict, {})

    # Initialize the variable storing the number of MPRs
    nMPRs = None

    # Consider whether the user wanted to use the Scores
    if useScores:

        # Filter values from the score dictionary so it can be passed to addScores
        for key in Score.keys():
            if key not in DTLReconGraph:
                del Score[key]

        DTLReconGraph, nMPRs = addScores(treeMin, DTLReconGraph, Score)  # Less efficient method of computing nMPRs

        # Make sure nMPRs is represented as an int - the result will always be an integer value
        nMPRs = int(nMPRs)

    # If not, suppress score from the output
    else:
        for mappingNode in DTLReconGraph:  # Search over all nodes in the recon graph
            for event in range(len(DTLReconGraph[mappingNode]) - 1):  # Search over all events for each node
                # Note we do len(... - 1) because the last element isn't an event, it's the min cost (defined above)

                # Filter out spurious entries from the return value
                DTLReconGraph[mappingNode][event] = DTLReconGraph[mappingNode][event][:-1]

                # The total number of MPRs for this optimal cost, calculated the more time-efficient way
                nMPRs = countMPRsWrapper(treeMin, DTLReconGraph)

    # The total cost of the best reconciliation
    bestCost = minCost[treeMin[0]]

    # Returns the graph, the optimal cost, and the number of MPRs
    return DTLReconGraph, bestCost, nMPRs


def preorderDTLsort(DTLReconGraph, ParasiteRoot):
    """
    :param DTLReconGraph: one of the outputs from DP, directly outputted by buildDTLReconGraph (see
    top of file for structure of the DTLReconGraph)
    :param ParasiteRoot: The root node of the parasite tree, represented as a string
    :return: an ordered list of tuples, in order of increasing level. Each element is of the form ((P, H), L),
    where P is a parasite node and H is a host node, and the parasite node is mapped onto the host node in
    the first tuple in the overarching tuple. The second element is the level in the tree at which the (P,H)
    tuple occurs. Note level 0 is the root and the highest level represents tips.
    """

    keysL = Greedy.orderDTL(DTLReconGraph, ParasiteRoot)
    orderedKeysL = []  # We could marginally improve efficiency here by locking list length, but we don't do that here
    levelCounter = 0
    while len(orderedKeysL) < len(keysL):
        toAdd = []
        for mapping in keysL:
            if mapping[-1] == levelCounter:
                    toAdd += [mapping]
        orderedKeysL += toAdd
        levelCounter += 1

    return orderedKeysL


def preorderCheck(preOrderList):
    """
    :param preOrderList: output from preorderDTLsort. See that function for the structure
    of the input.
    :return: The same ordered list inputted as preOrderList, except duplicate tuples
    in the list have been removed
    """

    # newList = [(a1, a2), b] where a is the map node and b is the depth level
    # Filtering for multiple instances of a with different b by keeping biggest
    # b instance. This is safe: ensures that all possible parents of a node will
    # be handled before a node to prevent considering duplicates in the
    # addScores function.
    # Correction by Jean Sung, July 2016

    newList = []
    preDict = {}
    for root in preOrderList:
        if root not in newList:
            newList.append(root)
    for x in range(len(newList)):
        currentRoot = newList[x][0]
        currentLevel = newList[x][1]
        if currentRoot in preDict:
            if preDict[currentRoot][0] > currentLevel:
                newList[x] = (None, None)
            else:
                location = preDict[currentRoot][1]
                newList[location] = (None, None)
        else:
            preDict[currentRoot] = (currentLevel, x)

    finalList = []
    for item in newList:
        node = item[0]
        depth = item[1]
        finalList.append(item)

        # Check for duplicate instances with smaller depth levels 
        for newItem in finalList:
            newNode = newItem[0]
            newDepth = newItem[1]
            if (node == newNode) and (depth > newDepth):
                finalList.remove(newItem)
                break
    return finalList


def addScores(treeMin, DTLReconGraph, ScoreDict):
    """
    :param treeMin: output from findBestRoots. It is a list of mapping
    nodes that could produce a MPR, with elements in the format (p, h),
    where p is the parasite node being mapped onto host node h
    :param DTLReconGraph: output from buildDTLReconGraph (see that
    function for the structure of this input)
    :param ScoreDict: the dictionary of frequency scores, where each key
    is a tuple (vp, vh), where vp is a parasite vertex and vh is a host
    vertex, and values are frequency scores. See the description of Score
    near the top of the function DP for more info on this input
    :return: a new DTLReconGraph with only frequency scores of each event
    modified, but it also returns the number of MPRs for the given
    DTLReconGraph
    """

    # Dictionary has mapping nodes with event node children
    newDTL = copy.deepcopy(DTLReconGraph)
    parentsDict = {}
    preOrder1 = preorderDTLsort(DTLReconGraph, treeMin[0][0])
    preOrder2 = preorderCheck(preOrder1)
    for root in preOrder2:
        if root != (None, None):
            vertex = root[0]

            if root[1] == 0:
                parentsDict[vertex] = ScoreDict[vertex]

            # LAST is garbage  value used in the dp construction   
            for n in range(len(DTLReconGraph[vertex]) - 1):
                _, child1, child2, oldScore = DTLReconGraph[vertex][n]

                # Copy the tuple but modify the last value - this is the only way to do it, since tuples are immutable
                newDTL[vertex][n] = (newDTL[vertex][n][0], newDTL[vertex][n][1], newDTL[vertex][n][2],
                                     parentsDict[vertex] * (oldScore / ScoreDict[vertex]))

                if child1 != (None, None):
                    if child1 in parentsDict:
                        parentsDict[child1] += newDTL[vertex][n][3]
                    else: 
                        parentsDict[child1] = newDTL[vertex][n][3] 
                if child2 != (None, None):
                    if child2 in parentsDict:
                        parentsDict[child2] += newDTL[vertex][n][3]
                    else: 
                        parentsDict[child2] = newDTL[vertex][n][3]
     
    normalize = newDTL[preOrder2[-1][0]][0][-1]  # preOrder2[-1][0] is the last value considered in the above loop

    # Adjust all values in the DTL
    for key in newDTL:

        # Again, last value is garbage 
        for event in newDTL[key][:-1]:

            newDTL[key][:-1][-1] = event[-1] / normalize

    return newDTL, normalize


def countMPRsWrapper(mappingNodeLst, DTLReconGraph):
    """
    :param mappingNodeLst: output from findBestRoots, a list of mapping
    nodes for the root of the parasite tree that could produce a MPR.
    See findBestRoots for more info on this input.
    :param DTLReconGraph: the output from buildDTLReconGraph. See that
    function for more info on this input.
    :return: this function uses the helper function countMPRs to loop over
    all of the minimum cost parasite root mappings and sum their MPR counts
    to find the total number of MPRs for the given DTLReconGraph. This
    number is returned as an integer
    """

    # Initialize the memo
    memo = dict()

    # Initialize the very start count, for the first call of countMPRs
    count = 0

    # Loop over all given minimum cost reconciliation roots
    for mappingNode in mappingNodeLst:
        count += countMPRs(mappingNode, DTLReconGraph, memo)

    return count


def countMPRs(mappingNode, DTLReconGraph, memo):
    """
    :param mappingNode: an individual mapping node that maps a node
    for the parasite tree onto a node of the host tree, in the format
    (p, h), where p is the parasite node and h is the host node
    :param DTLReconGraph: a DTLReconGraph, output from buildDTLReconGraph
    (see that function for more info on the format of this input)
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

    # Base case, occurs if being called on a child produced by a loss or contemporary evet
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
        count += countMPRs(mappingChild1, DTLReconGraph, memo) * countMPRs(mappingChild2, DTLReconGraph, memo)

    # Save the result in the memo
    memo[mappingNode] = count

    return count


def findBestRoots(parasiteTree, minCostDict):
    """
    :param parasiteTree: parasite tree in the format described at the
    top of this file
    :param minCostDict: a dictionary - keys representing mappings of
    parasite vertices onto host vertices (p, h) and values representing
    the minimum cost for a reconciliation based on that mapping. In DP,
    this dictionary was called minCost
    :return: returns a list of all of the mapping nodes involving the
    parasite root that can produce a MPR
    """
    treeTops = []
    for key in minCostDict:
        if key[0] == parasiteTree['pTop'][1]:
            treeTops.append(key)
    treeMin = []
    min_score = min([minCostDict[root] for root in treeTops])
    for pair in treeTops:
        if minCostDict[pair] == min_score:
            treeMin.append(pair)
    return treeMin


def buildDTLReconGraph(tupleList, eventDict, uniqueDict):
    """
     :param tupleList: a list of minimum cost reconciliation roots - see findBestRoots
     for more info on the format of this input
     :param eventDict: a dictionary representing events and the corresponding children
     for each node - see eventDict in DP for more info on the format of this input
     :param uniqueDict: a dictionary of unique vertex mappings, which initially
     starts empty and gets built up using eventDict and tupleList using recursion
     :return: the modified uniqueDict for this particular call of the function
     """
    for vertexPair in tupleList:
        if vertexPair not in uniqueDict:
            uniqueDict[vertexPair] = eventDict[vertexPair]
        for event in eventDict[vertexPair][:-1]:
            for location in event:
                if type(location) is tuple and location != (None, None):
                    buildDTLReconGraph([location], eventDict, uniqueDict)
    return uniqueDict


def reconcile(fileName, D, T, L, useScores=False):
    """
    :param fileName: the file in which the desired data set it stored, passed as
    a string. For Ran Libesking-Hadas's group, our data files were almost exclusively
    .newick files once we were sure our algorithm worked correctly, which needed to use
    the newick format reader to correctly read in the data.
    :param D: the cost associated with a duplication event
    :param T: the cost associated with a transfer event
    :param L: the cost associated with a loss event
    :param useScores: an boolean value that lets the user choose whether to use the
    more efficient method of finding the number of MPRs for a given data set (useScores=False)
    or to use the less efficient method that still outputs the number of MPRs but also
    gives frequency scores in the values for the DTLReconGraph
    :return: the host tree used, the parasite tree used, the DTLReconGraph, and the number of
    MPRs (as an int). See preceding functions for details on the format of the host and parasite
    trees as well as the DTLReconGraph
    """
    # Note: I have made modifications to the return statement to make Diameter.py possible without re-reconciling.
    host, paras, phi = newickFormatReader.getInput(fileName)
    graph, _, numRecon = DP(host, paras, phi, D, T, L, useScores)
    return host, paras, graph, numRecon
