# DP.py
# Ran Libeskind-Hadas, June 2015
# The basic DP algorithm for reconciling pairs of trees

# Altered and expanded by Carter Slocum and Annalise Schweickart
# Altered and expanded by Andrew Ramirez and Eli Zupke

# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".
# TODO: describe the structure of the DTL reconciliation graph, i.e. dictionary with keys representing...
# TODO: ensure trees called trees where necessary and DTL graphs called DTL recon graphs where necessary
# TODO: place all "housekeeping" functions in this file so the necessary files are self-contained

# TODO: use IDE's special format for docstrings (e.g. params, return, etc.)

# Edited by Annalise Schweickart and Carter Slocum, July 2015 to return
# the DTL reconciliation graph that uses frequency scoring, as well as the
# number of reconciliations of the host and parasite trees.

# Edited by Andrew Ramirez in late June 2017 to return the number of MPRs for a given data set
# with a time efficiency improved by about a factor of 2 relative to previous
# methods for returning the number of MPRs.

import newickFormatReader
import time  # TODO: delete this import statement once done testing

Infinity = float('inf')


def preorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a 
    list of the edges in that tree in preorder (high edges to low edges)"""

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
    """ Takes a tree as input (see format description above) and returns a 
    list of the edges in that tree in postorder (low edges to high edges)"""

    value = tree[rootEdgeName]
    _, _, leftChildEdgeName, rightChildEdgeName = value

    # Base case
    if leftChildEdgeName is None:  # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:  # Recursive call
        return postorder(tree, leftChildEdgeName) + \
               postorder(tree, rightChildEdgeName) + \
               [rootEdgeName]


def DP(hostTree, parasiteTree, phi, D, T, L):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, and
        duplication cost (D), transfer cost (T), and loss cost (L) and
        returns the DTL reconciliation graph in the form of a dictionary, as well as the
        total cost of the optimal Maximum Parsimony Reconciliation and
        the number of maximum parsimony reconciliations. The notation and 
        dynamic programming algorithm are explained in the tech report.
        Cospeciation is assumed to cost 0. """

    # TODO: for every dictionary, describe what the keys and their corresponding values represent, as well as the types and formats
    # A, C, O, and bestSwitch are all defined in tech report
    A = {}
    C = {}
    O = {}
    bestSwitch = {}

    eventsDict = {}  # Dictionary to keep track of events, children, and scores

    # TODO: rename Minimums minCost
    Minimums = {}  # Dictionary to keep track of minimum reconciliation cost for each (vp, vh)
    oBest = {}  # Dictionary to keep track of the lowest costing events in O
    bestSwitchLocations = {}  # Dictionary to keep track of switch locations
    Score = {}  # Dictionary to calculate the frequency scoring of each event
    # TODO: re-add functions using Score and add option to use it or not
    # TODO: note in the docstring that frequency ('Score') is included in the return value and describe it

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

                # TODO: look into omitting if/elif/else in favor of two if's with leq's instead of < and >
                # If ep2 switching has the lowest cost
                if (C[(ep1, eh)] + bestSwitch[(ep2, eh)]) < (C[(ep2, eh)] +
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

                # If ep1 switching has the lowest cost
                elif (C[(ep2, eh)] + bestSwitch[(ep1, eh)]) < (C[(ep1, eh)] +
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

                # If ep1 switching has the same cost as ep2 switching
                else:

                    # Same processes in this section as noted in previous 'if' blocks,
                    # except done for both cases

                    for location in bestSwitchLocations[(pChild2, vh)]:
                        currentLoc = location[1]
                        if currentLoc is not None:
                            switchList.append(("T", (pChild1, vh),
                                               (pChild2, currentLoc), (Score[(pChild1, vh)] *
                                                                       Score[(pChild2, currentLoc)])))
                        else:
                            switchList.append(("T", (pChild1, vh),
                                               (pChild2, currentLoc), Infinity))
                    for location in bestSwitchLocations[(pChild1, vh)]:
                        currentLoc = location[1]
                        if currentLoc is not None:
                            switchList.append(("T", (pChild2, vh),
                                               (pChild1, currentLoc), (Score[(pChild2, vh)] *
                                                                       Score[(pChild1, currentLoc)])))
                        else:
                            switchList.append(("T", (pChild1, vh),
                                               (pChild2, currentLoc), Infinity))
            else:  # vp is a tip
                SWITCHepeh = Infinity
                switchList = [Infinity]

            # Compute C[(ep, eh)] and add the event or events with that cost
            # to the dictionary eventsDict
            C[(ep, eh)] = min(A[(ep, eh)], DUPepeh, SWITCHepeh)

            # Add the minimum costs for the current edges to the Minimums dict
            Minimums[(vp, vh)] = C[(ep, eh)]

            # TODO: replace all of the min(...) statements with C[(ep, eh)] since they're the same anyhow
            # Find which events produce a minimum and add them to the event dict
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == DUPepeh:
                eventsDict[(vp, vh)].append(dupList)
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == SWITCHepeh:
                eventsDict[(vp, vh)].extend(switchList)
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == A[(ep, eh)]:
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
                    if type(event) is list:

                        # Frequency scores are the last element, and that's what this extracts
                        mapScore += event[-1]

                # Accumulate that and set it as the score for the current key
                Score[key] = mapScore

            # Remove all 'impossible' events from the options
            if Minimums[(vp, vh)] == Infinity:
                del Minimums[(vp, vh)]
                del eventsDict[(vp, vh)]

            # Compute oBest[(vp, vh)], the source of O(ep, eh)
            if vhIsATip: 
                O[(ep, eh)] = C[(ep, eh)]  
                oBest[(vp, vh)] = [(vp, vh)]              
            else:

                # Compute O(ep, eh) if vh is not a tip
                O[(ep, eh)] = min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)])

                # Finds the minimum switch locations for O
                # TODO: comment what oMin and oBest do, more detailed description
                oMin = [i for i, e in enumerate([C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)]])
                        if e == O[(ep, eh)]]
                if 0 in oMin:  # Corresponds to C
                    oBest[(vp, vh)].append((vp, vh))
                # TODO finish comment
                # Corresponds to O for both of the children
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
        eventsDict[key].append(Minimums[key])

    # Use buildEventGraph and findBestRoots to construct the DTL graph dictionary
    # TODO: finish commenting this section
    treeMin = findBestRoots(parasiteTree, Minimums)  # Construct list of min cost mapping nodes involving root of gene tree
    DTL = buildEventGraph(treeMin, eventsDict, {})  # TODO: rename DTLReconGraph

    # The total cost of the best reconciliation
    bestCost = Minimums[treeMin[0]]

    # The total number of MPRs for this optimal cost
    tStart = time.clock()
    nMPRs = countMPRs(True, treeMin, DTL)  # TODO: implement and use the simpler version of countMPRs written in lab with Ran and Eli
    end = time.clock()-tStart
    print "V1 gives {0} recs in {1} secs".format(nMPRs,end)
    tStart = time.clock()
    nMPRs = countMPRsWrapper(treeMin, DTL)
    end = time.clock()-tStart
    print "V2 gives {0} recs in {1} secs".format(nMPRs,end)

    # Returns the graph, the optimal cost, and the number of MPRs
    return DTL, bestCost, nMPRs


def preorderCheck(preOrderList):
    """This takes a list from preorderDTLsort and removes the duplicate tuples"""

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


def LRU(maxsize=None):
    """Memoization decorator that will help count MPRs in the graph
    (the following function). It takes a maxsize for the LRU cache
    and returns a function - the new countMPR function.
    LRU acts as a decorator that will act to greatly improve runtime
    of countMPRs via memoization. It uses the Least Recently Used
    cache model to optimize"""

    # Now define the memoizer, which takes the function we want to memoize
    def memoizer(func):

        # Create the cache data structure to store previous results
        mem_cache = []

        # Now get into the function
        # This function takes the input that countMPRs does, but only
        # really uses 'roots' and 'start'
        def use_memoizer(start, roots, eventGraph):

            # If we're running a brand new call, clear the cache so old results aren't detected from previous calls
            if start:
                del mem_cache[:]

            # Search the cache for our input
            for entry in range(len(mem_cache)):

                # Save the current result for cleaner code
                saved_val = (mem_cache[entry][0], mem_cache[entry][1])

                # If we've found a match in the cache
                if saved_val[0] == roots:

                    # Reorganize the cache to match the LRU model
                    mem_cache.remove(saved_val)
                    mem_cache.insert(0, saved_val)  # Place at front of list since it was recently used

                    # Now return the result
                    return saved_val[1]

            # If it's not in the cache, we'll need to calculate and add it

            # 'None' corresponds to no max size
            if maxsize is not None:

                # Consider whether we need more space
                if len(mem_cache) >= maxsize:

                    # Take off the least recently used element, which would be at the back of the list
                    dummy = mem_cache.pop()

            # Find the value for the given input
            result = func(start, roots, eventGraph)

            # Save it into the cache
            mem_cache.insert(0, (roots, result))

            # Now return the value, which was either calculated or taken from the cache
            return result

        # Return the internal function that implements the memoizer
        return use_memoizer

    # Return the memoizer, which will replace the function that the decorator is used on
    return memoizer


# Utilize the previously defined decorator
@LRU()
def countMPRs(start, roots, eventGraph):  # Don't forget rename DTLReconGraph
    """Takes a boolean value indicating whether the loop is just starting,
    minimum cost roots in a DTL reconciliation graph (output from findBestRoots)
     in a list, and an event graph (output from buildEventGraph).
    Each root should be represented as a tuple (e.g.('a', 'A')).
    This function recursively (with memoization) essentially
    finds the number of unique 'paths' through the solution
    space, and this is used as a proxy for the number of MPRs.
    It does this by starting at the initial 'best roots' (output from
    the findBestRoots function), and goes to those roots. It checks in the
    eventGraph dictionary for the options of a next node for an MPR from those roots.
    It cascades down the graph/tree until reaching a (None, None) map, which is counted
    as the base case. The results then flow back up to the initial roots, and these
    are added to get the total. It returns this number as an integer."""
    # ('L', ('a', 'A'), (None, None))
    # Roots == None is the main base case that works most simply with our algorithm
    if roots == (None, None):  # Signifies either a contemporary event or a single branch
        return 1

    # Initialize the count for the current set of roots
    count = 0

    # Check whether this is when the function was first invoked
    if start:

        # Search through all given starting roots
        for root in roots:

            # Search through all events applicable for that root and add
            # their counts
            # Note we index eventGraph[root][:-1] with -1 because the
            # last entry in an event list is a number, so we want to filter that out
            for event in eventGraph[root][:-1]:

                # Recursively add to the count, this time not calling the function as the first call (start = False)
                count += countMPRs(False, event[1], eventGraph) * countMPRs(False, event[2], eventGraph)

    # Go through the algorithm for the more frequent case where we aren't just starting the function
    else:

        # Loop over all events for the current node
        for event in eventGraph[roots][:-1]:

            # Again, recursively add to the count
            count += countMPRs(False, event[1], eventGraph) * countMPRs(False, event[2], eventGraph)

    return count

def countMPRsWrapper(mappingNodeLst, DTLReconGraph):

    memo = dict()

    count = 0

    for mappingNode in mappingNodeLst:
        count += countMPRs2(mappingNode, DTLReconGraph, memo)

    return count

def countMPRs2(mappingNode, DTLReconGraph, memo):

    if mappingNode in memo:
        return memo[mappingNode]

    if mappingNode == (None, None):
        return 1

    count = 0

    for eventNode in DTLReconGraph[mappingNode][:-1]:
        mappingChild1 = eventNode[1]
        mappingChild2 = eventNode[2]

        count += countMPRs2(mappingChild1, DTLReconGraph, memo) * countMPRs2(mappingChild2, DTLReconGraph, memo)

    memo[mappingNode] = count

    return count


 # TODO: if Parasite is a tree, incorporate that info into the var name
def findBestRoots(Parasite, MinimumDict):
    """Takes Parasite Tree and a dictionary of minimum reconciliation costs
    and returns a list of the minimum cost reconciliation tree roots"""
    treeTops = []
    for key in MinimumDict:
        if key[0] == Parasite['pTop'][1]:
            treeTops.append(key)
    treeMin = []
    min_score = min([MinimumDict[root] for root in treeTops])
    for pair in treeTops:
        if MinimumDict[pair] == min_score:
            treeMin.append(pair)
    return treeMin


# TODO: Rename buildDTLReconGraph, reformat docstring
def buildEventGraph(tupleList, eventDict, uniqueDict):
    """Takes as input: tupleList, a list of minimum reconciliation cost roots,
     eventDict, the dictionary of events and children for each node, and 
     uniqueDict, the dictionary of unique vertex mappings. This returns the 
     completed DTL graph as a dictionary
     :param tupleList:
     :param eventDict:
     :param uniqueDict:
     :return: """
    for vertexPair in tupleList:
        if vertexPair not in uniqueDict:
            uniqueDict[vertexPair] = eventDict[vertexPair]
        for event in eventDict[vertexPair][:-1]:
            for location in event:
                if type(location) is tuple and location != (None, None):
                    buildEventGraph([location], eventDict, uniqueDict)
    return uniqueDict


def reconcile(fileName, D, T, L):
    """Takes as input a newick file, FileName, a duplication cost, a transfer
    cost, and a loss cost. This uses newickFormatReader to extract the host
    tree, parasite tree and tip mapping from the file and then calls DP to
    return the DTL reconciliation graph of the provided newick file"""
    # Note: I have made modifications to the return statement to make Diameter.py possible without re-reconciling.
    host, paras, phi = newickFormatReader.getInput(fileName)
    graph, _, numRecon = DP(host, paras, phi, D, T, L)
    return host, paras, graph, numRecon
