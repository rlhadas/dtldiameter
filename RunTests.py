import NewDiameter
import csv
import time
import DTLReconGraph
import DTLMedian
import os

# Used for command line arguments:
import sys
import re
import traceback
import optparse


def write_to_csv(csv_file, costs, filename, mpr_count, gene_node_count, species_node_count,
                 DTLReconGraph_time_taken, properties):
    """Takes a large amount of information about a diameter solution and appends it as one row to the provided csv file.
    :param csv_file:                    The csv file to write to
    :param costs:                       A string representing the costs used to calculate the DTL recon graph
    :param filename:                    The name of the file that was reconciled
    :param mpr_count:                   The number of MPRS found
    :param gene_node_count:             The total number of nodes in the gene tree
    :param species_node_count:          The total number of nodes in the species tree
    :param DTLReconGraph_time_taken:    The amount of time that DTLReconGraph took to run
    :param properties:                  A list of tuples, where each tuple represents one property of the graph that was
                                         computed (for example, diameter), and has the following format:
                                         0: Name of the property
                                         1: Property computed
                                         2: Time taken
    """
    file_exists = os.path.isfile(csv_file)
    # If they only supply one argument, let's correct it for them.
    if isinstance(properties, tuple):
        properties = [properties]

    with open(csv_file, 'a') as output_file:
        writer = csv.writer(output_file)

        # Write the headers if we need to.
        if not file_exists:
            header = ["File Name", "Date Completed", "Costs", "MPR Count", "Gene Node Count", "Species Node Count",
                      "DTLReconGraph Computation Time"]
            for property in properties:
                header += [property[0], "{0} Computation Time".format(property[0])]
            writer.writerow(header)
        numbers = [filename, time.strftime("%c"), costs, mpr_count, gene_node_count, species_node_count,
                         DTLReconGraph_time_taken]
        for property in properties:
            numbers += [property[1], property[2]]
        writer.writerow(numbers)



def calculate_diameter_from_file(filename, D, T, L, log=None, debug=False, verbose=True, zero_loss=False, median=False,
                                 worst_median=False, median_cluster=0):

    """This function computes the diameter of space of MPRs in a DTL reconciliation problem,
     as measured by the symmetric set distance between the events of the two reconciliations of the pair
      that has the highest such difference.

      :param filename:      The path to the newick tree file to reconcile
      :param D:             The cost for duplication events
      :param T:             The cost for transfer events
      :param L:             The cost for loss events
      :param log:           The csv file to output results to (will create it if it does not exist)
      :param debug:         Whether to print out all of the tables
      :return:              Nothing, but we output results to a csv file."""

    # These statements check to make sure that all arguments were entered correctly.
    assert isinstance(log, (str, unicode)) or log is None
    assert isinstance(filename, (str, unicode))
    assert isinstance(D, (int, float))
    assert isinstance(T, (int, float))
    assert isinstance(L, (int, float))
    assert isinstance(debug, bool)

    # Record the time that DTLReconGraph starts
    start_time = time.clock()

    # Get everything we need from DTLReconGraph
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile(filename, D, T, L)

    # Record the time that this code starts

    # The gene tree needs to be in node format, not edge format, so we find that now.
    # (This also puts the gene_tree into postorder, as an ordered dict)
    gene_tree, gene_tree_root, gene_node_count = NewDiameter.reformat_tree(edge_gene_tree, "pTop")

    species_tree, species_tree_root, species_node_count = NewDiameter.reformat_tree(edge_species_tree, "hTop")

    # The DTL reconciliation graph as provided by DTLReconGraph has some extraneous numbers. We remove those here.
    NewDiameter.clean_graph(dtl_recon_graph, gene_tree_root)

    # And record the amount of time DTLReconGraph + cleaning up the graph took
    DTLReconGraph_time_taken = time.clock() - start_time

    if verbose:
        print "Reconciliation Graph Made in \033[33m\033[1m{0} seconds\033[0m".format(DTLReconGraph_time_taken)

    results = []

    start_time = time.clock()

    # Now we draw the rest of the owl
    diameter = NewDiameter.new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph, debug, False)

    # And record how long it took to compute the diameter.
    diameter_time_taken = time.clock() - start_time
    results += [("Diameter", diameter, diameter_time_taken)]

    if zero_loss:
        start_time = time.clock()
        zl_diameter = NewDiameter.new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph, debug, True)
        zl_diameter_time_taken = time.clock()-start_time
        results += [("Zero Loss Diameter", zl_diameter, zl_diameter_time_taken)]

    median_reconciliation = {}

    if median:
        start_time = time.clock()
        preorder_mapping_node_list = DTLMedian.mapping_node_sort(gene_tree, species_tree, dtl_recon_graph.keys())

        # Find the dictionary for frequency scores for the given mapping nodes and graph, as well as the given gene root
        scoresDict = DTLMedian.generate_scores(list(reversed(preorder_mapping_node_list)), dtl_recon_graph, gene_tree_root)

        median_reconciliation, n_meds, _ = DTLMedian.compute_median(dtl_recon_graph, scoresDict[0], preorder_mapping_node_list,
                                                                    best_roots)
        median_time_taken = time.clock()-start_time
        results += [("Median Count", n_meds, median_time_taken)]
    if median and worst_median:
        start_time = time.clock()
        worst_median_diameter = NewDiameter.new_diameter_algorithm(species_tree, gene_tree, gene_tree_root,
                                                                   median_reconciliation, dtl_recon_graph, debug, False)
        worst_median_diameter_time_taken = time.clock()-start_time
        results += [("Worst Median Diameter", worst_median_diameter, worst_median_diameter_time_taken)]

    if median_cluster > 0:

        start_time = time.clock()
        avg = 0.0
        _, file_log_name = os.path.split(filename)
        file_log_name, _ = os.path.splitext(file_log_name)
        file_log_path = os.path.splitext(log)[0] + "/" + file_log_name + ".csv"
        costs = "D: {0} T: {1} L: {2}".format(D, T, L)

        # Every time this loop repeats, we calculate another random median and find its diameter
        for i in range(0, median_cluster):
            start_sub_time = time.clock()

            random_median = {}

            median_hash = hash(frozenset(random_median.items()))

            end_random_time = time.clock() - start_sub_time
            start_sub_time = time.clock()

            random_median_diameter = NewDiameter.new_diameter_algorithm(species_tree, gene_tree, gene_tree_root,
                                                                      random_median, dtl_recon_graph, debug, False)

            end_sub_time = time.clock() - start_sub_time
            sub_results = [("Random Median", median_hash, end_random_time),
                           ("Random Median Diameter", random_median_diameter, end_sub_time)]
            avg += random_median_diameter

            if log is not None:
                # TODO make filename not include folder or extension for here
                write_to_csv(file_log_path , costs, filename, mpr_count, gene_node_count, species_node_count,
                             DTLReconGraph_time_taken, sub_results)

        avg /= median_cluster
        random_median_diameter_time_taken = time.clock()-start_time
        results += [("Best Median Diameter", avg, random_median_diameter_time_taken)]


    if median and worst_median:
        start_time = time.clock()
        worst_median_diameter = NewDiameter.new_diameter_algorithm(species_tree, gene_tree, gene_tree_root, median_reconciliation,
                                                         dtl_recon_graph, debug, False)
        worst_median_diameter_time_taken = time.clock() - start_time
        results += [("Worst Median Diameter", worst_median_diameter, worst_median_diameter_time_taken)]

    #if verbose:
     #   print "The diameter of the given reconciliation graph is \033[33m\033[1m{0}\033[0m, (or \033[33m\033[1m{1}\033[0m if losses do not affect the diameter)".format(diameter, zl_diameter)

    # Timing data is inaccurate in debug mode (print statements take too long), so we only give it to the user in non-
    # debug mode.
    if not debug and verbose:
        print "Diameter found in \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken)
        print "Total time: \033[33m\033[1m{0} seconds\033[0m".format(diameter_time_taken + DTLReconGraph_time_taken)

    # Now, we write our results to a csv file.
    if log is not None:
        costs = "D: {0} T: {1} L: {2}".format(D, T, L)
        write_to_csv(log + ".csv", costs, filename, mpr_count, gene_node_count, species_node_count, DTLReconGraph_time_taken,
                     results)

    # And we're done.
    return


def repeatedly_calculate_diameter(file_pattern, start, end, d, t, l, log=None, debug=False, verbose=True, loud=False, cluster=0):
    """Iterates over a lot of input files and finds the diameter of all of them.
    :param file_pattern: A string contains the name of the files to be used, with the counting number replaced with #'s
    :param start:       Numbered file to start on
    :param end:         Numbered file to end with
    :param d:           Duplication event cost
    :param t:           Transfer event cost
    :param l:           Loss event cost
    :param log:         csv file to log results to
    :param debug:       Whether to print out every DP table made (not recommended)
    :return:
    """
    match = re.match("([^#]*)(#+)([^#]*)", file_pattern)
    if not match:
        print "Filepath '" + file_pattern + "' not understood. Please enter the path to your files, with repeated hash marks" \
                                       "(#'s) in place of sequential numbering."
        return
    fill = len(match.group(2))
    if fill < len(str(end-1)) or fill < len(str(start)):
        print "Starting or ending number is larger than '{1}' supports ({0})!".format((10 ** fill) -1, file_pattern)
        return
    print "Running {4} sequential jobs on files '{3}' with costs D = {0}, T = {1}, and L = {2}".format(d, t, l, file_pattern, end - start)
    for i in range(start, end):
        cur_file = "{0}{1}{2}".format(match.group(1), str(i).zfill(fill), match.group(3))
        print "Reconciling {0}".format(cur_file)
        try:
            calculate_diameter_from_file(cur_file, d, t, l, log, debug, verbose, median_cluster=cluster)
        except IOError:
            print "(File Not Found)"
        except (KeyboardInterrupt, SystemExit):
            raise  # Don't prevent the user from exiting the program.
        except:
            if loud:
                print "\07"
            if verbose:
                print traceback.print_exc(sys.exc_traceback)
            print "Could not reconcile file '{0}'. Continuing, but please make sure the file was formatted correctly!"\
                .format(cur_file)



def main():
    """Processes command line arguments"""
    usage = "usage: %prog [options] file d t l"
    p = optparse.OptionParser(usage=usage)
    p.add_option("-l", "--log", dest="logfile", help="writes a logfile in CSV format to LOGFILE", metavar="LOGFILE")
    p.add_option("-i", "--iterate", dest="count", action="store", nargs=2, help="calculates every file matching a "
                                                                                "pattern (defined in the file argument)"
                                                                                " with number MIN to MAX.",
                 metavar="MIN MAX", type=int)
    p.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help= "prints out every DP table with size less"
                                                                           "than 30x30")
    p.add_option("-q", "--quiet", dest="verbose", action="store_false", default=True,
                 help="suppresses (most) text output")
    p.add_option("-c", "--cluster", dest="cluster", action="store", default=0,
                 help="suppresses (most) text output")
    p.add_option("-L", "--loud", dest="loud", action="store_true", default=False,
                 help="print the bell character after each failed file")

    (options, args) = p.parse_args()
    if len(args) != 4:
        p.error("4 arguments must be provided: file, d, t, and l")
    file = args[0]
    d = float(args[1])
    t = float(args[2])
    l = float(args[3])
    log = options.logfile
    debug = options.debug
    verbose = options.verbose
    loud = options.loud
    cluster = options.cluster
    if not (log or debug or verbose):
        p.error("some form of output must be specified! (-l or -d must be used when -q is used)")
    elif options.count is not None:
        rep = options.count
        repeatedly_calculate_diameter(file, rep[0], rep[1], d, t, l, log, debug, verbose, loud, cluster)
    else:
        calculate_diameter_from_file(file, d, t, l, log, debug, verbose, median_cluster=cluster)

if __name__ == "__main__":
    main()
