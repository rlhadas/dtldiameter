# DataAnalysis.py
# Written by Eli Zupke, June 2017
# This is a rather messy file for analyzing the data returned by Diameter. It is not intended for release.


import csv
import numpy
import matplotlib.pyplot as plt
import os
import optparse

def displayListValues(list, name):
    print name + ": "
    if not isinstance(list[0],(int, float)):
        print "(Not Number)"
    else:
        print "\tMin:\t{0}".format(min(list))
        print "\tMax:\t{0}".format(max(list))
        print "\tMedian:\t{0}".format(numpy.median(list))
        print "\tMean:\t{0}".format(numpy.mean(list))
        print ""

def read_file(csv_file):

    properties = {}
    column = {}

    with open(csv_file) as file:
        reader = csv.reader(file)
        header = reader.next()

        # Get the index of each property to allow us to search that row.
        for i, element in enumerate(header):
            properties[element] = []
            column[i] = element

        for row in reader:
            for i, property in enumerate(row):
                properties[column[i]] += [property]
    length = len(properties[column[0]])
    return properties, length



def make_plot(file, zero_loss, non_normalized, timings, gene_count_list, diameter_list, diameter_over_gene_list, mpr_list, DP_timings,
         diameter_timings, total_timings):
    size = 4
    color = 'black'
    if zero_loss:
        diameter_ylim_b = -0.05
        diameter_ylim_t = 1.05
    else:
        diameter_ylim_b = -0.1
        diameter_ylim_t = 2.1

    gene_xlim = max(gene_count_list) * 1.05

    name = ""
    if zero_loss:
        name = " (Zero Loss)"


    if non_normalized:

        fig, ax = plt.subplots(ncols=3, nrows=1)
        fig.canvas.set_window_title("{0} Plots{1}".format(file, name))

        diameter = ax[1]
        diameter_hist = ax[0]
        mpr_diameter = ax[2]
        ax[0].set_ylabel("Diameter")
        diameter.scatter(gene_count_list, diameter_list, c=color, s=size)
        diameter.set_xlabel("Gene Tree Size")
        diameter.set_xlim(0, gene_xlim)
        diameter.set_title("Diameter vs. Gene Count")
        diameter.grid()

        diameter_hist.hist(diameter_list, 100, orientation='horizontal')
        # diameter_hist.set_ylabel("Diameter")
        diameter_hist.set_xlabel("Number of Gene Families")
        diameter_hist.set_title("Diameter")
        diameter_hist.grid()

        mpr_diameter.scatter(mpr_list, diameter_list, c=color, s=size)
        mpr_diameter.set_xlabel("MPR Count")
        mpr_diameter.set_title("Diameter vs. MPR Count")
        mpr_diameter.grid()
        mpr_diameter.set_xscale('log')


    fig, ax = plt.subplots(ncols=3, nrows=1)
    fig.canvas.set_window_title("{0} Normalized Plots{1}".format(file, name))
    norm_diameter_hist = ax[0]
    norm_mpr_diameter = ax[2]
    norm_diameter = ax[1]
    ax[0].set_ylabel("Normalized Diameter")

    norm_diameter.scatter(gene_count_list, diameter_over_gene_list, c=color, s=size)
    norm_diameter.set_xlabel("Gene Tree Size")
    norm_diameter.set_xlim(0, gene_xlim)
    # norm_diameter.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter.set_ylim(diameter_ylim_b, diameter_ylim_t)
    norm_diameter.set_title("Normalized Diameter vs. Gene Tree Size")
    norm_diameter.grid()

    norm_diameter_hist.grid()
    norm_diameter_hist.hist(diameter_over_gene_list, 100, orientation='horizontal')
    # norm_diameter_hist.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter_hist.set_xlabel("Number of Gene Families")
    norm_diameter_hist.set_title("Normalized Diameter Counts")
    norm_diameter_hist.set_ylim(diameter_ylim_b, diameter_ylim_t)

    norm_mpr_diameter.scatter(mpr_list, diameter_over_gene_list, c=color, s=size)
    norm_mpr_diameter.set_ylim(diameter_ylim_b, diameter_ylim_t)
    norm_mpr_diameter.set_xlabel("MPR Count")
    norm_mpr_diameter.set_title("Normalized Diameter vs. MPR Count")
    norm_mpr_diameter.grid()
    norm_mpr_diameter.set_xscale('log')

    # plt.show()
    # return
    if timings:
        fig, ax = plt.subplots(ncols=3, nrows=1)
        fig.canvas.set_window_title("{0} Running Times{1}".format(file, name))
        DP_time = ax[0]
        diameter_time = ax[1]
        total_time = ax[2]
        DP_time.scatter(gene_count_list, DP_timings, c=color, s=size)
        DP_time.set_xlabel("Gene Tree Size")
        DP_time.set_ylabel("Time (seconds)")
        DP_time.set_title("Computing Reconciliation Graph")
        DP_time.grid()
        DP_time.set_yscale('log')
        DP_time.set_xscale('log')
        diameter_time.scatter(gene_count_list, diameter_timings, c=color, s=size)
        diameter_time.set_xlabel("Gene Tree Size")
        diameter_time.set_ylabel("Time (seconds)")
        diameter_time.set_title("Computing Diameter")
        diameter_time.grid()
        diameter_time.set_yscale('log')
        diameter_time.set_ylim(0.001,10**4)
        diameter_time.set_xscale('log')
        total_time.scatter(gene_count_list, total_timings, c=color, s=size)
        total_time.set_xlabel("Gene Tree Size")
        total_time.set_ylabel("Time (seconds)")
        total_time.set_title("Total Running Time")
        total_time.grid()
        total_time.set_yscale('log')
        total_time.set_xscale('log')

def findExtrema(csv_file, properties, non_normalized, timings, plot, latex):
    """Finds the minimums, maximums, medians, and means of the:
        MPR Count
        Diameter
        Gene Count
        MPR Count/Diameter
        Diameter/Gene Count
        MPR Count/(Diameter/Gene Count)
    Of a csv file created by Diameter.py"""

    if plot:
        plt.rc('text', usetex=latex)
        if latex:
            plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
            plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

    mpr_list = []
    DP_timings = []
    number_list = []
    number = 0
    DTL = "n/a"

    diameter_present = False
    zero_loss_present = False

    file_props, length = read_file(csv_file)
    DTL = file_props["Costs"][0]
    mpr_list = file_props["MPR Count"]
    mpr_list = map(lambda e: float(e), mpr_list)
    gene_count_list = file_props["Gene Node Count"]
    gene_count_list= map(lambda e: float(e), gene_count_list)
    species_count_list = file_props["Species Node Count"]
    species_count_list = map(lambda e: float(e), species_count_list)
    DP_timings = file_props["DTLReconGraph Computation Time"]
    DP_timings = map(lambda e: float(e), DP_timings)
    prop_list_dict = {}
    timing_list_dict = {}

    for property in properties:
        prop_list_dict[property] = file_props[property]
        if "Diameter" in property:
            prop_list_dict[property] = map(lambda e: float(e), prop_list_dict[property])
            timing_list_dict[property + " Computation Time"] = file_props[property + " Computation Time"]

    total_timings = DP_timings[:]

    for timing_list in timing_list_dict:
        timing_list_dict[timing_list] = map(lambda e: float(e), timing_list_dict[timing_list])
        cur_list = timing_list_dict[timing_list]
        for i, time in enumerate(cur_list):
            total_timings[i] += time


    displayListValues(mpr_list, "MPR Count")
    displayListValues(gene_count_list, "Gene Tree Size")
    displayListValues(gene_count_list, "Species Tree Size")
    #displayListValues(mpr_over_d_list, "MPR Count/Diameter")
    #displayListValues(mpr_over_normalized_d_list, "MPR Count/(Diameter/Gene Count)")
    if timings:
        displayListValues(DP_timings, "Reconciliation Running Time (seconds)")
        displayListValues(total_timings, "Total Running Time (seconds)")

    for property_list in prop_list_dict:
        displayListValues(prop_list_dict[property_list], property_list)



    filepath, extension = os.path.splitext(csv_file)
    zero_loss_file = filepath + "_zl" + extension
    if plot:
        if "Diameter" in prop_list_dict:

            normalized = map(lambda i: prop_list_dict["Diameter"][i[0]] / gene_count_list[i[0]], enumerate(prop_list_dict["Diameter"]))
            make_plot(csv_file, False, non_normalized, timings, gene_count_list, prop_list_dict["Diameter"], normalized,
             mpr_list, DP_timings, timing_list_dict["Diameter Computation Time"], total_timings)
        if "Zero Loss Diameter" in prop_list_dict:
            normalized = map(lambda i: prop_list_dict["Zero Loss Diameter"][i[0]] / gene_count_list[i[0]],
                             enumerate(prop_list_dict["Zero Loss Diameter"]))
            make_plot(csv_file, True, non_normalized, timings, gene_count_list, prop_list_dict["Zero Loss Diameter"], normalized,
             mpr_list, DP_timings, timing_list_dict["Zero Loss Diameter Computation Time"], total_timings)


def findSpecific(col, value, csv_file="COG_Pilot_Log_02.csv", tol=0):
    with open(csv_file) as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] != "File Name":
                if col == -1:
                    if value-tol <= (float(row[2])/float(row[3])) <= value+tol:
                        print row
                elif value-tol <= float(row[col]) <= value+tol:
                    print row


def check_files(log, path):
    if path[-1] != "/":
        path = path + "/"

    data_files = map(lambda x: path + x, os.listdir(path))
    duplicates = []
    with open(log) as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] not in data_files:
                duplicates += [row[0]]
            else:
                data_files.remove(row[0])
        print "Duplicates ({0}):\t {1}".format(len(duplicates[1:]), duplicates[1:])
        print "Missing ({0}):\t {1}".format(len(data_files), data_files)


def compare_logs(log1="New_COG_01.csv", log2="New_COG_02_zl.csv"):
    """"""
    log1_diams = []
    log2_diams = []
    filenames = []
    mismatches = 0
    difference = 0.0
    count = 0
    with open(log1) as file:
        reader = csv.reader(file)
        for row in reader:
            filenames += [row[0]]
            log1_diams += [row[3]]
    with open(log2) as file:
        reader = csv.reader(file)
        for row in reader:
            log2_diams += [row[3]]
    for i in range(0, len(log1_diams)):
        count += 1
        if log1_diams[i] != log2_diams[i]:
            difference += (float(log1_diams[i]) - int(log2_diams[i]))/int(log2_diams[i])
            print "Mismatch in {0}: {1} vs. {2}".format(filenames[i], log1_diams[i], log2_diams[i])
            mismatches += 1
        #else:
            #print "Match in {0}: {1} vs. {2}".format(filenames[i], log1_diams[i], log2_diams[i])
    print "{0} mismatches, or {0}/{2} = {1}%".format(mismatches,mismatches/(float(count))*100,count)


def compare_log_ratios(log1="New_COG_02.csv", log2="COG_Median_med.csv", output_log="Compare.csv"):
    """"""
    log1_diams = []
    log2_diams = []
    gene_node_count = []
    mprs = []
    filenames = []
    count = 0
    with open(log1) as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if i != 0:
                filenames += [row[0]]
                log1_diams += [row[3]]
                gene_node_count += [int(row[4])]
                mprs += [int(row[2])]
    with open(log2) as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if i != 0:
                log2_diams += [row[3]]
    ratios = []
    nonzero_gene = []
    nonzero_mprs = []
    with open(output_log, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(["Filename","Regular","Median","Ratio"])
        for i in range(0, len(log1_diams)):
            if int(log1_diams[i]) * int(log2_diams[i]) is not 0:
                ratios += [float(log2_diams[i])/float(log1_diams[i])]
                nonzero_gene += [gene_node_count[i]]
                nonzero_mprs += [mprs[i]]
                writer.writerow([filenames[i], log1_diams[i], log2_diams[i], ratios[-1]])
    displayListValues(ratios, "Ratio (Median Diameter Over Diameter)")
    size = 4
    color = 'black'

    gene_xlim = max(nonzero_gene) * 1.05

    name = ""
    fig, ax = plt.subplots(ncols=3, nrows=1)
    fig.canvas.set_window_title("Ratios")
    norm_diameter_hist = ax[0]
    norm_mpr_diameter = ax[2]
    norm_diameter = ax[1]
    ax[0].set_ylabel("Ratio (Median Diameter Over Diameter)")

    ratio_ylim_b = 0.5
    ratio_ylim_t = 1

    norm_diameter.scatter(nonzero_gene, ratios, c=color, s=size)
    norm_diameter.set_xlabel("Gene Tree Size")
    norm_diameter.set_xlim(0, gene_xlim)
    # norm_diameter.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter.set_ylim(ratio_ylim_b, ratio_ylim_t)
    norm_diameter.set_title("Ratio vs. Gene Tree Size")
    norm_diameter.grid()

    norm_diameter_hist.grid()
    norm_diameter_hist.hist(ratios, 100, orientation='horizontal')
    # norm_diameter_hist.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter_hist.set_xlabel("Number of Gene Families")
    norm_diameter_hist.set_title("Ratio Count")
    norm_diameter_hist.set_ylim(ratio_ylim_b, ratio_ylim_t)

    norm_mpr_diameter.scatter(nonzero_mprs, ratios, c=color, s=size)
    norm_mpr_diameter.set_ylim(ratio_ylim_b, ratio_ylim_t)
    norm_mpr_diameter.set_xlabel("MPR Count")
    norm_mpr_diameter.set_title("Ratio vs. MPR Count")
    norm_mpr_diameter.grid()
    norm_mpr_diameter.set_xscale('log')

    plt.show()



def main():
    """Processes command line arguments"""
    usage = "usage: %prog [options] file"
    p = optparse.OptionParser(usage=usage)
    p.add_option("-p", "--plot", dest="plot", action="store_true", default=False,
                 help="outputs some plots!")
    p.add_option("-z", "--zero-loss", dest="zero_loss", action="store_true", default=False ,
                 help="also plot related zero-loss logfiles")
    p.add_option("-l", "--use-latex", dest="use_latex", action="store_true", default=False,
                 help="use LaTeX for plot text rendering (you must have LaTeX installed on your system!)")
    p.add_option("-t", "--timings", dest="timings", action="store_true", default=False,
                 help="includes algorithm timing plot")
    p.add_option("-n", "--non-normalized", dest="non_normalized", action="store_true", default=False,
                 help="includes plot with non-normalized diameter")
    p.add_option("-c", "--compare", dest="compare_file", help="compare the diameters of this logfile to another, "
                                                              "and report any mismatches between the two",
                 metavar="COMPARE_FILE")
    p.add_option("-k", "--check-files", dest="check_path", help="Compare the files in the logfile with the files in the"
                                                                "directory, and report any duplicate files in the log,"
                                                                "and any files in the directory but not in the log.",
                 metavar="CHECK_PATH")


    (options, args) = p.parse_args()
    if len(args) != 1:
        p.error("1 argument must be provided: file")
    file = args[0]
    zero_loss = options.zero_loss
    latex = options.use_latex
    plot = options.plot
    compare_file = options.compare_file
    check = options.check_path
    timings = options.timings
    non_normalized = options.non_normalized
    if latex and not plot:
        print "Warning: option '-l' (--use-latex) has no effect without option '-p' (--plot)!"
    if not os.path.isfile(file):
        p.error("File not found, '{0}'. Please be sure you typed the name correctly!".format(file))
    else:
        findExtrema(file, ["Diameter", "Zero Loss Diameter"] if zero_loss else ["Diameter"], non_normalized, timings, plot, latex)
        if compare_file is not None:
            compare_logs(file, compare_file)
        if check is not None:
            check_files(file, check)
        if plot:
            plt.show()


if __name__ == "__main__":
    main()