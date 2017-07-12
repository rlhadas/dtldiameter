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

def read_file(csv_file, mpr_strip=0, mpr_equals_median_strip=False):

    properties = {}
    column = {}
    ignore_mprs = []
    if mpr_strip > 0:
        ignore_mprs = map(str, range(0, mpr_strip))

    with open(csv_file) as file:
        reader = csv.reader(file)
        header = reader.next()

        MPR_count_column = 0
        median_count_column = 0

        # Get the index of each property to allow us to search that row.
        for i, element in enumerate(header):
            properties[element] = []
            column[i] = element
            if element == "MPR Count":
                MPR_count_column = i
            elif element == "Median Count":
                median_count_column = i

        for row in reader:

            if row[MPR_count_column] not in ignore_mprs and (not mpr_equals_median_strip or row[MPR_count_column] !=
                                                             row[median_count_column]):
                for i, property in enumerate(row):
                    properties[column[i]] += [property]

    length = len(properties[column[0]])
    name_to_row = {v: k for k, v in column.iteritems()}
    return properties, length, name_to_row

def set_label(axis, label, letter, latex):
    if latex:
        axis.set_xlabel(label+"\n"r"{\fontsize{30pt}{3em}\selectfont{}("+letter+")}", linespacing=2.5,
                                 labelpad=20)
    else:
        axis.set_xlabel(label)

def make_plot(file, zero_loss, non_normalized, timings, gene_count_list, diameter_list, diameter_name,
              normalized_diameter, normalized_diameter_name, mpr_list, DP_timings,
              diameter_timings, total_timings, name, latex):
    size = 4
    color = 'black'
    if zero_loss:
        diameter_ylim_b = -0.05
        diameter_ylim_t = 1.05
    else:
        diameter_ylim_b = -0.1
        diameter_ylim_t = 2.1

    gene_xlim = max(gene_count_list) * 1.05

    if non_normalized:
        y_max = max(diameter_list)
        y_min = min(diameter_list)
        padding = (y_max - y_min) * 0.05
        y_bottom = y_min - padding
        y_top = y_max + padding

        fig, ax = plt.subplots(ncols=3, nrows=1)
        fig.canvas.set_window_title("{0} Plots {1}".format(file, name))

        diameter = ax[1]
        diameter_hist = ax[0]
        mpr_diameter = ax[2]
        ax[0].set_ylabel(diameter_name)
        diameter.scatter(gene_count_list, diameter_list, c=color, s=size)
        set_label(diameter, "Gene Tree Size", "b", latex)
        diameter.set_yscale('log')
        diameter.set_xlim(0, gene_xlim)
        diameter.set_title("{0} vs. Gene Count".format(diameter_name))
        diameter.grid()

        diameter_hist.hist(diameter_list, orientation='horizontal', bins=numpy.logspace(numpy.log10(y_bottom), 35, 100))
        # diameter_hist.set_ylabel("Diameter")
        set_label(diameter_hist, "Number of Gene Families", "a", latex)
        diameter_hist.set_title(diameter_name)
        diameter_hist.set_yscale('log')
        diameter_hist.grid()

        mpr_diameter.scatter(mpr_list, diameter_list, c=color, s=size)
        set_label(mpr_diameter, "MPR Count", "c", latex)
        mpr_diameter.set_title("{0} vs. MPR Count".format(diameter_name))
        mpr_diameter.grid()
        mpr_diameter.set_xscale('log')
        mpr_diameter.set_yscale('log')


    fig, ax = plt.subplots(ncols=3, nrows=1)
    fig.canvas.set_window_title("{0} Normalized Plots {1}".format(file, name))
    norm_diameter_hist = ax[0]
    norm_mpr_diameter = ax[2]
    norm_diameter = ax[1]
    ax[0].set_ylabel(normalized_diameter_name)

    norm_diameter.scatter(gene_count_list, normalized_diameter, c=color, s=size)
    set_label(norm_diameter, "Gene Tree Size", "b", latex)
    norm_diameter.set_xlim(0, gene_xlim)
    # norm_diameter.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter.set_ylim(diameter_ylim_b, diameter_ylim_t)
    norm_diameter.set_title("{0} vs. Gene Tree Size".format(normalized_diameter_name))
    norm_diameter.grid()

    norm_diameter_hist.grid()
    norm_diameter_hist.hist(normalized_diameter, 100, orientation='horizontal')
    # norm_diameter_hist.set_ylabel("Diameter (normalized to gene node count)")
    set_label(norm_diameter_hist, "Number of Gene Families", "a", latex)
    norm_diameter_hist.set_title("{0} Counts".format(normalized_diameter_name))
    norm_diameter_hist.set_ylim(diameter_ylim_b, diameter_ylim_t)

    norm_mpr_diameter.scatter(mpr_list, normalized_diameter, c=color, s=size)
    norm_mpr_diameter.set_ylim(diameter_ylim_b, diameter_ylim_t)
    set_label(norm_mpr_diameter, "MPR Count", "c", latex)
    norm_mpr_diameter.set_title("{0} vs. MPR Count".format(normalized_diameter_name))
    norm_mpr_diameter.grid()
    norm_mpr_diameter.set_xscale('log')

    # plt.show()
    # return
    if timings:
        fig, ax = plt.subplots(ncols=3, nrows=1)
        fig.canvas.set_window_title("{0} Running Times {1}".format(file, name))
        DP_time = ax[0]
        diameter_time = ax[1]
        total_time = ax[2]
        DP_time.scatter(gene_count_list, DP_timings, c=color, s=size)
        set_label(DP_time, "Gene Tree Size", "a", latex)
        DP_time.set_ylabel("Time (seconds)")
        DP_time.set_title("Computing Reconciliation Graph")
        DP_time.grid()
        DP_time.set_yscale('log')
        DP_time.set_xscale('log')
        diameter_time.scatter(gene_count_list, diameter_timings, c=color, s=size)
        set_label(diameter_time, "Gene Tree Size", "b", latex)
        diameter_time.set_ylabel("Time (seconds)")
        diameter_time.set_title("Computing {0}".format(diameter_name))
        diameter_time.grid()
        diameter_time.set_yscale('log')
        diameter_time.set_ylim(0.001,10**4)
        diameter_time.set_xscale('log')
        total_time.scatter(gene_count_list, total_timings, c=color, s=size)
        set_label(total_time, "Gene Tree Size", "c", latex)
        total_time.set_ylabel("Time (seconds)")
        total_time.set_title("Total Running Time")
        total_time.grid()
        total_time.set_yscale('log')
        total_time.set_xscale('log')

def findExtrema(csv_file, properties, non_normalized, timings, plot, latex, strip_mprs, strip_equal):
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

    diameter_present = False
    zero_loss_present = False

    file_props, length, _ = read_file(csv_file, strip_mprs, strip_equal)
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
        if "Diameter" in property or "Count" in property or "Number" in property:
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

    name = ""
    if strip_mprs == 2:
        name += " >1 MPR"
    elif strip_mprs > 2:
        name += " >{0} MPRs".format(strip_mprs-1)

    if strip_equal:
        name += " MPR != Median"

    filepath, extension = os.path.splitext(csv_file)
    if plot:
        if "Diameter" in prop_list_dict:

            normalized = map(lambda i: prop_list_dict["Diameter"][i[0]] / gene_count_list[i[0]], enumerate(prop_list_dict["Diameter"]))
            make_plot(csv_file, False, non_normalized, timings, gene_count_list, prop_list_dict["Diameter"], "Diameter", normalized,
            "Normalized Diameter", mpr_list, DP_timings, timing_list_dict["Diameter Computation Time"], total_timings, "Diameter" + name, latex)
        if "Zero Loss Diameter" in prop_list_dict:
            normalized = map(lambda i: prop_list_dict["Zero Loss Diameter"][i[0]] / gene_count_list[i[0]],
                             enumerate(prop_list_dict["Zero Loss Diameter"]))
            make_plot(csv_file, True, non_normalized, timings, gene_count_list, prop_list_dict["Zero Loss Diameter"], "Zero Loss Diameter", normalized,
             "Normalized Zero Loss Diameter", mpr_list, DP_timings, timing_list_dict["Zero Loss Diameter Computation Time"], total_timings, "Zero Loss Diameter " +name, latex)
        if "Median Count" in prop_list_dict:
            normalized = map(lambda i: prop_list_dict["Median Count"][i[0]] / mpr_list[i[0]],
                             enumerate(prop_list_dict["Median Count"]))
            make_plot(csv_file, True, non_normalized, timings, gene_count_list, prop_list_dict["Median Count"], "Median Count", normalized,
            "Normalized Median Count", mpr_list, DP_timings, timing_list_dict["Median Count Computation Time"], total_timings, "Median Count " +name, latex)
        if "Worst Median Diameter" in prop_list_dict:
            normalized = map(lambda i: prop_list_dict["Worst Median Diameter"][i[0]] / gene_count_list[i[0]],
                             enumerate(prop_list_dict["Worst Median Diameter"]))
            make_plot(csv_file, False, non_normalized, timings, gene_count_list, prop_list_dict["Worst Median Diameter"],
                      "Worst Median Diameter", normalized,
                      "Normalized Worst Median Diameter", mpr_list, DP_timings,
                      timing_list_dict["Worst Median Diameter Computation Time"], total_timings, "Worst Median Diameter " +name, latex)


def find_specific(csv_file="COG_Median_13.csv"):
    file_props, length, column_lookup = read_file(csv_file, False, False)

    mprs_1 = 0
    mprs_2 = 0
    mprs_gt_2 = 0
    mprs_gt_2_not = 0
    mprs_gt_2_eq = 0
    mpr_eq_med = 0
    mpr_not_med = 0
    mpr_lt = [0]*13
    mpr_gt = 0

    for i in range(0, length):
        mpr_count = file_props["MPR Count"][i]
        median = file_props["Median Count"][i]
        if mpr_count == "1":
            mprs_1 += 1
        elif mpr_count == "2":
            mprs_2 += 1
        else:
            mprs_gt_2 += 1
            if mpr_count == median:
                mprs_gt_2_not += 1
            else:
                mprs_gt_2_eq += 1
        if mpr_count == median:
            mpr_eq_med += 1
        else:
            mpr_not_med += 1
        med_count = int(median)
        placed = False
        for i in range(0, len(mpr_lt)):
            lower = 2**i
            upper = 2**(i+1)
            if lower <= med_count < upper:
                mpr_lt[i] += 1
                placed = True
        if not placed:
            mpr_gt += 1

    print "{0} families observed.".format(length)
    print "MPR count == 1: {0}".format(mprs_1)
    print "MPR count == 2: {0}".format(mprs_2)
    print "MPR count > 2: {0}".format(mprs_gt_2)
    print "MPR == median count: {0}".format(mpr_eq_med)
    print "MPR != median count: {0}".format(mpr_not_med)
    print "MPR count > 2 & MPR == median count: {0}".format(mprs_gt_2_not)
    print "MPR count > 2 & MPR != median count: {0}".format(mprs_gt_2_eq)
    print ""
    print "Families per median range:"
    for i, count in enumerate(mpr_lt):
        lower = 2 ** i
        upper = 2 ** (i+1)
        print "[{0}, {1}): \t\t{2}".format(lower, upper, count)
    print ">{0}:\t\t{1}".format(2**len(mpr_lt), mpr_gt)


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
                 help="plot the diameter when losses are not included")
    p.add_option("-D", "--no-diameter", dest="diameter", action="store_false", default=True,
                 help="don't plot the diameter")
    p.add_option("-l", "--use-latex", dest="use_latex", action="store_true", default=False,
                 help="use LaTeX for plot text rendering (you must have LaTeX installed on your system!)")
    p.add_option("-t", "--timings", dest="timings", action="store_true", default=False,
                 help="for every plot, include a complimentary timing plot")
    p.add_option("-n", "--non-normalized", dest="non_normalized", action="store_true", default=False,
                 help="includes plot with non-normalized y-axes")
    p.add_option("-m", "--median-count", dest="median_count", action="store_true", default=False,
                 help="plot the number of medians")
    p.add_option("-M", "--worst-median", dest="worst_median", action="store_true", default=False,
                 help="plot the worst medians")
    p.add_option("-s", "--strip-mprs", dest="strip_mprs", help="ignore any file with fewer mprs than this number",
                 metavar="MPR-MIN", default="0")
    p.add_option("-e", "--strip-equal", dest="strip_equal", help="ignore any file where mprs == median",
                 action="store_true", default=False)
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
    strip_mprs = int(options.strip_mprs)
    diameter = options.diameter
    median_count = options.median_count
    worst_median = options.worst_median
    strip_equal = options.strip_equal

    plot_types = []

    if diameter:
        plot_types += ["Diameter"]
    if zero_loss:
        plot_types += ["Zero Loss Diameter"]
    if median_count:
        plot_types += ["Median Count"]
    if worst_median:
        plot_types += ["Worst Median Diameter"]


    if latex and not plot:
        print "Warning: option '-l' (--use-latex) has no effect without option '-p' (--plot)!"
    if not os.path.isfile(file):
        p.error("File not found, '{0}'. Please be sure you typed the name correctly!".format(file))
    else:
        findExtrema(file, plot_types, non_normalized, timings, plot, latex, strip_mprs, strip_equal)
        if compare_file is not None:
            compare_logs(file, compare_file)
        if check is not None:
            check_files(file, check)
        if plot:
            plt.show()


if __name__ == "__main__":
    main()