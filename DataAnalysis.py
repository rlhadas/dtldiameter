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
    print "\tMin:\t{0}".format(min(list))
    print "\tMax:\t{0}".format(max(list))
    print "\tMedian:\t{0}".format(numpy.median(list))
    print "\tMean:\t{0}".format(numpy.mean(list))
    print ""


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


def findExtrema(csv_file, zero_loss, non_normalized, timings, plot, latex):
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
    diameter_list = []
    gene_count_list = []
    mpr_over_d_list = []
    diameter_over_gene_list = []
    mpr_over_normalized_d_list = []
    diameter_timings = []
    DP_timings = []
    total_timings = []
    number_list = []
    number = 0
    DTL = "n/a"

    with open(csv_file) as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) > 0 and row[0] != "File Name" and row[3] != "0":
                number += 1
                number_list += [number]
                DTL = row[1]
                mpr = float(row[2])
                diameter = int(row[3])
                gene_count = int(row[4])
                gene_count = int(row[5]) #TODO reverse this
                DP_timings += [float(row[6])]
                diameter_timings += [float(row[7])]
                total_timings += [DP_timings[-1] + diameter_timings[-1]]
                mpr_list += [mpr]
                diameter_list += [diameter]
                gene_count_list += [gene_count]
                #mpr_over_d_list += [mpr/float(diameter)]
                diameter_over_gene_list += [diameter/float(gene_count)]
                #mpr_over_normalized_d_list += [mpr/(diameter/float(gene_count))]

        displayListValues(mpr_list, "MPR Count")
        displayListValues(diameter_list, "Diameter")
        displayListValues(gene_count_list, "Gene Tree Size")
        #displayListValues(mpr_over_d_list, "MPR Count/Diameter")
        displayListValues(diameter_over_gene_list, "Normalized Diameter")
        #displayListValues(mpr_over_normalized_d_list, "MPR Count/(Diameter/Gene Count)")
        displayListValues(diameter_timings, "Diameter Running Time (seconds)")
        displayListValues(DP_timings, "Reconciliation Running Time (seconds)")
        displayListValues(total_timings, "Total Running Time (seconds)")

    filepath, extension = os.path.splitext(csv_file)
    zero_loss_file = filepath + "_zl" + extension
    if plot:
        make_plot(csv_file, False, non_normalized, timings, gene_count_list, diameter_list, diameter_over_gene_list,
             mpr_list, DP_timings,
             diameter_timings, total_timings)
    if zero_loss:
        mpr_list = []
        diameter_list = []
        gene_count_list = []
        mpr_over_d_list = []
        diameter_over_gene_list = []
        mpr_over_normalized_d_list = []
        diameter_timings = []
        DP_timings = []
        total_timings = []
        number_list = []
        number = 0
        DTL = "n/a"
        with open(zero_loss_file) as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) > 0 and row[0] != "File Name" and row[3] != "0":
                    number += 1
                    number_list += [number]
                    DTL = row[1]
                    mpr = float(row[2])
                    diameter = int(row[3])
                    gene_count = int(row[4])
                    species_count = int(row[5])
                    DP_timings += [float(row[6])]
                    diameter_timings += [float(row[7])]
                    total_timings += [DP_timings[-1] + diameter_timings[-1]]
                    mpr_list += [mpr]
                    diameter_list += [diameter]
                    gene_count_list += [gene_count]
                    # mpr_over_d_list += [mpr/float(diameter)]
                    diameter_over_gene_list += [diameter / float(gene_count)]
                    # mpr_over_normalized_d_list += [mpr/(diameter/float(gene_count))]
        if plot:
            make_plot(csv_file, True, non_normalized, timings, gene_count_list, diameter_list, diameter_over_gene_list,
                 mpr_list, DP_timings,
                 diameter_timings, total_timings)


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

    ratio_ylim_b = -0.1
    ratio_ylim_t = 1.1

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
        findExtrema(file, zero_loss, non_normalized, timings, plot, latex)
        if compare_file is not None:
            compare_logs(file, compare_file)
        if check is not None:
            check_files(file, check)
        if plot:
            plt.show()


if __name__ == "__main__":
    main()