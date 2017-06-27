import csv
import numpy
import matplotlib.pyplot as plt

def displayListValues(list, name):
    print name + ": "
    print "\tMin:\t{0}".format(min(list))
    print "\tMax:\t{0}".format(max(list))
    print "\tMedian:\t{0}".format(numpy.median(list))
    print "\tMean:\t{0}".format(numpy.mean(list))
    print ""

def findExtrema(csv_file):
    """Finds the minimums, maximums, medians, and means of the:
        MPR Count
        Diameter
        Gene Count
        MPR Count/Diameter
        Diameter/Gene Count
        MPR Count/(Diameter/Gene Count)
    Of a csv file created by Diameter.py"""

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

    with open(csv_file) as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] != "File Name": #and row[2] != "0":
                number += 1
                number_list += [number]
                mpr = float(row[1])
                diameter = int(row[2])
                gene_count = int(row[3])
                DP_timings += [float(row[4])]
                diameter_timings += [float(row[5])]
                total_timings += [DP_timings[-1] + diameter_timings[-1]]
                mpr_list += [mpr]
                diameter_list += [diameter]
                gene_count_list += [gene_count]
                #mpr_over_d_list += [mpr/float(diameter)]
                diameter_over_gene_list += [diameter/float(gene_count)]
                #mpr_over_normalized_d_list += [mpr/(diameter/float(gene_count))]

        displayListValues(mpr_list, "MPR Count")
        displayListValues(diameter_list, "Diameter")
        displayListValues(gene_count_list, "Gene Node Count")
        #displayListValues(mpr_over_d_list, "MPR Count/Diameter")
        displayListValues(diameter_over_gene_list, "Diameter/Gene Node Count")
        #displayListValues(mpr_over_normalized_d_list, "MPR Count/(Diameter/Gene Count)")
        displayListValues(diameter_timings, "Diameter Time (seconds)")
        displayListValues(DP_timings, "DP Time (seconds)")
        displayListValues(total_timings, "Total Time (seconds)")

    fig, ax = plt.subplots(ncols=2, nrows=2)
    fig.canvas.set_window_title("{0} Diameters Calculated Main".format(len(diameter_list)))
    diameter_hist = ax[0][1]
    mpr_diameter = ax[0][0]
    diameter = ax[1][0]
    norm_diameter = ax[1][1]
    diameter.scatter(gene_count_list, diameter_list, c=gene_count_list)
    diameter.set_xlabel("Gene Node Count")
    diameter.set_ylabel("Diameter")
    diameter.set_title("Diameter/Gene Count")
    diameter.grid()
    norm_diameter.scatter(gene_count_list, diameter_over_gene_list, c=gene_count_list)
    norm_diameter.set_xlabel("Gene Node Count")
    norm_diameter.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter.set_title("Normalized Diameter/Gene Count")
    norm_diameter.grid()
    #ax[1][0].scatter(gene_count_list, DP_timings, c=gene_count_list)
    #ax[1][0].set_xlabel("Gene Count")
    #ax[1][0].set_ylabel("DP Time (seconds)")
    #ax[1][0].grid()

    diameter_hist.hist(diameter_over_gene_list, 100, orientation='horizontal')
    diameter_hist.set_ylabel("Diameter (normalized to gene node count)")
    diameter_hist.set_xlabel("Number of Results")
    diameter_hist.set_title("Normalized Diameter Counts")
    diameter_hist.set_ylim(-0.1, 2.1)
    diameter_hist.grid()
    #ax[1][1].scatter(diameter_over_gene_list, diameter_timings, c=gene_count_list)
    #ax[1][1].set_xlabel("Diameter (normalized to gene count)")
    #ax[1][1].set_ylabel("Diameter Time (seconds)")
    #ax[1][1].grid()
    mpr_diameter.scatter(mpr_list, diameter_over_gene_list, c=gene_count_list)
    mpr_diameter.set_ylim(-0.1, 2.1)
    mpr_diameter.set_xlabel("MPR Count")
    mpr_diameter.set_ylabel("Diameter (normalized to gene node count)")
    mpr_diameter.set_title("Normalized Diameter/MPR Count")
    mpr_diameter.grid()
    mpr_diameter.set_xscale('log')

    plt.show()

    fig, ax = plt.subplots(ncols=3, nrows=1)
    fig.canvas.set_window_title("{0} Diameters Calculated Time Complexity".format(len(diameter_list)))
    DP_time = ax[0]
    diameter_time = ax[1]
    total_time = ax[2]
    diameter_time.scatter(gene_count_list, diameter_timings, c=diameter_over_gene_list)
    diameter_time.set_xlabel("Gene Node Count")
    diameter_time.set_ylabel("Diameter Time (seconds)")
    diameter_time.set_title("Diameter Time Complexity")
    diameter_time.grid()
    diameter_time.set_ylim(0.01, (10**5))
    diameter_time.set_yscale('log')
    DP_time.scatter(gene_count_list, DP_timings, c=diameter_over_gene_list)
    DP_time.set_xlabel("Gene Node Count")
    DP_time.set_ylabel("DP Time (seconds)")
    DP_time.set_title("DP Time Complexity")
    DP_time.grid()
    DP_time.set_yscale('log')
    DP_time.set_xscale('log')
    total_time.scatter(gene_count_list, total_timings, c=diameter_over_gene_list)
    total_time.set_xlabel("Gene Node Count")
    total_time.set_ylabel("Total Time (seconds)")
    total_time.set_title("Total Time Complexity")
    total_time.grid()
    total_time.set_ylim(0.01, (10**5))
    total_time.set_yscale('log')

    plt.show()

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

def t():
    findExtrema("COG_Pilot_Log_02.csv")
