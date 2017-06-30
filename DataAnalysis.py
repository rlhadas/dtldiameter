import csv
import numpy
import matplotlib.pyplot as plt
import os

def displayListValues(list, name):
    print name + ": "
    print "\tMin:\t{0}".format(min(list))
    print "\tMax:\t{0}".format(max(list))
    print "\tMedian:\t{0}".format(numpy.median(list))
    print "\tMean:\t{0}".format(numpy.mean(list))
    print ""

def findExtrema(csv_file, zero_loss):
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
                DP_timings += [float(row[5])]
                diameter_timings += [float(row[6])]
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

    size = 4
    color = 'black'
    if zero_loss:
        diameter_ylim_b = -0.05
        diameter_ylim_t = 1.05
    else:
        diameter_ylim_b = -0.1
        diameter_ylim_t = 2.1

    name = "Regular Loss "
    if zero_loss:
        name = "Zero Loss "

    fig, ax = plt.subplots(ncols=3, nrows=2)
    fig.canvas.set_window_title("{2}{0} Diameters {1} Main".format(len(diameter_list), DTL, name))
    norm_diameter_hist = ax[1][0]
    norm_mpr_diameter = ax[1][2]
    norm_diameter = ax[1][1]
    diameter = ax[0][1]
    diameter_hist = ax[0][0]
    mpr_diameter = ax[0][2]
    ax[1][0].set_ylabel("Normalized Diameter Counts")
    ax[0][0].set_ylabel("Diameter Counts")
    diameter.scatter(gene_count_list, diameter_list, c=color, s=size)
    diameter.set_xlabel("Gene Tree Size")
    diameter.set_title("Diameter Vs. Gene Count")
    diameter.set_ylim(-60, 1260)
    diameter.set_xlim(-100,2100)
    diameter.grid()

    norm_diameter.scatter(gene_count_list, diameter_over_gene_list, c=color, s=size)
    norm_diameter.set_xlabel("Gene Tree Size")
    #norm_diameter.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter.set_ylim(diameter_ylim_b, diameter_ylim_t)
    norm_diameter.set_xlim(-100,2100)
    norm_diameter.set_title("Normalized Diameter Vs. Gene Count")
    norm_diameter.grid()

    norm_diameter_hist.hist(diameter_over_gene_list, 100, orientation='horizontal')
    #norm_diameter_hist.set_ylabel("Diameter (normalized to gene node count)")
    norm_diameter_hist.set_xlabel("Number of Gene Families")
    norm_diameter_hist.set_title("Normalized Diameter Counts")
    norm_diameter_hist.set_ylim(diameter_ylim_b, diameter_ylim_t)
    norm_diameter_hist.grid()

    diameter_hist.hist(diameter_list, 100, orientation='horizontal')
    # diameter_hist.set_ylabel("Diameter")
    diameter_hist.set_xlabel("Number of Gene Families")
    diameter_hist.set_title("Diameter Counts")
    diameter_hist.set_ylim(-60, 1260)
    diameter_hist.grid()

    #ax[1][1].scatter(diameter_over_gene_list, diameter_timings, c=gene_count_list)
    #ax[1][1].set_xlabel("Diameter (normalized to gene count)")
    #ax[1][1].set_ylabel("Diameter Time (seconds)")
    #ax[1][1].grid()
    norm_mpr_diameter.scatter(mpr_list, diameter_over_gene_list, c=color, s=size)
    norm_mpr_diameter.set_ylim(diameter_ylim_b, diameter_ylim_t)
    norm_mpr_diameter.set_xlabel("MPR Count")
    norm_mpr_diameter.set_title("Normalized Diameter Vs. MPR Count")
    norm_mpr_diameter.grid()
    norm_mpr_diameter.set_xscale('log')

    mpr_diameter.scatter(mpr_list, diameter_list, c=color, s=size)
    mpr_diameter.set_ylim(-60, 1260)
    mpr_diameter.set_xlabel("MPR Count")
    mpr_diameter.set_title("Diameter Vs. MPR Count")
    mpr_diameter.grid()
    mpr_diameter.set_xscale('log')



    plt.show()
    fig, ax = plt.subplots(ncols=3, nrows=1)
    fig.canvas.set_window_title("{2}{0} Diameters {1} Time Complexity".format(len(diameter_list), DTL, name))
    DP_time = ax[0]
    diameter_time = ax[1]
    total_time = ax[2]
    diameter_time.scatter(gene_count_list, diameter_timings, c=diameter_over_gene_list, s=size)
    diameter_time.set_xlabel("Gene Tree Size")
    diameter_time.set_ylabel("Diameter Time (seconds)")
    diameter_time.set_title("Diameter Time Complexity")
    diameter_time.grid()
    diameter_time.set_ylim(0.01, (10**5))
    diameter_time.set_yscale('log')
    diameter_time.set_xscale('log')
    DP_time.scatter(gene_count_list, DP_timings, c=diameter_over_gene_list, s=size)
    DP_time.set_xlabel("Gene Tree Size")
    DP_time.set_ylabel("DP Time (seconds)")
    DP_time.set_title("DP Time Complexity")
    DP_time.grid()
    DP_time.set_yscale('log')
    DP_time.set_xscale('log')
    total_time.scatter(gene_count_list, total_timings, c=diameter_over_gene_list, s=size)
    total_time.set_xlabel("Gene Tree Size")
    total_time.set_ylabel("Total Time (seconds)")
    total_time.set_title("Total Time Complexity")
    total_time.grid()
    total_time.set_ylim(0.01, (10**5))
    total_time.set_yscale('log')
    total_time.set_xscale('log')

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
    findExtrema("COG_Pilot_Log_03.csv",False)
    findExtrema("ZL_COG_Pilot_Log_03.csv", True)


def check_files():

    data_files = map(lambda x: "TreeLifeData/" + x, os.listdir("TreeLifeData/"))
    with open("COG_Pilot_Log_02.csv") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] not in data_files:
                print "Duplicate: {0}".format(row[0])
            else:
                data_files.remove(row[0])
        print data_files

