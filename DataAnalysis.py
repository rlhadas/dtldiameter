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
    number_list = []
    number = 0

    with open(csv_file) as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] != "File Name": #and row[2] != "0":
                number += 1
                number_list += [number]
                mpr = int(row[1])
                diameter = int(row[2])
                gene_count = int(row[3])
                DP_timings += [float(row[4])]
                diameter_timings += [float(row[5])]
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


    fig, ax = plt.subplots(ncols=2, nrows=2)
    fig.canvas.set_window_title("{0} Diameters Calculated".format(len(diameter_list)))
    ax[0][0].scatter(gene_count_list, diameter_list, c=gene_count_list)
    ax[0][0].set_xlabel("Gene Count")
    ax[0][0].set_ylabel("Diameter")
    ax[0][0].set_title("{0} Diameters Calculated".format(len(diameter_list)))
    ax[0][0].grid()
    ax[1][0].scatter(gene_count_list, diameter_over_gene_list, c=gene_count_list)
    ax[1][0].set_xlabel("Gene Count")
    ax[1][0].set_ylabel("Diameter (normalized to gene count)")
    ax[1][0].grid()
    #ax[1][0].scatter(gene_count_list, DP_timings, c=gene_count_list)
    #ax[1][0].set_xlabel("Gene Count")
    #ax[1][0].set_ylabel("DP Time (seconds)")
    #ax[1][0].grid()

    ax[0][1].hist(diameter_over_gene_list, 100)
    ax[0][1].set_xlabel("Diameter (normalized to gene count)")
    ax[0][1].set_ylabel("Number of Results")
    ax[0][1].grid()
    #ax[1][1].scatter(diameter_over_gene_list, diameter_timings, c=gene_count_list)
    #ax[1][1].set_xlabel("Diameter (normalized to gene count)")
    #ax[1][1].set_ylabel("Diameter Time (seconds)")
    #ax[1][1].grid()
    ax[1][1].scatter(diameter_list, diameter_over_gene_list, c=gene_count_list)
    ax[1][1].set_xlabel("Diameter")
    ax[1][1].set_ylabel("Diameter (normalized to gene count)")
    ax[1][1].grid()


    plt.show()

def t():
    findExtrema("COG_Pilot_Log_01.csv")
