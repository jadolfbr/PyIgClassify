import datetime
import sys
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plot
from matplotlib.dates import MonthLocator, DateFormatter, WeekdayLocator

from src.tools.path import *


#Just a quick script for plotting increase in structures/data over time in matplotlib

def plot_date_stamp_data_simple(in_path, only_redundant=True):

    print "In File: "+in_path

    if only_redundant:
        types = ["redundant"]
    else:
        types = ["redundant", "nr_by_cdr"]

    out_dir = os.path.os.path.dirname(in_path)
    out_dir = out_dir+"/plots"
    print out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    INFILE = open_file(in_path, 'r')

    data = defaultdict(lambda: defaultdict(dict))

    PDBIDs = 1
    renum_ab_chains = 2
    total_chains = 3
    lambda_chains = 4
    kappa_chains = 5
    heavy_chains = 6
    cdrs = 7
    cdrs_known = 8
    cdrs_known_dih_cutoff = 9

    to_string = dict([(1, "PDBIDs"), (2, "renum_structures"), (3, "ab_chains"), (4, "lambda_chains"), (5, "kappa_chains"), (6, "heavy_chains"), (7, "total_cdrs"),
        (8, "total_cdrs_with_distances"), (9, "total_cdrs_within_cutoffs")])

    #data_types = ["PDBIDs", "renum_ab chains", "lambda", "kappa" "heavy" "cdrs" "cdrs_known" "cdrs_known_dih_cutoff"]
    for line in INFILE:
        line = line.strip()
        if not line: continue
        if line[0] == "#":continue

        lineSP = line.split()
        if lineSP[-1] in types:
            for i in range(1, 10):
                data[lineSP[-1]][i][lineSP[0]] = lineSP[i]

    #Plot Data:
    for type in data:

            for value_type in data[type]:

                out_full = out_dir+"/n_over_time_"+type+"."+to_string[value_type]+".pdf"
                print out_full

                dates_dt, dates_mpl, values = get_values(data, type, value_type)
                print repr(dates_mpl)

                print "Creating Fig: "+out_full
                create_fig(dates_dt, values, outpath=out_full, label = to_string[value_type])


            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
            all_types = defaultdict()

            #Setup each combined plot:

            #SuperImposed plot of lambda, kappa, all.
            #types1 = [lambda_chains, kappa_chains, heavy_chains, total_chains]
            types1 = [ heavy_chains, kappa_chains, lambda_chains, "PyIgClassify Antibody Chains"]
            all_types["comparison_chains"] = types1

            #SuperImposed plot of CDRs, CDRs known, CDRs with outliers
            types2 = [cdrs, cdrs_known, cdrs_known_dih_cutoff, "PyIgClassify CDR Entries"]
            all_types["comparison_total_cdrs"] = types2

            #SuperImposed plot of Total PDBIds, Total ab structures, Total chains.
            types3 = [total_chains, renum_ab_chains, PDBIDs, "PyIgClassify Entries" ]
            all_types["comparison_general"] = types3

            for comp in all_types:

                outpath = out_dir+"/"+type + "_"+comp+".pdf"

                fig=None
                ax=None
                for i in range(0, len(all_types[comp]) - 1):
                    value_type = all_types[comp][i]
                    color = colors[i]
                    print repr(value_type)
                    print repr(to_string)
                    named_type = to_string[value_type]

                    dates_dt, dates_mpl, values = get_values(data, type, value_type)
                    fig, ax = create_fig(dates_dt, values, color=color, fig=fig, ax=ax, label=named_type, ylabel="total")
                    #ax.legend(loc='best', fancybox=True)


                    # Shrink current axis by 20%
                    #box = ax.get_position()
                    #ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])

                    # Put a legend to the right of the current axis
                    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))




                #plot.savefig(outpath, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')

                print "Creating Fig: "+outpath
                lgd = plot.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
                plot.title(all_types[comp][-1])
                plot.savefig(outpath, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')


def create_fig(dates_dt, values, linestyle="--", marker="^", color="b", fig=None, ax=None, outpath=None, label=None, ylabel=None):


            label = label.replace("_", " ")
            if not fig:
                fig, ax = plot.subplots()

            if not ylabel:
                plot.ylabel(label, fontweight="bold")
            else:
                plot.ylabel(ylabel, fontweight="bold")

            plot.xlabel("date", fontweight="bold")
            plot.xticks(dates_dt, rotation = 'vertical' )

            ax.plot_date(dates_dt, values, linestyle=linestyle, marker=marker, color=color, label = label)

            fmt = DateFormatter("%b '%y")

            min_dt = min(dates_dt)
            max_dt = max(dates_dt)

            #print repr(min_dt)
            #print repr(max_dt)

            left_date = datetime.datetime(min_dt.year, min_dt.month, 1)
            right_date = datetime.datetime(max_dt.year, max_dt.month+1, 1)

            print repr(left_date)
            print repr(right_date)

            locator = MonthLocator()
            min_locator = WeekdayLocator()
            ax.xaxis.set_major_locator(locator)
            ax.set_xlim(left_date, right_date)
            #ax.set_ylim(0, max(values)+1500)
            #ax.xaxis.set_minor_locator(min_locator)


            ax.xaxis.set_major_formatter(fmt)
            ax.autoscale_view()

            #plot.axis('tight')
            fig.autofmt_xdate()

            if outpath:
                fig.savefig(outpath, dpi=300)

            return fig, ax


            #plot.close()

def get_values(data, db_name, value_type):

        dates_dt = []
        for s in sorted( data[db_name][value_type], key=lambda x: datetime.datetime.strptime(x, '%m/%d/%Y') ):
            print s
            if data[db_name][value_type][s] == "NA": continue
            sSP = s.split("/")
            dates_dt.append(datetime.datetime(int(sSP[2]), int(sSP[0]), int(sSP[1])))

        dates_mpl = mpl.dates.date2num(dates_dt)

        values = []
        for key in sorted( data[db_name][value_type], key=lambda x: datetime.datetime.strptime(x, '%m/%d/%Y') ):
            if data[db_name][value_type][key] == "NA": continue
            values.append(int(data[db_name][value_type][key]))

        return dates_dt, dates_mpl, values

if __name__ == "__main__":

    #Example Command-line: python plot_logs_in_matplotlib.py ../../DBOUT/logs/reported_overall_totals.txt 1

    in_path = sys.argv[1]

    only_redundant = int(sys.argv[2])


    plot_date_stamp_data_simple(in_path, only_redundant)