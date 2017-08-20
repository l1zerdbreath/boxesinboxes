"""
The purpose of this script is to create nested box plots
to display statistics of a combined dataset in conjunction 
with the statistcs of the individual datasets that
contribute to it.

Command line call:
python nested_box_plots.py

Output:
Image of specified format in specified output directory
Format is one of the file extensions supported by the
    active backend. Most backends support png, pdf, 
    ps, eps, and svg file extensions
    
Configurables (listed under __main__):
Input and output data and locations
Figure size
Plot characteristics including scatter size and transparency,
box spacing, axis names, ranges, and font style

Corresponding author for code:
Heather Cronk <heather.q.cronk@gmail.com>

Corresponding author for science:

"""


import pandas as pd
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import setp
from collections import OrderedDict
import math

def setBoxColors(bp, idx, color):

    """
    Subroutine for setting box plot element colors
    """
    setp(bp['boxes'][idx], color=color)
    setp(bp['caps'][2 * idx], color=color)
    setp(bp['caps'][2 * idx + 1], color=color)
    setp(bp['whiskers'][2 * idx], color=color)
    setp(bp['whiskers'][2 * idx + 1], color=color)
    setp(bp['medians'][idx], color=color)


if __name__ == "__main__":

    """
    Configurables
    
    """
    
    ### Data & Output Information ###
    data_dir = "/home/oco2/heather/liz/data"
    data_name = "dataforHeather.csv"
    output_dir = "/home/oco2/heather/liz/output"
    plot_name = "plot_ex1.png"
    outfile = os.path.join(output_dir, plot_name)
    
    ### Plot Infomation ###
    figsize = (15,5) # (x,y)
    out_format = "png"
    
    transparency_of_scatter = 1 #value between 0 and 1 where 0 is fully transparent and 1 is opaque
    scatter_size = 22
    
    y_axis_name = "# presynaptic structures"
    yrange = [] #[min,max], script will get it from the data if left empty
    
    first_small_box_pos = 0.75
    space_between_small_boxes = 0.5
    space_between_small_and_big_boxes = 0.25
    space_between_clusters = 1.25
    
    xtick_fontsize = 15
    ylabel_fontsize = 15
    xlabel_fontsize = 15
    ytick_fontsize = 12       
    
    #font weight options: 'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
    ylabel_fontweight = 'bold'
    ytick_fontweight = 'bold'
    xtick_fontweight = 'bold'
    

    """
    Prepare Data. 

    Dictionary variables:
    1) xidx is the position of the center of the box plot and corresponding scatter on the x axis
    2) width is the width of the box plot
    3) color is the hex color of the box and scatter (an rgb tuple or a named Python color would also work here)
    4) data is the data for the given population from the file
    """
    
    read_data = pd.read_csv(os.path.join(data_dir, data_name))
    
    cols = list(read_data.columns)
    data_name = list(set(cols) - set(["Animal", "Channel", "Geno"]))[0]
    
    number_of_small_boxes = len(set(read_data["Channel"])) * len(set(read_data["Animal"]))
    number_of_big_boxes = len(set(read_data["Channel"])) * len(set(read_data["Geno"]))
    total_boxes = number_of_small_boxes + number_of_big_boxes
    
    first_big_box_pos = number_of_big_boxes / len(set(read_data["Channel"]))
    series_addition = total_boxes / number_of_big_boxes
    idx_big_boxes = [first_big_box_pos + series_addition * n for n in xrange(0, number_of_big_boxes)]
    idx_big_boxes_plus_one = [i + 1 for i in idx_big_boxes]
    
    idx_start_clusters = [series_addition * n for n in xrange(0, number_of_big_boxes)]
    
    xidx = [first_small_box_pos] * total_boxes
    

    for n in xrange(1, len(xidx)):
        if n in idx_big_boxes or n in idx_big_boxes_plus_one:
	    xidx[n] = xidx[n-1] + space_between_small_and_big_boxes
	elif n in idx_start_clusters:
	    xidx[n] = xidx[n-1] + space_between_clusters
	else:
	    xidx[n] = xidx[n-1] + space_between_small_boxes


    vertical_lines_positions = []
    xtick_positions = []

    for i in idx_start_clusters[1::2]:
	xtick_positions.append(xidx[i - 1] + (xidx[i] - xidx[i-1])/2)
    for i in idx_start_clusters[2::2]:
	vertical_lines_positions.append(xidx[i - 1] + (xidx[i] - xidx[i-1])/2)
    
    mega_data_dict = OrderedDict([("VGLUT2",
                                             OrderedDict([("wt", 
                                                             OrderedDict([("A1",
		                                                                OrderedDict([("width", 0.2),
                                                                                            ("color", "#a63603"),
                                                                                            ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 1)]])),
                                                                                            ]
											   )
									   ),
									   ("A2",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#e6550d"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 2)]])),
                                                                                            ]
											   )
									   ),
									   ("All",
									        OrderedDict([("width", 2.),
                                                                                             ("color", "#000000"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Geno"] == "wt")]])),
                                                                                            ]
											   )
									   ),
									   ("A3",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#fd8d3c"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 3)]])),
                                                                                            ]
											   )
									   ),
									   ("A4",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#fdb385"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 4)]])),
                                                                                            ]
											   )
									   )
						                          ]
									 )
									 ),
							
							("ko", 
                                                             OrderedDict([("A0",
		                                                                OrderedDict([("width", 0.2),
                                                                                            ("color", "#54278f"),
                                                                                            ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 0)]])),
                                                                                            ]
											   )
									   ),
									   ("A5",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#756bb1"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 5)]])),
                                                                                            ]
											   )
									   ),
									   ("All",
									        OrderedDict([("width", 2.),
                                                                                             ("color", "#000000"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Geno"] == "ko")]])),
                                                                                            ]
											   )
									   ),
									   ("A6",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#9e9ac8"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 6)]])),
                                                                                            ]
											   )
									   ),
									   ("A7",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#acbfef"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "VGLUT2") & (read_data["Animal"] == 7)]])),
                                                                                            ]
											   )
									   )
						                          ])	  
						             )
							])		  
						       ),
				  ("GlyT2",
                                             OrderedDict([("wt", 
                                                             OrderedDict([("A1",
		                                                                OrderedDict([("width", 0.2),
                                                                                            ("color", "#a63603"),
                                                                                            ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 1)]])),
                                                                                            ]
											   )
									   ),
									   ("A2",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#e6550d"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 2)]])),
                                                                                            ]
											   )
									   ),
									   ("All",
									        OrderedDict([("width", 2.),
                                                                                             ("color", "#000000"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Geno"] == "wt")]])),
                                                                                            ]
											   )
									   ),
									   ("A3",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#fd8d3c"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 3)]])),
                                                                                            ]
											   )
									   ),
									   ("A4",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#fdb385"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 4)]])),
                                                                                            ]
											   )
									   )
						                          ]
									 )
									 ),
							
							("ko", 
                                                             OrderedDict([("A0",
		                                                                OrderedDict([("width", 0.2),
                                                                                            ("color", "#54278f"),
                                                                                            ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 0)]])),
                                                                                            ]
											   )
									   ),
									   ("A5",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#756bb1"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 5)]])),
                                                                                            ]
											   )
									   ),
									   ("All",
									        OrderedDict([("width", 2.),
                                                                                             ("color", "#000000"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Geno"] == "ko")]])),
                                                                                            ]
											   )
									   ),
									   ("A6",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#9e9ac8"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 6)]])),
                                                                                            ]
											   )
									   ),
									   ("A7",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#acbfef"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GlyT2") & (read_data["Animal"] == 7)]])),
                                                                                            ]
											   )
									   )
						                          ])	  
						             )
							])		  
						       ),
				  ("GAD67",
                                             OrderedDict([("wt", 
                                                             OrderedDict([("A1",
		                                                                OrderedDict([("width", 0.2),
                                                                                            ("color", "#a63603"),
                                                                                            ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 1)]])),
                                                                                            ]
											   )
									   ),
									   ("A2",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#e6550d"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 2)]])),
                                                                                            ]
											   )
									   ),
									   ("All",
									        OrderedDict([("width", 2.),
                                                                                             ("color", "#000000"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Geno"] == "wt")]])),
                                                                                            ]
											   )
									   ),
									   ("A3",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#fd8d3c"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 3)]])),
                                                                                            ]
											   )
									   ),
									   ("A4",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#fdb385"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 4)]])),
                                                                                            ]
											   )
									   )
						                          ]
									 )
									 ),
							
							("ko", 
                                                             OrderedDict([("A0",
		                                                                OrderedDict([("width", 0.2),
                                                                                            ("color", "#54278f"),
                                                                                            ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 0)]])),
                                                                                            ]
											   )
									   ),
									   ("A5",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#756bb1"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 5)]])),
                                                                                            ]
											   )
									   ),
									   ("All",
									        OrderedDict([("width", 2.),
                                                                                             ("color", "#000000"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Geno"] == "ko")]])),
                                                                                            ]
											   )
									   ),
									   ("A6",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#9e9ac8"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 6)]])),
                                                                                            ]
											   )
									   ),
									   ("A7",
									        OrderedDict([("width", 0.2),
                                                                                             ("color", "#acbfef"),
                                                                                             ("data", list(read_data[data_name][read_data.index[(read_data["Channel"] == "GAD67") & (read_data["Animal"] == 7)]])),
                                                                                            ]
											   )
									   )
						                          ])	  
						             )
							])		  
						       ),
				 ]
				)
    
    

    """
    Iterate over dictionary to get all the box plot info ready
    """
    data_to_plot = []
    color_of_plot = []
    widths = []
    xlabels = list(set(read_data["Channel"]))
    
    vertical_lines_ymin = 1000000.
    vertical_lines_ymax = -999999.
    
    for channel in xlabels:
	for geno in mega_data_dict[channel].keys():
            for animal in mega_data_dict[channel][geno].keys():
		data_to_plot.append(mega_data_dict[channel][geno][animal]["data"])
		vertical_lines_ymax = max(vertical_lines_ymax, max(mega_data_dict[channel][geno][animal]["data"]))
		vertical_lines_ymin = min(vertical_lines_ymin, min(mega_data_dict[channel][geno][animal]["data"]))
		color_of_plot.append(mega_data_dict[channel][geno][animal]["color"])
		widths.append(mega_data_dict[channel][geno][animal]["width"])
    
    vertical_lines_ymin = min(0, math.floor(vertical_lines_ymin))
    vertical_lines_ymax = math.ceil(vertical_lines_ymax)
    
    
    """
    Plot the box plots
    """

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)

    bp = ax.boxplot(data_to_plot, positions=xidx, widths=widths)
    
    
    """
    Set colors of box plots
    """
    for idx, color in enumerate(color_of_plot):
	setBoxColors(bp, idx, color)

    
    """
    Iterate over the dictionary again to do the scatter plots
    """
    cnt = 0
    for channel in xlabels:
	for geno in mega_data_dict[channel].keys():
            for animal in mega_data_dict[channel][geno].keys():
		if animal == "All":
		    #Skip scatter for the full channel + geno boxes
		    cnt += 1
	            continue
		y = mega_data_dict[channel][geno][animal]["data"]
		x = [xidx[cnt]] * len(y)
		color = mega_data_dict[channel][geno][animal]["color"]
		ax.scatter(x, y, c=color, edgecolor="None", alpha=transparency_of_scatter, s=scatter_size)
		cnt += 1

    
    """
    Set other plot parameters
    """
    
    if not yrange or len(yrange) != 2:
        ymin = vertical_lines_ymin
	ymax = vertical_lines_ymax
    else:
        ymin = yrange[0]
	ymax = yrange[1]
    ax.vlines(vertical_lines_positions, ymin, ymax, linestyle='dotted')
    mpl.rcParams['font.sans-serif']='Arial'
    
    ax.set_ylabel(y_axis_name, fontweight=ylabel_fontweight, fontsize=ylabel_fontsize)
    for t in ax.get_yticklabels():
        t.set_weight(ytick_fontweight)
	t.set_fontsize(ytick_fontsize)
    
    ax.set_xlabel("")
    plt.xticks(xtick_positions, xlabels)
    for t in ax.get_xticklabels():
        t.set_weight(xtick_fontweight)
	t.set_fontsize(xtick_fontsize)    
    
    """
    Display the plot
    """
    #plt.show()
    
    
    """
    Save the plot
    """
    plt.savefig(outfile, format=out_format, dpi=1000)
