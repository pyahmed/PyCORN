#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import AutoMinorLocator
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
from pycorn import pc_res3

pcscript_version = 0.13

parser = argparse.ArgumentParser(
    description = "Extract data from UNICORN .res files to .csv/.txt and plot them (matplotlib required)",
    epilog = "Make it so!")
parser.add_argument("-c", "--check", 
                    help = "Perform simple check if file is supported",
                    action = "store_true")
parser.add_argument("-n", "--info", 
                    help = "Display entries in header",
                    action = "store_true")
parser.add_argument("-i", "--inject", type = int, default = None, 
                    help = "Set injection number # as zero retention, use -t to find injection points",
                    metavar="#")
parser.add_argument("-r", "--reduce", type = int, default = 1,
                    help = "Write/Plot only every n sample",
                    metavar="#")
parser.add_argument("-t", "--points", 
                    help = "Display injection points",
                    action = "store_true")

group0 = parser.add_argument_group('Extracting', 'Options for writing csv/txt files')
group0.add_argument("-e", "--extract", 
                    help = "Write csv file for supported data blocks",
                    action = "store_true")

group1 = parser.add_argument_group('Plotting', 'Options for plotting')
group1.add_argument("-p", "--plot", 
                    help = 'Plot curves',
                    action = "store_true")
group1.add_argument("--no_fractions", 
                    help="Disable plotting of fractions",
                    action = "store_false")
group1.add_argument("--xmin", type = float, default=None,
                    help="Lower bound on the x-axis",
                    metavar="#")
group1.add_argument("--xmax", type = float, default=None,
                    help="Upper bound on the x-axis",
                    metavar="#")
group1.add_argument("--par1", type = str, default='Cond',
                    help="Data for 2nd y-axis (Default=Cond)")
group1.add_argument("--par2", type = str, default=None,
                    help="Data for 3rd y-axis (Default=None)")                 
group1.add_argument('-f', '--format', type = str,
                    choices=['svg','svgz','tif','tiff','jpg','jpeg',
                    'png','ps','eps','raw','rgba','pdf','pgf'],
                    default = 'pdf',
                    help = "File format of plot files (default: pdf)")
group1.add_argument('-d', '--dpi', default=300, type=int, 
					help="DPI (dots per inch) for raster images (png, jpg, etc.). Default is 300.")
parser.add_argument("-u", "--user", 
                    help = "Show stored user name",
                    action = "store_true")
parser.add_argument('--version', action='version', version=str(pcscript_version))
parser.add_argument("inp_res",
                    help="Input .res file(s)",
                    nargs='+',
                    metavar="<file>.res")

args = parser.parse_args()

def mapper(min_val, max_val, perc):
    '''
    calculate relative position in delta min/max
    '''
    x = abs(max_val - min_val) * perc
    if min_val < 0:
        return (x - abs(min_val))
    else:
        return (x + min_val)


def expander(min_val, max_val, perc):
    '''
    expand -/+ direction of two values by a percentage of their delta
    '''
    delta = abs(max_val - min_val)
    x = delta * perc
    return (min_val - x, max_val + x)


def xy_data(inp):
    '''
    Takes a data block and returns two lists with x- and y-data
    '''
    x_data = [x[0] for x in inp]
    y_data = [x[1] for x in inp]
    return x_data, y_data


def uvdata(inp):
    '''
    helps in finding the useful data
    '''
    UV_blocks = [i for i in inp if i.startswith('UV') or i.endswith('nm')]
    for i in UV_blocks:
        if i.endswith("_0nm"):
            UV_blocks.remove(i)


def smartscale(inp):
    '''
    input is the entire fdata block
    checks user input/fractions to determine scaling of x/y-axis
    returns min/max for x/y
    '''
    UV_blocks = [i for i in inp.keys() if i.startswith('UV') and not i.endswith('_0nm')]
    uv1_data = inp[UV_blocks[0]]['data']
    uv1_x, uv1_y = xy_data(uv1_data)
    try:
        uv2_data = inp[UV_blocks[1]]['data']
        uv2_x, uv2_y = xy_data(uv2_data)
        uv3_data = inp[UV_blocks[2]]['data']
        uv3_x, uv3_y = xy_data(uv3_data)
    except:
        KeyError
        uv2_data = None
        uv3_data = None
    frac_delta = []
    try:
        frac_data = inp['Fractions']['data']
        frac_x, frac_y = xy_data(frac_data)
        frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
        frac_delta.append(frac_delta[-1])
    except:
        KeyError
        frac_data = None
    if args.xmin:
        plot_x_min = args.xmin
    else:
        if frac_data:
            plot_x_min = frac_data[0][0]
        else:
            plot_x_min = uv1_x[0]
    if args.xmax:
        plot_x_max = args.xmax
    else:
        if frac_data:
            plot_x_max = frac_data[-1][0] + frac_delta[-1]*2 # recheck
        else:
            plot_x_max = uv1_x[-1]
    if plot_x_min > plot_x_max:
        print("Warning: xmin bigger than xmax - adjusting...")
        plot_x_min = uv1_x[0]
    if plot_x_max < plot_x_min:
        print("Warning: xmax smaller than xmin - adjusting...")
        plot_x_max = uv1_x[-1]
    # optimize y_scaling
    min_y_values = []
    max_y_values = []
    for i in UV_blocks:
        tmp_x, tmp_y = xy_data(inp[i]['data'])
        range_min_lst = [abs(a - plot_x_min) for a in tmp_x]
        range_min_idx = range_min_lst.index(min(range_min_lst))
        range_max_lst = [abs(a - plot_x_max) for a in tmp_x]
        range_max_idx = range_max_lst.index(min(range_max_lst))
        values_in_range = tmp_y[range_min_idx:range_max_idx]
        min_y_values.append(min(values_in_range))
        max_y_values.append(max(values_in_range))
    plot_y_min_tmp = min(min_y_values)
    plot_y_max_tmp = max(max_y_values)
    plot_y_min, plot_y_max = expander(plot_y_min_tmp, plot_y_max_tmp, 0.085)
    return plot_x_min, plot_x_max, plot_y_min, plot_y_max

def plotterX(inp,fname):
    plot_x_min, plot_x_max, plot_y_min, plot_y_max = smartscale(inp)
    host = host_subplot(111, axes_class=AA.Axes)
    host.set_xlabel("Elution volume (ml)")
    host.set_ylabel("Absorbance (mAu)")
    host.set_xlim(plot_x_min, plot_x_max)
    host.set_ylim(plot_y_min, plot_y_max)
    for i in inp.keys():
        if i.startswith('UV') and not i.endswith('_0nm'):
            x_dat, y_dat = xy_data(inp[i]['data'])
            print("Plotting: " + inp[i]['data_name'])
            stl = styles[i[:4]]
            p0, = host.plot(x_dat, y_dat, label=inp[i]['data_name'], color=stl['color'],
                            ls=stl['ls'], lw=stl['lw'],alpha=stl['alpha'])
    if args.par1:
        par1_inp = args.par1
        par1 = host.twinx()
        par1_data = inp[par1_inp]
        stl = styles[par1_inp]
        par1.set_ylabel(par1_data['data_name'] + " (" + par1_data['unit'] + ")", color=stl['color'])
        x_dat_p1, y_dat_p1 = xy_data(par1_data['data'])
        p1_ymin, p1_ymax = expander(min(y_dat_p1), max(y_dat_p1), 0.085)
        par1.set_ylim(p1_ymin, p1_ymax)
        print("Plotting: " + par1_data['data_name'])
        p1, = par1.plot(x_dat_p1, y_dat_p1, label=par1_data['data_name'], color=stl['color'],
                        ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'])
    if args.par2:
        par2_inp = args.par2
        par2 = host.twinx()
        offset = 60
        new_fixed_axis = par2.get_grid_helper().new_fixed_axis
        par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))  
        par2.axis["right"].toggle(all=True)
        par2_data = inp[par2_inp]
        stl = styles[par2_inp]
        par2.set_ylabel(par2_data['data_name'] + " (" + par2_data['unit'] + ")", color=stl['color'])
        x_dat_p2, y_dat_p2 = xy_data(par2_data['data'])
        p2_ymin, p2_ymax = expander(min(y_dat_p2), max(y_dat_p2), 0.075)
        par2.set_ylim(p2_ymin, p2_ymax)
        print("Plotting: " + par2_data['data_name'])
        p2, = par2.plot(x_dat_p2, y_dat_p2, label=par2_data['data_name'], color=stl['color'],
                        ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'])
    if args.no_fractions:
        try:
            frac_data = inp['Fractions']['data']
            frac_x, frac_y = xy_data(frac_data)
            frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
            frac_delta.append(frac_delta[-1])
            frac_y_pos = mapper(host.get_ylim()[0], host.get_ylim()[1], 0.015)
            for i in frac_data:
                host.axvline(x=i[0], ymin=0.065, ymax=0.0, color='r', linewidth=0.85)
                host.annotate(str(i[1]), xy=(i[0] + frac_delta[frac_data.index(i)] * 0.55, frac_y_pos),
                         horizontalalignment='center', verticalalignment='bottom', size=8, rotation=90)
        except:
            KeyError
    host.set_xlim(plot_x_min, plot_x_max)
    host.legend(fontsize=8, fancybox=True, labelspacing=0.4, loc='upper right')
    host.xaxis.set_minor_locator(AutoMinorLocator())
    host.yaxis.set_minor_locator(AutoMinorLocator())
    plt.title(fname, loc='left')
    internal_run_name = str(inp['Logbook']['run_name'])
    plot_file = fname[:-4] + "_" + internal_run_name + "_plot." + args.format
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.dpi)
    print("Plot saved to: " + plot_file)
    plt.clf()

def data_writer1(fname, inp):
    '''
    writes sensor/run-data to csv-files
    '''
    run_name = inp['Logbook']['run_name']
    for i in inp.keys():
        print("Extracting: " + inp[i]['data_name'])
        outfile_base = fname[:-4] + "_" + run_name + "_" + inp[i]['data_name']
        type = inp[i]['data_type']
        if type == 'meta':
            data = inp[i]['data']
            data_to_write = data.encode('utf-8')
            ext = '.txt'
            sep = '\t'
            with open(outfile_base + ext, 'wb') as fout:
                fout.write(data_to_write)
        else:
            x_dat,y_dat = xy_data(inp[i]['data'])
            ext = '.csv'
            sep = ','
            with open(outfile_base + ext, 'wb') as fout:
                for x,y in zip(x_dat,y_dat):
                    dp = str(x) + sep + str(y) + str('\r\n')
                    data_to_write = dp.encode('utf-8')
                    fout.write(data_to_write)

                    
styles = {'UV':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV1_':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV2_':{'color': '#e51616', 'lw': 1.4, 'ls': "-", 'alpha':1.0},
'UV3_':{'color': '#c73de6', 'lw': 1.2, 'ls': "-", 'alpha':1.0},
'Cond':{'color': '#62181A', 'lw': 1.4, 'ls': "-", 'alpha':0.75},
'Conc':{'color': '#0F990F', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Pres':{'color': '#3E1719', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Temp':{'color': '#e04730', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Inje':{'color': '#d56d9d', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'pH':{'color': '#0C7F7F', 'lw': 1.0, 'ls': "-", 'alpha':0.75},}


def main2():
    for fname in args.inp_res:
        if args.inject == None:
            args.inject = -1
        fdata = pc_res3(fname, reduce = args.reduce, inj_sel=args.inject)
        fdata.load()
        if args.extract:
            data_writer1(fname, fdata)
        if args.check:
            fdata.input_check(show=True)
        if args.info:
            fdata.showheader()
        if args.points:
            fdata.inject_det(show=True)
        if args.user:
            user = fdata.get_user()
            print("User: " + user)
        if args.plot:
            plotterX(fdata, fname)

main2()
