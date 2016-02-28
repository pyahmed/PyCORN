#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PyCORN - script to extract data from .res (results) files generated
by UNICORN Chromatography software supplied with ÄKTA Systems
(c)2014-2016 - Yasar L. Ahmed
v0.18
'''

import argparse
from pycorn import pc_res3
from pycorn import pc_uni6

try:
    from mpl_toolkits.axes_grid1 import host_subplot
    from matplotlib.ticker import AutoMinorLocator
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt
    plotting = True
except:
    ImportError
    print("WARNING: Matplotlib not found - Plotting disabled!")
    plotting = False

try:
    import xlsxwriter
    xlsx = True
except:
    ImportError
    print("WARNING: xlsxwriter not found - xlsx-output disabled!")
    xlsx = False

pcscript_version = 0.14

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
group0.add_argument("-e", "--extract", type=str, choices=['csv','xlsx'],
                    help = "Write data to csv or xlsx file for supported data blocks")

group1 = parser.add_argument_group('Plotting', 'Options for plotting')
group1.add_argument("-p", "--plot", 
                    help = 'Plot curves',
                    action = "store_true")
group1.add_argument("--no_fractions", 
                    help="Disable plotting of fractions",
                    action = "store_true")
group1.add_argument("--no_inject", 
                    help="Disable plotting of inject marker(s)",
                    action = "store_true")
group1.add_argument("--no_legend", 
                    help="Disable legend for plot",
                    action = "store_true")
group1.add_argument("--no_title", 
                    help="Disable title for plot",
                    action = "store_true")
group1.add_argument("--xmin", type = float, default=None,
                    help="Lower bound on the x-axis",
                    metavar="#")
group1.add_argument("--xmax", type = float, default=None,
                    help="Upper bound on the x-axis",
                    metavar="#")
group1.add_argument("--par1", type = str, default='Cond',
                    help="Data for 2nd y-axis (Default=Cond), to disable 2nd y-axis, use --par1 None")
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
#args.no_inject
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
    try:
        frac_data = inp['Fractions']['data']
        frac_x, frac_y = xy_data(frac_data)
        frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
        frac_delta.append(frac_delta[-1])
    except:
        KeyError
        frac_data = None
    if args.xmin != None:
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
    if args.par1 == 'None':
        args.par1 = None
    if args.par1:
        try:
            par1_inp = args.par1
            par1 = host.twinx()
            par1_data = inp[par1_inp]
            stl = styles[par1_inp[:4]]
            par1.set_ylabel(par1_data['data_name'] + " (" + par1_data['unit'] + ")", color=stl['color'])
            x_dat_p1, y_dat_p1 = xy_data(par1_data['data'])
            p1_ymin, p1_ymax = expander(min(y_dat_p1), max(y_dat_p1), 0.085)
            par1.set_ylim(p1_ymin, p1_ymax)
            print("Plotting: " + par1_data['data_name'])
            p1, = par1.plot(x_dat_p1, y_dat_p1, label=par1_data['data_name'], 
            color=stl['color'], ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'])
        except:
            KeyError
            if par1_inp != None:
                print("Warning: Data block chosen for par1 does not exist!")
    if args.par2:
        try:
            par2_inp = args.par2
            par2 = host.twinx()
            offset = 60
            new_fixed_axis = par2.get_grid_helper().new_fixed_axis
            par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))  
            par2.axis["right"].toggle(all=True)
            par2_data = inp[par2_inp]
            stl = styles[par2_inp[:4]]
            par2.set_ylabel(par2_data['data_name'] + " (" + par2_data['unit'] + ")", color=stl['color'])
            x_dat_p2, y_dat_p2 = xy_data(par2_data['data'])
            p2_ymin, p2_ymax = expander(min(y_dat_p2), max(y_dat_p2), 0.075)
            par2.set_ylim(p2_ymin, p2_ymax)
            print("Plotting: " + par2_data['data_name'])
            p2, = par2.plot(x_dat_p2, y_dat_p2, label=par2_data['data_name'], 
            color=stl['color'],ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'])
        except:
            KeyError
            if par2_inp != None:
                print("Warning: Data block chosen for par2 does not exist!")
    if not args.no_fractions:
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
    if inp.inject_vol != 0.0:
        injections = inp.injection_points
        host.axvline(x=0, ymin=0.10, ymax=0.0, color='#FF3292',
                     ls ='-', marker='v', markevery=2, linewidth=1.5, alpha=0.85, label='Inject')
    host.set_xlim(plot_x_min, plot_x_max)
    if not args.no_legend:
        host.legend(fontsize=8, fancybox=True, labelspacing=0.4, loc='upper right', numpoints=1)
    host.xaxis.set_minor_locator(AutoMinorLocator())
    host.yaxis.set_minor_locator(AutoMinorLocator())
    if not args.no_title:
        plt.title(fname, loc='left', size=9)
    plot_file = fname[:-4] + "_" + inp.run_name + "_plot." + args.format
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.dpi)
    print("Plot saved to: " + plot_file)
    plt.clf()

def data_writer1(fname, inp):
    '''
    writes sensor/run-data to csv-files
    '''
    for i in inp.keys():
        print("Writing: " + inp[i]['data_name'])
        outfile_base = fname[:-4] + "_" + inp.run_name + "_" + inp[i]['data_name']
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

def generate_xls(inp, fname):
    '''
    Input = pycorn object
    output = xlsx file
    '''
    xls_filename = fname[:-4] + "_" + inp.run_name + ".xlsx"
    workbook = xlsxwriter.Workbook(xls_filename)
    worksheet = workbook.add_worksheet()
    writable_blocks = [inp.Fractions_id, inp.Fractions_id2, inp.SensData_id, inp.SensData_id2]
    d_list = []
    for i in inp.keys():
        if inp[i]['magic_id'] in writable_blocks:
            d_list.append(i)
    for i in d_list:
        dat = inp[i]['data']
        try:
            unit = inp[i]['unit']
        except:
            KeyError
            unit = 'Fraction'
        header1 = (inp[i]['data_name'], '')
        header2 = ('ml', unit)
        dat.insert(0, header1)
        dat.insert(1, header2)
        row = 0
        col = d_list.index(i) *2
        print("Writing: " + i)
        for x_val, y_val in (dat):
            worksheet.write(row, col, x_val)
            worksheet.write(row, col + 1, y_val)
            row += 1
    workbook.close()
    print("Data written to: " + xls_filename) 

                    
styles = {'UV':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV1_':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV2_':{'color': '#e51616', 'lw': 1.4, 'ls': "-", 'alpha':1.0},
'UV3_':{'color': '#c73de6', 'lw': 1.2, 'ls': "-", 'alpha':1.0},
'UV 1':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV 2':{'color': '#e51616', 'lw': 1.4, 'ls': "-", 'alpha':1.0},
'UV 3':{'color': '#c73de6', 'lw': 1.2, 'ls': "-", 'alpha':1.0},
'Cond':{'color': '#FF7C29', 'lw': 1.4, 'ls': "-", 'alpha':0.75},
'Conc':{'color': '#0F990F', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Pres':{'color': '#C0CBBA', 'lw': 1.0, 'ls': "-", 'alpha':0.50},
'Temp':{'color': '#b29375', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Inje':{'color': '#d56d9d', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'pH':{'color': '#0C7F7F', 'lw': 1.0, 'ls': "-", 'alpha':0.75},}


def main2():
    for fname in args.inp_res:
        if args.inject == None:
            args.inject = -1
        if (fname[-3:]).lower() == "zip":
            fdata = pc_uni6(fname)
            fdata.load()
            fdata.xml_parse()
            fdata.clean_up()
        if (fname[-3:]).lower() == "res":
            fdata = pc_res3(fname, reduce = args.reduce, inj_sel=args.inject)
            fdata.load()
        if args.extract == 'csv':
            data_writer1(fname, fdata)
        if args.extract == 'xlsx' and xlsx == True:
            generate_xls(fdata, fname)
        if args.check:
            fdata.input_check(show=True)
        if args.info:
            fdata.showheader()
        if args.points:
            fdata.inject_det(show=True)
        if args.user:
            user = fdata.get_user()
            print("User: " + user)
        if args.plot and plotting:
            plotterX(fdata, fname)

main2()
