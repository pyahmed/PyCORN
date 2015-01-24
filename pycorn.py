#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PyCORN - script to extract data from .res (results) files generated
by UNICORN Chromatography software supplied with ÄKTA Systems
(c)2014 - Yasar L. Ahmed
v0.12
'''

from __future__ import print_function

from collections import OrderedDict
import struct
import argparse
import codecs
import os
import time
try:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
except ImportError:
    plt = None

parser = argparse.ArgumentParser(
    description = "Extract data from UNICORN .res files to .csv/.txt and plot them (matplotlib required)",
    epilog = "Make it so!")
parser.add_argument("-c", "--check",
                    help = "Perform simple check if file is supported",
                    action = "store_true")
parser.add_argument("-n", "--info",
                    help = "Display entries in header",
                    action = "store_true")
parser.add_argument("-i", "--inject", type = int, default = -1,
                    help = "Set injection number # as zero retention, use -t to find injection points",
                    metavar="#")
parser.add_argument("-r", "--reduce", type = int, default = 1,
                    help = "Write/Plot only every n sample",
                    metavar="#")
parser.add_argument("-t", "--points",
                    help = "Display injection points",
                    action = "store_true")
parser.add_argument("-e", "--extract",
                    help = "Extract supported data blocks",
                    action = "store_true")
group1 = parser.add_argument_group('Plotting', 'Options for plotting')
group1.add_argument("-p", "--plot",
                    help = 'Plot curves',
                    action = "store_true")
group1.add_argument("-b", "--begin", type = float, default=None,
                    help = "Start point for plotting (in ml)",
                    metavar="#")
group1.add_argument("-s", "--finish", type = float, default=None,
                    help = "End point for plotting (in ml)",
                    metavar="#")
group1.add_argument('-f', '--format', type = str,
                    choices=['svg','svgz','tif','tiff','jpg','jpeg','png','ps','eps','raw','rgba','pdf','pgf'],
                    default = 'pdf',
                    help = "File format of plot files (default: pdf)")
parser.add_argument("-u", "--user",
                    help = "Show stored user name",
                    action = "store_true")
parser.add_argument("inp_res",
                    help = "Input .res file",
                    metavar = "<file>.res")




class Pycorn(OrderedDict):
    """A class for holding the PyCORN data.

    A subclass of `dict`, with the form `data_name`: `data`.
    """

    # first, some magic numbers
    RES_magic_id = b'\x11\x47\x11\x47\x18\x00\x00\x00\xB0\x02\x00\x00\x20\x6C\x03\x00'
    CNotes_id = b'\x00\x00\x01\x00\x02\x00\x03\x22'
    Methods_id = b'\x00\x00\x01\x00\x02\x00\x01\x02'
    Logbook_id = b'\x00\x00\x01\x00\x04\x00\x48\x04'
    Logbook_id2 = b'\x00\x00\x01\x00\x04\x00\x49\x04'
    SensData_id = b'\x00\x00\x01\x00\x04\x00\x01\x14'
    SensData_id2 = b'\x00\x00\x01\x00\x04\x00\x02\x14'
    Fractions_id = b'\x00\x00\x01\x00\x04\x00\x44\x04'
    Fractions_id2 = b'\x00\x00\x01\x00\x04\x00\x45\x04'
    Inject_id = b'\x00\x00\x01\x00\x04\x00\x46\x04'
    Inject_id2 = b'\x00\x00\x01\x00\x04\x00\x47\x04'
    LogBook_id = b'\x00\x00\x01\x00\x02\x00\x01\x13' # capital B!

    def __init__(self, file_name, reduce=1, default_inject=1):
        OrderedDict.__init__(self)
        self.file_name = file_name
        self.file_base = file_name[:-4]
        self.reduce = reduce
        self.default_inject = default_inject
        self.inj_sel = 0.0   # injection point is by default 0.0 ml
        self.header_read = False

        with open(self.file_name, 'rb') as f:
            self.raw_data = f.read()

    def input_check(self, show=False):
        '''
        Checks if input file is a supported res file
        x = magic number, y = version string, z = file size/EOF

        Returns True or False
        '''
        if show: print((" ---- \n Input file: {0}").format(self.file_name))

        x=self.raw_data.find(self.RES_magic_id, 0, 16)
        y=self.raw_data.find(b'UNICORN 3.10',16,36)
        z=struct.unpack("i", self.raw_data[16:20])

        if (x,y) == (0,24):
            if show: print(" Input is UNICORN 3.10 file!")
            x,y = (0,0)
        else:
            if show: print(" Input is not UNICORN 3.10 file!")
            x,y = (1,1)

        if z[0] == os.path.getsize(self.file_name):
            if show: print(" File size check - OK")
            z = 0
        else:
            if show: print(" File size mismatch - file corrupted?")
            z = 1
        if (x,y,z) != (0,0,0):
            if show: print("\n File not supported - stop!")
            return False
        else:
            if show: print("\n Alles safe! - Go go go!")
            return True


    def readheader(self):
        '''
        Extracts all the entries/declarations in the header (starts at position 686)
        '''

        # we only need to do this once
        if self.header_read: return
        self.header_read = True

        fread = self.raw_data
        header_end = fread.find(self.LogBook_id) + 342
        for i in range(686,header_end, 344):
            decl = struct.unpack("8s296s4i", fread[i:i+320])
            full_label = codecs.decode(decl[1], 'iso8859-1').rstrip("\x00")
            if full_label.find(':') == -1:
                r_name = ''
                d_name = full_label
            else:
                r_name = full_label[:full_label.find(':')]
                d_name = full_label[full_label.find('_')+1:]
            x = dict(magic_id = decl[0],
                     run_name = r_name,
                     data_name = d_name,
                     d_size = decl[2],
                     off_next = decl[3],
                     adresse = decl[4],
                     off_data = decl[5],
                     d_start = decl[4] + decl[5],
                     d_end = decl[4] + decl[2])
            name = x['data_name']
            dat = self.get(name, dict())
            dat.update(x)
            self[name] = dat

    def showheader(self, full=True):
        '''
        Prints content of header
        '''
        self.readheader()
        print((" ---- \n Header of {0}: \n").format(self.file_name))
        if full:
            print("  MAGIC_ID, ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
        else:
            print(" ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
        for name, dat in self.items():
            if full:
                print(" ", i['magic_id'],i['data_name'],i['d_size'],i['off_next'],i['adresse'],i['off_data'])
            else:
                print("",i['data_name'], i['d_size'], i['off_next'], i['adresse'], i['off_data'])


    def get_user(self):
        '''
        Show stored user name
        '''
        fread = self.raw_data[:512]
        u = struct.unpack("40s", fread[118:158])
        dec_u = codecs.decode(u[0], 'iso8859-1').rstrip("\x00")
        return dec_u


    def dataextractor(self, dat, show=False):
        '''
        Identify data type by comparing magic id, then run appropriate
        function to extract data, update orig. dict to include new data
        '''
        meta1 = [
            self.Logbook_id, self.Logbook_id2,
            self.Inject_id, self.Inject_id2,
            self.Fractions_id, self.Fractions_id2]
        meta2 = [self.CNotes_id, self.Methods_id]
        sensor = [self.SensData_id, self.SensData_id2]
        if dat['d_size'] == 0:
            pass
        elif dat['magic_id'] in meta1:
            dat.update(data=self.meta1_read(dat, show=show))
            return dat
        elif dat['magic_id'] in meta2:
            dat.update(data=self.meta2_read(dat, show=show))
            return dat
        elif dat['magic_id'] in sensor:
            values,unit=self.sensor_read(dat, show = show)
            dat.update(data=values, unit=unit)
            return dat


    def meta1_read(self, dat, show=False):
        '''
        Extracts meta-data/type1, Logbook, fractions and Inject marks
        for a specific datum
        '''
        if show:
            print((" Reading: {0}").format(dat['data_name']))
        final_data = []

        for i in range(dat['d_start'], dat['d_end'], 180):
            dp = struct.unpack("dd158s", self.raw_data[i:i+174])
            #acc_time = dp[0] # not used atm
            acc_volume = round(dp[1]-self.inj_sel,4)
            label = (codecs.decode(dp[2], 'iso8859-1')).rstrip('\x00')
            merged_data=acc_volume,label
            final_data.append(merged_data)
        return(final_data)


    def meta2_read(self, dat, show=False):
        '''
        Extracts meta-data/type2, Method/Program used in the run
        '''
        if show: print((" Reading: {0}").format(dat['data_name']))
        start, size = dat['d_start'], dat['d_size']
        tmp_data = self.raw_data[start:start+size]
        size = tmp_data.rfind(b'\n') # declared block-size in header is always off
                                     # by a few bytes, hence it is redetermined here
        if show and size != len(tmp_data):
            print('meta2: reevaluated size {} -> {}'.format(size, len(tmp_data)))

        raw_data = codecs.decode(self.raw_data[start:start+size], 'iso8859-1')
        if '\r' in raw_data:
            data = raw_data
        else:
            data = raw_data.replace('\n','\r\n')
        return data


    def sensor_read(self, dat, show=False):
        '''
        extracts sensor/run-data and applies correct division
        '''
        final_data = []
        if "UV" in dat['data_name'] or "Cond" == dat['data_name'] or "Flow" == dat['data_name']:
            sensor_div = 1000.0
        elif "Pressure" in dat['data_name']:
            sensor_div = 100.0
        else:
            sensor_div = 10.0
        if show: print((" Reading: {0}").format(dat['data_name']))

        fread = self.raw_data
        for i in range(dat['adresse']+207, dat['adresse']+222, 15):
            s_unit = struct.unpack("15s", fread[i:i+15])
            s_unit_dec = (codecs.decode(s_unit[0], 'iso8859-1')).rstrip('\x00')
        for i in range(dat['d_start'], dat['d_end'], 8):
            sread = struct.unpack("ii", fread[i:i+8])
            data=round((sread[0]/100.0)-self.inj_sel,4),sread[1]/sensor_div
            final_data.append(data)
        return(final_data[0::self.reduce],s_unit_dec)


    def inject_det(self, show=False):
        '''
        Finds injection points - required for adjusting retention volume
        '''
        injection_points = []
        self.readheader()
        inject_ids = [self.Inject_id, self.Inject_id2]
        for i in self.values():
            if i['magic_id'] in inject_ids:
                injection = self.meta1_read(i,show=show)[0][0]
                injection_points.append(injection)
                if injection != 0.0:
                    injection_points.insert(0,0.0)
        if injection_points == []:
            injection_points = [0.0]
            if show:
                print((" ---- \n Injection points: \n # \t ml \n 0 \t {0}").format(injection_points[0]))
            return(injection_points)
        else:
            if show:
                print(" ---- \n Injection points: \n # \t ml")
                for x,y in enumerate(injection_points):
                    print((" {0} \t {1}").format(x,y))
            return(injection_points)


    def load(self, show=False):
        '''
        extract all data and store in list
        '''
        injection_points = self.inject_det(show=show)
        try:
            self.inj_sel = injection_points[self.default_inject]
        except IndexError:
            if show:
                print("\n ERROR - Injection point does not exist! Selected default.\n")
            self.inj_sel = injection_points[-1]
        self.readheader()
        for name,dat in list(self.items()):
            dat = self.dataextractor(dat, show=show)
            if dat is not None:
                self[name] = dat
            else:
                #TODO: Maybe we should keep this around?
                del self[name]


    def writer(self, dat, show=False):
        '''
        general writer, picks the correct writer function based on data type
        '''
        if dat['magic_id'] == self.Methods_id or dat['magic_id'] == self.CNotes_id:
            self.meta_writer(dat, show=show)
        else:
            self.data_writer(dat, show=show)


    def meta_writer(self, dat, show=False, file_name = None):
        '''
        writes meta-data to txt-files
        '''
        if file_name is None:
            ext = "_"+dat['data_name']+".txt"
            file_name = self.file_base+ext
        with open(file_name,'wb') as fout:
            if show: print((" Writing {0}").format(dat['data_name']))
            data = (dat['data']).encode('utf-8')
            fout.write(data)


    def data_writer(self, dat, show=False):
        '''
        writes sensor/run-data to csv-files
        '''
        val_x = [str(x[0]) for x in dat['data']]
        if dat['data_name'] == 'Logbook':
            ext = '.txt'
            sep = '\t'
            val_y = [x[1] for x in dat['data']]
        else:
            ext = '.csv'
            sep = ','
            val_y = [str(x[1]) for x in dat['data']]
        fname = self.file_base+"_"+dat['run_name']+"_"+dat['data_name']+ext
        with open(fname, 'wb') as fout:
            if show: print((" Writing {0}").format(dat['data_name']))
            for x,y in zip(val_x,val_y):
                dp = (x + sep + y + str('\r\n')).encode('utf-8')
                fout.write(dp)

    @staticmethod
    def mapper(min_val,max_val,perc):
        '''
        calculate relative position in delta min/max
        '''
        x = abs(max_val-min_val)*perc
        if min_val < 0:
            return(x - abs(min_val))
        else:
            return(x + min_val)

    @staticmethod
    def expander(min_val,max_val,perc):
        '''
        expand -/+ direction of two values by a percentage of their delta
        '''
        delta = abs(max_val-min_val)
        x = delta*perc
        return(min_val - x, max_val + x)

    def plotter(self, dat, fractions, show=False,
                    begin=None, finish=None,
                    file_name=None, format='pdf'):
        '''
        plots a data series with fractions if present
        x-scaling is based on first/last fraction (if present)
        y-scaling is based on biggest/smallest y-value within x-range
        '''
        x_val = [x[0] for x in dat['data']]
        y_val = [x[1] for x in dat['data']]
        if fractions == None:
            plot_x_min = x_val[0]
            plot_x_max = x_val[-1]
        else:
            frac_vol = [x[0] for x in fractions]
            frac_delta = [abs(a - b) for a,b in zip(frac_vol,frac_vol[1:])]
            frac_delta.append(frac_delta[-1])
            plot_x_min = fractions[0][0]
            plot_x_max = fractions[-1][0] + frac_delta[-1]
        if begin != None:
            plot_x_min = begin
        if finish != None:
            plot_x_max = finish
        if plot_x_min != x_val[0]:
            range_min = x_val.index(min(x_val, key= lambda x:abs(x-plot_x_min)))
        else:
            range_min = 0
        if plot_x_max != x_val[-1]:
            range_max = x_val.index(min(x_val, key= lambda x:abs(x-plot_x_max)))
        else:
            range_max = -1
        x_val = x_val[range_min:range_max]
        y_val = y_val[range_min:range_max]
        plot_y_min = self.expander(min(y_val),max(y_val),0.085)[0]
        plot_y_max = self.expander(min(y_val),max(y_val),0.025)[1]
        if type(y_val[0]) == float or type(y_val[0]) == int:
            plt.xlim(xmin = plot_x_min, xmax = plot_x_max)
            plt.ylim(ymin = plot_y_min, ymax = plot_y_max)
            ax = plt.gca()
            frac_y_pos = self.mapper(ax.get_ylim()[0],ax.get_ylim()[1],0.015)
            plt.title(self.file_name, size=11)
            plt.ylabel(dat['unit'])
            plt.xlabel('ml')
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            plt.tick_params(axis='x',top='off',direction='out',length=6)
            plt.tick_params(axis='y',right='off',direction='out',length=6)
            plt.plot(x_val, y_val,linewidth=1.5, color='b')
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(axis='x',top='off',which='minor', direction='out', length=4, color='0.17')
            ax.tick_params(axis='y',right='off',which='minor', direction='out', length=4, color='0.17')
            if fractions != None:
                for i in fractions:
                    plt.axvline(x=i[0], ymin=0.065, ymax=0.0,color='r', linewidth=0.85)
                    plt.annotate(i[1], xy=(i[0]+frac_delta[fractions.index(i)]*0.5,frac_y_pos), horizontalalignment='center', verticalalignment='bottom', size=8, rotation=90)
            if show: print(" Plotting " + dat['data_name'])
            if file_name is None:
                ext = "." + format
                file_name = self.file_base + "_Plot_" + dat['run_name'] + "_" + dat['data_name'] + ext
            plt.savefig(file_name, papertype='a4',dpi=300)
            plt.clf()

def main():
    if plt is not None:
        print("\n Matplotlib found - Printing enabled\n")
    else:
        print("\n Matplotlib not found. Printing will not work!\n")

    args = parser.parse_args()
    pycorn = Pycorn(args.inp_res, reduce=args.reduce, default_inject=args.inject)

    plottable = [pycorn.SensData_id, pycorn.SensData_id2]
    start = time.time()
    if args.user:
        print((" ---- \n User: {0}").format(pycorn.showuser()))
    if args.check:
        pycorn.input_check(show=True)
    if args.info:
        pycorn.showheader(full=False)
    if args.points:
        pycorn.inject_det(show=True)
    if args.plot or args.extract:
        if pycorn.input_check() == True:
            pycorn.load(show=True)
        if args.extract:
            for i in pycorn.values():
                pycorn.writer(i, show=True)
        if args.plot and plt !=None:
            for dat in pycorn.values():
                if dat['magic_id'] in (pycorn.Fractions_id, pycorn.Fractions_id2):
                    fractions = dat['data']
                    break
            else:
                fractions = None

            for i in pycorn.values():
                if i['magic_id'] in plottable:
                    pycorn.plotter(i,fractions, show=True,
                        begin=args.begin, finish=args.finish,
                        format=args.format)

    end = time.time()
    runtime = str(end - start)
    print("\n ----")
    print((" Runtime: {0} seconds!").format(runtime[:5]))

if __name__ == '__main__':
    main()
