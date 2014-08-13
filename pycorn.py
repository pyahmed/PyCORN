# -*- coding: utf-8 -*-
'''
PyCORN - script to extract data from .res (results) files generated 
by UNICORN Chromatography software supplied with ÄKTA Systems
(c)2014 - Yasar L. Ahmed
v0.1
'''

import struct
import argparse
import codecs
import os
import time

start = time.time()

parser = argparse.ArgumentParser(
    description = "Extract data from UNICORN .res files to .csv/.txt",
    epilog = "Make it so!")
parser.add_argument("-c", "--check", 
                    help = "Perform simple check if file is supported",
                    action = "store_true")
parser.add_argument("-n", "--info", 
                    help = "Display entries in header",
                    action = "store_true")
parser.add_argument("-i", "--inject", type = int, default = 0, 
                    help = "Set injection number # as zero retention, use -t to find injection points",
                    metavar="#")
parser.add_argument("-t", "--points", 
                    help = "Display injection points",
                    action = "store_true")
parser.add_argument("-e", "--extract", 
                    help = "Extract supported data blocks",
                    action = "store_true")
parser.add_argument("-p", "--plot", 
                    help = 'Plot curves',
                    action = "store_true")
parser.add_argument("-u", "--user", 
                    help = "Show stored user name",
                    action = "store_true")
parser.add_argument("inp_res",
                    help = "Input .res file",
                    metavar = "<file>.res",
                    type = argparse.FileType('rb'))

args = parser.parse_args()
file_in = args.inp_res.name # file handling is not done by argparse
file_base = args.inp_res.name[:-4]
args.inp_res.close() # see comment above

'''
Magic numbers
'''
RES_magic_id = b'\x11\x47\x11\x47\x18\x00\x00\x00\xB0\x02\x00\x00\x20\x6C\x03\x00'
CNotes_id = b'\x00\x00\x01\x00\x02\x00\x03\x22'
Methods_id = b'\x00\x00\x01\x00\x02\x00\x01\x02'
Logbook_id = b'\x00\x00\x01\x00\x04\x00\x48\x04'
SensData_id = b'\x00\x00\x01\x00\x04\x00\x01\x14'
Fractions_id = b'\x00\x00\x01\x00\x04\x00\x44\x04'
Inject_id = b'\x00\x00\x01\x00\x04\x00\x46\x04'
LogBook_id = b'\x00\x00\x01\x00\x02\x00\x01\x13' #capital B!

inj_sel = 0.0   #injection point is by default 0.0 ml


def input_check(inp):
    '''
    Checks if input file is a supported res file
    x = magic number, y = version string, z = file size/EOF
    '''
    print((" ---- \n Input file: {0}").format(inp))
    with open(inp, 'rb') as fo:
        file_read = fo.read()
        x=file_read.find(RES_magic_id, 0, 16)
        y=file_read.find(b'UNICORN 3.10',16,36)
        z=struct.unpack("i", file_read[16:20])
    if (x,y) == (0,24):
        print(" Input is UNICORN 3.10 file!")
        x,y = (0,0)
    else:
        print(" Input not UNICORN 3.10 file!")
        x,y = (1,1)
    if z[0] == os.path.getsize(file_in):
        print(" File size check - OK")
        z = 0
    else:
        print(" File size mismatch - file corrupted?")
        z = 1
    return(x,y,z)


def readheader(inp):
    '''
    Extracts all the entries in the header (starts at position 686)
    '''
    tmp_list = []
    with open(inp, 'rb') as fo:
        fread = fo.read()
        header_end = fread.find(LogBook_id) + 342
        for i in range(686,header_end, 344):
            decl = struct.unpack("8s296s4i", fread[i:i+320])
            x = decl[0], codecs.decode(decl[1], 'cp1250').rstrip("\x00"),decl[2],decl[3],decl[4],decl[5]
            tmp_list.append(x)
    return(tmp_list)


def showheader(inp,full="true"):
    '''
    Prints content of header
    '''
    header = readheader(file_in)
    print((" ---- \n Header of {0}: \n").format(inp))
    print("  ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
    for i in header:
        if full == "true":
            print(i)
        else:
            print(i[1:])


def showuser(inp):
    '''
    Show stored user name
    '''
    with open(inp, 'rb') as fo:
        fread = fo.read(512)
        u = struct.unpack("40s", fread[118:158])
        dec_u = codecs.decode(u[0], 'cp1250').rstrip("\x00")
        print((" ---- \n User: {0}").format(dec_u))


def dataextractor(in_tup):
    '''
    Identify data type by comparing magic id, then run appropriate
    function to extract data, return data as list
    '''
    inp_id = in_tup[0]
    data_tup = in_tup[1:]
    #skip empty blocks
    if (in_tup[2],in_tup[3]) == (0, 0):
        pass
    elif inp_id == CNotes_id or inp_id == Methods_id:
        x = meta2_read(data_tup)
        x.insert(0,inp_id)
        return x
    elif inp_id == Logbook_id or inp_id == Inject_id or Fractions_id == inp_id:
        return meta1_read(data_tup)
    elif inp_id == SensData_id:
        return sensor_read(data_tup)
    else:
        pass


def meta1_read(in_tup, silent="false"):
    '''
    Extracts meta-data/type1, Logbook, fractions and Inject marks
    '''
    if silent == "true":
        pass
    else:
        print((" Extracting: {0}").format(in_tup[0]))
    finaldata = []
    start = start = in_tup[3] + in_tup[4]
    end = in_tup[3] + in_tup[1]
    with open(file_in, 'rb') as fo:
        fread = fo.read()
        for i in range(start, end, 180):
            dp = struct.unpack("dd158s", fread[i:i+174])
            #acc_time = dp[0] # not used atm
            acc_volume = round(dp[1]-inj_sel,4)
            label = codecs.decode(dp[2], 'cp1250')
            merged_data=acc_volume,label.rstrip('\x00')#.encode("utf-8")
            finaldata.append(merged_data)
    finaldata.insert(0, in_tup[0])
    return(finaldata)


def meta2_read(in_tup):
    '''
    Extracts meta-data/type2, Method used in the run
    '''
    print((" Extracting: {0}").format(in_tup[0]))
    start = in_tup[3] + in_tup[4]
    with open(file_in,'rb') as fo:
        fo.seek(start)
        tmp_data = fo.read(in_tup[1])
        size = tmp_data.rfind(b'\n') #declared block-size in header is always off
        fo.seek(start)               #by a few bytes, hence it is redetermined here
        data = fo.read(size)
    output = [in_tup[0],data]
    return(output)


def sensor_read(in_tup):
    '''
    extracts sensor/run-data and applies correct division
    '''
    tmp_list0 = []
    # set positions
    start = in_tup[3] + in_tup[4]
    end = in_tup[3] + in_tup[1]
    # get the labels
    sep0 = in_tup[0].find(":")
    sep1 = in_tup[0].find("_")
    #run_label = (in_tup[0])[:sep0] # not used atm
    sensor_label = (in_tup[0])[sep1+1:]
    # set the dividor
    if "UV" in sensor_label or "Cond" == sensor_label or "Flow" == sensor_label:
        sensor_div = 1000.0
    elif "Pressure" in sensor_label:
        sensor_div = 100.0
    else:
        sensor_div = 10.0
    print((" Extracting: {0}").format(in_tup[0]))
    with open(file_in, 'rb') as fo:
        fread = fo.read()
        for i in range(start, end, 8):
            sread = struct.unpack("ii", fread[i:i+8])
            data=round((sread[0]/100.0)-inj_sel,4),sread[1]/sensor_div
            tmp_list0.append(data)
        tmp_list0.insert(0, in_tup[0])
    return(tmp_list0)


def inject_det(inp,show="false"):
    '''
    Finds injection points - required for adjusting retention
    '''
    injection_points = []
    header = readheader(inp)
    for i in header:
        if i[0] == Inject_id:
            injection = meta1_read(i[1:],silent="true")[1][0]
            injection_points.append(injection)
            if injection != 0.0:
                injection_points.insert(0,0.0)
    if injection_points == []:
        injection_points = [0.0]
        if show == "true":
            print((" ---- \n Injection points: \n # \t ml \n 0 \t {0}").format(injection_points[0]))
        return(injection_points)
    else:
        if show == "true":
            print(" ---- \n Injection points: \n # \t ml")
            for x,y in enumerate(injection_points):
                print((" {0} \t {1}").format(x,y))
        return(injection_points)


def store_in_list(inp):
    '''
    extract all data and store in list
    '''
    global inj_sel
    inj_sel = inject_det(file_in,show="false")[args.inject]
    header = readheader(inp)
    for i in header:
        x = dataextractor(i)
        if x == None:
            pass
        else:
            data_storage.append(x)


def writer(in_list):
    '''
    general writer, picks the correct writer function based on data type
    '''
    if in_list[0] == Methods_id or in_list[0] == CNotes_id:
        meta_writer(in_list[1:])
    else:
        data_writer(in_list)


def meta_writer(in_list):
    '''
    writes meta-data to txt-files
    '''
    ext = "_"+in_list[0]+".txt"
    inp_data = codecs.decode(in_list[1], 'cp1250')#.replace("\r","")
    content = inp_data.encode('utf-8')
    with open(file_base+ext,'wb') as fout:
        print(" Writing %s" %in_list[0])
        fout.write(content)

             
def data_writer(in_list):
    '''
    writes sensor/run-data to csv-files
    '''
    sep = in_list[0].find(":")+1
    label = in_list[0][sep:]
    run_name = in_list[0][:sep-1]
    fname = file_base+"_"+run_name+"_"+label+".csv"
    if "Logbook" in label:
        pass
    else:
        with open(fname, 'wb') as fout:
            print((" Writing {0}").format(label))
            for i in in_list[1:]:
                dp = (str(i[0]) + "," + str(i[1]) + str("\n")).encode('utf-8')
                fout.write(dp)


def plotter(in_list):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(" Matplotlib not found - Plotting will not work!")
        endscript()
    else:
        print(" Matplotlib found - WARNING - PLOTTING IS EXPERIMENTAL")
        sep = in_list[0].find(":")+1
        label = in_list[0][sep:]
        run_name = in_list[0][:sep-1]
        f_list = in_list[1:] #filter the list
        x_val=[x[0] for x in f_list]
        y_val=[x[1] for x in f_list]
        plot_x_min = in_list[1][0]
        plot_x_max = in_list[-1][0]
        if type(y_val[0]) == float:
            print(" Plotting %s" %label)
            plt.xlim(xmin = plot_x_min, xmax = plot_x_max)
            plt.title("%s: " %file_in + run_name + " - " "%s" %label)
            plt.xlabel('ml')
            plt.grid(which='major', axis='both')
            plt.plot(x_val, y_val,linewidth=1.5, alpha=0.85,color='b')
            ext = '.pdf'
            fname = file_base + "_Plot_" + run_name + "_" + label + ext
            plt.savefig(fname, dpi=300)
            plt.clf()


def endscript():
    print("stopping")
    exit

#--- run time---#

if args.user:
    showuser(file_in)
if args.check:
    input_check(file_in)
if args.info:
    showheader(file_in,full="false")
if args.points:
    inject_det(file_in,show="true")
if args.extract:
    data_storage = []
    if input_check(file_in) != (0,0,0):
        print("\n Gotta stop!")
    else:
        print("\n Go go go!")
        store_in_list(file_in)
        for i in data_storage:
            writer(i)
            if args.plot:
                plotter(i)

end = time.time()
runtime = str(end - start)
print("\n ----")
print((" Runtime: {0} seconds!").format(runtime[:5]))
