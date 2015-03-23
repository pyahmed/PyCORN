import sys
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

dr = 'PyCORN'
if dr not in sys.path: sys.path.append(dr)
from pycorn import Pycorn

blue = '#377EB8'
green = '#4DAF4A'
parser = argparse.ArgumentParser()

parser.add_argument("inp_res",
                    help = "Input .res file(s)",
                    nargs='+',
                    metavar = "<file>.res")
parser.add_argument("-e", "--ext", default='pdf',
                    help = "Image type to use, e.g. 'jpg', 'png', 'eps', or 'pdf' (default: pdf)")
parser.add_argument("--xmin", default=None, type=float,
                    help = "Lower bound on the x-axis")
parser.add_argument("--xmax", default=None, type=float,
                    help = "Upper bound on the x-axis")
parser.add_argument("--ymin", default=0, type=float,
                    help = "Lower bound on the y-axis")
parser.add_argument("--ymax", default=None, type=float,
                    help = "Upper bound on the y-axis")
parser.add_argument("--ymin2", default=None, type=float,
                    help = "Lower bound on the right y-axis")
parser.add_argument("--ymax2", default=None, type=float,
                    help = "Upper bound on the right y-axis")
parser.add_argument("--dpi", default=None, type=int,
                    help = "DPI (dots per inch) for raster images (png, jpg, etc.)")

args = parser.parse_args()

for fname in args.inp_res:
    path = Path(fname)
    pc = Pycorn(fname)
    pc.load()
    UVx, UV = np.array(pc['UV']['data']).T
    ix = np.isfinite(UVx)
    if args.xmin is not None:
        ix = ix & (UVx >= args.xmin)
    if args.xmax is not None:
        ix = ix & (UVx <= args.xmax)
    UV -= min(UV[ix])
    xmin = min(UVx) if args.xmin is None else args.xmin
    xmax = max(UVx) if args.xmax is None else args.xmax
    condx, cond = np.array(pc['Cond']['data']).T
    
    for name in 'general', 'sizeex':
        fig, ax = plt.subplots(figsize=(4,3))
        ax.plot(UVx, UV, color=blue, lw=2)

        ax2 = ax.twinx()
        ax2.plot(condx, cond, color=green, lw=2)

        ax.set_xlabel('Elution Volume (ml)')
        ax.set_ylabel('Absorbance (mAu)', color=blue)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(args.ymin, args.ymax)
        ax2.set_ylabel(
            r'Conductivity $\left(\frac{\mathrm{mS}}{\mathrm{cm}}\right)$',
            color=green)

        for tl in ax.get_yticklabels():
            tl.set_color(blue)

        for tl in ax2.get_yticklabels():
            tl.set_color(green)
        
        ymax2 = args.ymax2
        if name == 'sizeex':
            mx = max(cond)
            exp = int(np.log10(mx/2.))
            ymax2 = round(mx*2, -exp)*2
        ax2.set_ylim(args.ymin2, ymax2)
        
        outputpath = path.with_name('{}-{}.{}'.format(
                    path.stem, name, args.ext))
        plt.savefig(str(outputpath), bbox_inches='tight', dpi=args.dpi)
        print(outputpath)
