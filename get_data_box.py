"""
Script to extract sub-domains from full 3D flow field of
the stably stratified turbulence database. A sample plotting
routine is also provided to visualize the extracted sub-domain.

The filename is of format [var]_[time-step]
The time-step is denoted with 6 digits of precision.
E.g.: for velocity in x-direction at time-step 10.15, the
filename will be `u_10.150000`

The data is stored in binary files with 32-bit, little-endian.
The dimension of the stored data is in the format [(nx+2) x ny x nz], 
with two zeros padded at the end of the first (x) direction.
So, in the `nx` argument, add 2 along with true nx value of
the dataset. 
E.g.: if the dataset has nx=512, use `--nx 514`

In the `method` argument, `memmap` is the most efficient one 
for extracting large sub-domains. It uses memory-mapped files 
for accessing small segments of large files on disk, 
without reading the entire file into memory.

Sample usage: 
python get_data_box.py --rawdatadir ./rawData/ --rawdatasnap 10.000000 --var u --nx 514 --ny 256 --nz 512 --gravity y --Lv 0.5 --method memmap --nxsl 64 --nysl 64 --nzsl 64 --nxoffset 32 --nyoffset 32 --nzoffset 32 --nxskip 2 --nyskip 2 --nzskip 2 --saveFig --saveData --saveDataDir ./saveData/

Authors: 
Murali Gopalakrishnan Meena (Oak Ridge National Laboratory)
Steve de Bruyn Kops (University of Massachusetts Amherst)
"""

# Import libraries
import numpy as np
from funcs_data import SSTData, get_1Dgrid
from funcs_plotting import plot_contour_box
import matplotlib.pyplot as plt
import argparse

# Args
parser = argparse.ArgumentParser()
# Raw data specifications
parser.add_argument("--var", required=True, help="Variable name: 'r' 'u' 'v' 'w'")
parser.add_argument("--rawdatadir", type=str, required=True, help="Dir of raw data")
parser.add_argument("--rawdatasnap", type=str, required=True, help="Snapshot of raw data")
parser.add_argument("--delimiter", type=str, required=True, help="Delimiter to seperate var and snapshot. E.g. use `_` for `rho_15000`")
parser.add_argument("--nx", type=int, required=True, help="Number of grid points in x dir for full data. Need to add 2: nx+2")
parser.add_argument("--ny", type=int, required=True, help="Number of grid points in y dir for full data")
parser.add_argument("--nz", type=int, required=True, help="Number of grid points in z dir for full data")
parser.add_argument("--drhobar", type=float, default=-577.31, required=False, help="drhobar for computing ambient density. See global.hpp")
parser.add_argument("--gravity", type=str, default="y", required=False, help="Direction of gravity: 'y' or 'z'. Can be understood based on ny and nz.")
parser.add_argument("--Lh", type=float, default=1.0, required=False, help="Horizontal length of full box. See global.hpp")
parser.add_argument("--Lv", type=float, default=0.5, required=False, help="Vertical length of full box. See global.hpp")
# Method to extract sub-domain
data_methods = ['seek', 'memmap', 'fromfile']
parser.add_argument("--method", type=str, default="memmap", required=False, choices=data_methods, help="Data loading method")
# Subcube specifications
parser.add_argument("--nxsl", type=int, default=32, required=False, help="Number of grid points in x dir for sampled data")
parser.add_argument("--nysl", type=int, default=32, required=False, help="Number of grid points in y dir for sampled data")
parser.add_argument("--nzsl", type=int, default=32, required=False, help="Number of grid points in z dir for sampled data")
parser.add_argument("--nxskip", type=int, default=1, required=False, help="Subsampling rate in x dir (for full resolution, use value 1)")
parser.add_argument("--nyskip", type=int, default=1, required=False, help="Subsampling rate in y dir (for full resolution, use value 1)")
parser.add_argument("--nzskip", type=int, default=1, required=False, help="Subsampling rate in z dir (for full resolution, use value 1)")
parser.add_argument("--nxoffset", type=int, default=0, required=False, help="Offset these many samples in each direction in x dir to set corner of the sampled box")
parser.add_argument("--nyoffset", type=int, default=0, required=False, help="Offset these many samples in each direction in x dir to set corner of the sampled box")
parser.add_argument("--nzoffset", type=int, default=0, required=False, help="Offset these many samples in each direction in x dir to set corner of the sampled box")
# Plotting and saving data 
parser.add_argument("--saveFig", default=False, action='store_true', help="Save figure")
parser.add_argument("--saveData", default=False, action='store_true', help="Save data")
parser.add_argument("--saveDataDir", type=str, default="./saveData/", required=False, help="Dir to save data")
args = parser.parse_args()

if __name__ == '__main__':
    # ===============
    # Data info
    # save fig & data
    saveFig = args.saveFig
    figname = f'{args.var}_PVsample-nx{int(args.nxsl)}ny{int(args.nysl)}nz{int(args.nzsl)}_nskip{args.nxskip}_t{args.rawdatasnap}_{args.method}'
    dpi_val = 200
    saveData = args.saveData

    # ===============
    # Get data
    print(f'Loading snapshot: {args.var} {args.rawdatasnap}')
    data_sample = SSTData(args)
    # check-sum for data file 
    data_sample._check_data()
    # get subdomain
    datacube = data_sample._get_data()

    # ===============
    # Grid
    x = get_1Dgrid(data_sample.Lh, data_sample.nx, data_sample.nxoffset, data_sample.nxsl, data_sample.nxskip)
    if args.gravity == 'y':
        y = get_1Dgrid(data_sample.Lv, data_sample.ny, data_sample.nyoffset, data_sample.nysl, data_sample.nyskip)
        z = get_1Dgrid(data_sample.Lh, data_sample.nz, data_sample.nzoffset, data_sample.nzsl, data_sample.nzskip)
        if args.var == 'r' or args.var == 'rho':
            # add ambient density (y * drhobardz) to perturbation density (r)
            scalar = y * args.drhobar
            datacube = datacube * scalar.reshape(1, args.nysl, 1)
    elif args.gravity == 'z':
        y = get_1Dgrid(data_sample.Lh, data_sample.ny, data_sample.nyoffset, data_sample.nysl, data_sample.nyskip)
        z = get_1Dgrid(data_sample.Lv, data_sample.nz, data_sample.nzoffset, data_sample.nzsl, data_sample.nzskip)
        if args.var == 'r' or args.var == 'rho':
            # add ambient density (y * drhobardz) to perturbation density (r)
            scalar = z * args.drhobar
            datacube = datacube * scalar.reshape(1, 1, args.nzsl)
    else:
        raise Exception("Gravity should be defined")
    print(f"data shape: {datacube.shape}")

    # ===============
    # Plotting

    # Rotated box
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    plot_contour_box(ax, x, y, z, datacube, args.gravity)
    if saveFig: 
        plt.savefig(f'Figs/{figname}_box_plot.png', bbox_inches='tight', dpi=dpi_val)


    # ===============
    # # Save data
    if saveData: 
        np.savez(f'{args.saveDataDir}{figname}.npz', datacube=datacube, x=x, y=y, z=z, args=args)

