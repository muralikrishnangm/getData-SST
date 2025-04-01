"""
Functions and class for data extraction.
"""

# Import libraries
import struct
import numpy as np
import time

def check_data(loadpath, nx, ny, nz, nbyte):
    """
    Function to check the data specifications with the content 
    of the binary file.
  
    Args:
        loadpath: Path to raw data
        nx: # grid points of raw data in x-direction
        ny: # grid points of raw data in y-direction
        nz: # grid points of raw data in z-direction
        nbyte: # bytes per data in file. 4 bytes = 8 hex characters = 32-bit little endian
    """
    print('Checking data file...')
    # read in test binary and check number of samples
    binary = open(loadpath, 'rb')
    binary.seek(0,2) ## seeks to the end of the file (needed for getting number of bytes)
    num_bytes = binary.tell() ## how many bytes are in this file is stored as num_bytes
    print(f'Number of bytes in file =\t{num_bytes:,}')
    print(f'Number of counted samples =\t{int(num_bytes/nbyte):,}')
    print(f'Number of actual samples =\t{nx*ny*nz:,}')
    if int(num_bytes/nbyte)==nx*ny*nz:
        num_samp = nx*ny*nz
        print(f'Number of samples counted == actual. Check complete.')
    else:
        raise Exception(f'Number of samples counted != actual')
    binary.close()

def get_1Dgrid(Lh, nx, nxoffset, nxsl, nxskip):
    """
    Function to compute grid coordinates for subdomain/box.
    
    Args:
        Lh: Physical length of 1D grid
        nx: # grid points of raw data in x-direction 
            Can be any direction based on the 1D grid of interest
        nxoffset: Offset these many samples in x dir to set corner of the sampled domain
        nxsl: # grid points in x dir of sampled data
        nxskip: Subsampling rate used in x dir

    Returns:
        x: 1D grid (location of grid points) for given grid resolution specifications
    """
    dx = Lh/nx
    xin = 0 + (dx*nxoffset)
    xfi = xin + dx*nxsl*nxskip
    x = np.linspace(xin, xfi, nxsl)
    return x

def get_data(method, loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip, nbyte):
    """
    Function to extract subdomain using a specified method.
    
    Args:
        method: Data loading method
        loadpath: Path to raw data
        nx: # grid points of raw data in x-direction 
        ny: # grid points of raw data in y-direction 
        nz: # grid points of raw data in z-direction 
        nxsl: # grid points in x dir of sampled data
        nysl: # grid points in y dir of sampled data
        nzsl: # grid points in z dir of sampled data
        nxoffset: Offset these many samples in x dir to set corner of the sampled domain
        nyoffset: Offset these many samples in y dir to set corner of the sampled domain
        nzoffset: Offset these many samples in z dir to set corner of the sampled domain
        nxskip: Subsampling rate used in x dir
        nyskip: Subsampling rate used in y dir
        nzskip: Subsampling rate used in z dir
        nbyte: # bytes per data in file. 4 bytes = 8 hex characters = 32-bit little endian

    Returns:
        datacube: 3D subdomain from raw data for given grid resolution specifications
    """
    if method == "seek":
        # using seek-type operation
        datacube = get_data_single_points(loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip, nbyte)
    elif method == "memmap":
        # using memmap
        datacube = get_data_memmap(loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip)
    elif method == "fromfile":
        # using fromfile loader
        datacube = get_data_fromfile(loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip, nbyte)
    else:
        raise NotImplementedError
    return datacube
   
# Different functions to get data from binary file
def get_data_single_points(loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip, nbyte):
    """
    Get data using Python's `open` file function and `seek` operation.
    NOTE: A single item is read at a time using this procedure. Not
    efficient to extract large subcubes.
    
    Args:
        loadpath: Path to raw data
        nx: # grid points of raw data in x-direction 
        ny: # grid points of raw data in y-direction 
        nz: # grid points of raw data in z-direction 
        nxsl: # grid points in x dir of sampled data
        nysl: # grid points in y dir of sampled data
        nzsl: # grid points in z dir of sampled data
        nxoffset: Offset these many samples in x dir to set corner of the sampled domain
        nyoffset: Offset these many samples in y dir to set corner of the sampled domain
        nzoffset: Offset these many samples in z dir to set corner of the sampled domain
        nxskip: Subsampling rate used in x dir
        nyskip: Subsampling rate used in y dir
        nzskip: Subsampling rate used in z dir
        nbyte: # bytes per data in file. 4 bytes = 8 hex characters = 32-bit little endian

    Returns:
        datacube: 3D subdomain from raw data for given grid resolution specifications
    """
    # skip these many samples in each direction
    nxoff = nxoffset * (nbyte)
    nyoff = nyoffset * (nx*nbyte)
    nzoff = nzoffset * (nx*ny*nbyte)
    # initial corner of the cuboid
    init = nzoff + nyoff + nxoff
    m = 0 # for sample counting
    datacube = np.zeros(nxsl*nysl*nzsl)
    t = time.time()
    binary = open(loadpath, 'rb')
    nyshift = 0   # to shift a slice of nx*ny
    for k in range(0,nzsl):
        nxshift = 0  # to shift a row of nx. Reset to 0 after 1 slice of nx*ny
        for j in range(0,nysl):
            for i in range(0,nxsl):
                pos = init + ((i*nxskip)*nbyte) + nxshift + nyshift
                binary.seek(pos,0)
                datacube[m] = struct.unpack("<f", binary.read(nbyte))[0]
                m += 1
            nxshift += (nx*nyskip)*nbyte  # shift 1 row of nx for every ny
        nyshift += ((nx*ny)*nzskip)*nbyte  # shift 1 slice of nx*ny for every nz
    binary.close()
    elpsdt = time.time() - t
    print(f'Time elapsed for loading datacube: {int(elpsdt/60)} min {elpsdt%60:.2f} sec')
    # reshape 1-D to 3-D array
    datacube = datacube.reshape(nzsl,nysl,nxsl)
    # transpose x & z dimensions
    datacube = datacube.transpose(2,1,0)
    # Print the shape of the sub-cube
    print(f'Shape of the sub-cube: {datacube.shape}')
    return datacube

def get_data_memmap(loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip):
    """
    Get data using numpy memmap function.
    NOTE: memmap outputs a file pointer accessing the binary data 
    in the shape provided. Need to copy data to local memory.
    
    Args:
        loadpath: Path to raw data
        nx: # grid points of raw data in x-direction 
        ny: # grid points of raw data in y-direction 
        nz: # grid points of raw data in z-direction 
        nxsl: # grid points in x dir of sampled data
        nysl: # grid points in y dir of sampled data
        nzsl: # grid points in z dir of sampled data
        nxoffset: Offset these many samples in x dir to set corner of the sampled domain
        nyoffset: Offset these many samples in y dir to set corner of the sampled domain
        nzoffset: Offset these many samples in z dir to set corner of the sampled domain
        nxskip: Subsampling rate used in x dir
        nyskip: Subsampling rate used in y dir
        nzskip: Subsampling rate used in z dir
        nbyte: # bytes per data in file. 4 bytes = 8 hex characters = 32-bit little endian

    Returns:
        datacube: 3D subdomain from raw data for given grid resolution specifications
    """
    # Memory-map the binary file
    t = time.time()
    data_memmap = np.memmap(loadpath, dtype=np.float32, mode='r', shape=(nz, ny, nx)) # NOTE: data is stored [z, y, x]
    elpsdt = time.time() - t
    print(f'Time elapsed for memmap: {int(elpsdt/60)} min {elpsdt%60:.4f} sec')
    # Extract the sub-cube
    t = time.time()
    sub_cube = data_memmap[ nzoffset:nzoffset+(nzsl*nzskip):nzskip, # start from `nzoffset` location and get `nzsl` points, but skip every `nzskip` point
                          nyoffset:nyoffset+(nysl*nyskip):nyskip, 
                          nxoffset:nxoffset+(nxsl*nxskip):nxskip] 
    elpsdt = time.time() - t
    print(f'Time elapsed for slice: {int(elpsdt/60)} min {elpsdt%60:.4f} sec')
    # Copy the sub-cube to a new array to avoid memory-mapping issues when processing
    t = time.time()
    datacube = sub_cube.copy().transpose(2,1,0) # transposing data to be [x, y, z]
    elpsdt = time.time() - t
    print(f'Time elapsed for copying data: {int(elpsdt/60)} min {elpsdt%60:.4f} sec')
    data_memmap._mmap.close()
    del data_memmap, sub_cube
    # Print the shape of the sub-cube
    print(f'Shape of the sub-cube: {datacube.shape}')
    return datacube

def get_data_fromfile(loadpath, nx, ny, nz, nxsl, nysl, nzsl, nxoffset, nyoffset, nzoffset, nxskip, nyskip, nzskip, nbyte):
    """
    Get data using numpy fromfile function. 
    NOTE: fromfile reads the `count` # of items (of type `dtype`) 
    in a continous manner of memory storage.
    
    Args:
        loadpath: Path to raw data
        nx: # grid points of raw data in x-direction 
        ny: # grid points of raw data in y-direction 
        nz: # grid points of raw data in z-direction 
        nxsl: # grid points in x dir of sampled data
        nysl: # grid points in y dir of sampled data
        nzsl: # grid points in z dir of sampled data
        nxoffset: Offset these many samples in x dir to set corner of the sampled domain
        nyoffset: Offset these many samples in y dir to set corner of the sampled domain
        nzoffset: Offset these many samples in z dir to set corner of the sampled domain
        nxskip: Subsampling rate used in x dir
        nyskip: Subsampling rate used in y dir
        nzskip: Subsampling rate used in z dir
        nbyte: # bytes per data in file. 4 bytes = 8 hex characters = 32-bit little endian

    Returns:
        datacube: 3D subdomain from raw data for given grid resolution specifications
    """
    # skip these many samples in each direction
    nxoff = nxoffset * (nbyte)
    nyoff = nyoffset * (nx*nbyte)
    nzoff = nzoffset * (nx*ny*nbyte)
    # initial corner of the cuboid
    init = nzoff + nyoff + nxoff
    datacube = np.zeros((nxsl, nysl, nzsl))
    t = time.time()
    nyshift = 0   # to shift a slice of nx*ny
    for k in range(0,nzsl):
        nxshift = 0  # to shift a row of nx. Reset to 0 after 1 slice of nx*ny
        for j in range(0,nysl):
            pos = init + nxshift + nyshift
            datacube[:, j, k] = np.fromfile(loadpath, dtype=np.float32, count=nxsl*nxskip, offset=pos)[::nxskip]
            nxshift += (nx*nyskip)*nbyte  # shift 1 row of nx for every ny
        nyshift += ((nx*ny)*nzskip)*nbyte  # shift 1 slice of nx*ny for every nz
    elpsdt = time.time() - t
    print(f'Time elapsed for loading datacube: {int(elpsdt/60)} min {elpsdt%60:.2f} sec')
    print(f'Shape of the sub-cube: {datacube.shape}')
    return datacube

# Data Class using functions from above
class SSTData():
    def __init__(self, args):
        """
        Data class to store specification and extract sub-domain 
        using functions from above.
        
        Args:
            args: Argument parser variable with all the specifications
        """
        self.var = args.var
        self.nx = args.nx
        self.ny = args.ny
        self.nz = args.nz
        self.nbyte = 4 # 4 bytes = 8 hex characters = 32-bit little endian
        
        # domain size
        self.Lh = args.Lh 
        self.Lv = args.Lv
        
        # Sample cuboid
        # number of samples in each direction
        self.nxsl = args.nxsl
        self.nysl = args.nysl
        self.nzsl = args.nzsl
        # subsampling rate (for full resolution, use value 1)
        self.nxskip = args.nxskip
        self.nyskip = args.nyskip
        self.nzskip = args.nzskip
        # offset these many samples in each direction
        self.nxoffset = args.nxoffset
        self.nyoffset = args.nyoffset
        self.nzoffset = args.nzoffset
        
        self.loadpath = f'{args.rawdatadir}{self.var}{args.delimiter}{args.rawdatasnap}'
        self.method = args.method
        
    # Check data file
    def _check_data(self):
        check_data(self.loadpath, self.nx, self.ny, self.nz, self.nbyte)

    # Extract data
    def _get_data(self):
        datacube = get_data(self.method, self.loadpath, self.nx, self.ny, self.nz, 
                            self.nxsl, self.nysl, self.nzsl, self.nxoffset, self.nyoffset, self.nzoffset, 
                            self.nxskip, self.nyskip, self.nzskip, self.nbyte)
        return datacube
      
