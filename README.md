# Get SST Data

Scripts to access and process data from the stably stratified turbulence (SST) database.

# Installation

1. Create custom conda environment.
    <details><summary>Notes for OLCF:</summary>

      * Make sure to install your custom env in either `/ccs` or `/lustre/orion` (recommended). This is required to seamlessly run the plotting routines on OLCF JupyterLab.
      * Follow steps in [OLCF Docs](https://docs.olcf.ornl.gov/software/python/index.html) for loading base conda env or the following for your own mini-conda installation

      ```
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      bash Miniconda3-latest-Linux-x86_64.sh -b -p /lustre/orion/proj-shared/[project-id]/[dir-to-install-mini-conda]
      source /lustre/orion/proj-shared/[project-id]/[dir-to-install-mini-conda]/bin/activate
      ```
      
    </details>

    * Make custom conda env
      ```
      conda create --name env_getdata-sst python=3.11
      conda activate env_getdata-sst
      ```
2. Install libraries
    ```
    pip install -r requirements.txt --no-cache-dir
    ```

# Usage

* [get_data_box.py](get_data_box.py) - Extract a sub-domain from the full data and visualize it.
```
python get_data_box.py --rawdatadir [data-dir] --rawdatasnap [snapshot-time-step] --var r --nx 514 --ny 256 --nz 512 --saveFig
```

# Authors
* Murali Gopalakrishnan Meena (Oak Ridge National Laboratory)
* Steve de Bruyn Kops (University of Massachusetts Amherst)

# Cite this code and database

* [SST database](https://stratified-turbulence.github.io/web/)

**Data References:**
* Riley, J. J., & De Bruyn Kops, S. M. (2003). [Dynamics of turbulence strongly influenced by buoyancy. Physics of Fluids, 15(7), 2047-2059.](https://doi.org/10.1063/1.1578077)
* de Bruyn Kops, S. M., & Riley, J. J. (2019). [The effects of stable stratification on the decay of initially isotropic homogeneous turbulence. Journal of Fluid Mechanics, 860, 787-821.](https://doi.org/10.1017/jfm.2018.888)
* Riley, J. J., Couchman, M. M., & de Bruyn Kops, S. M. (2023). [The effect of Prandtl number on decaying stratified turbulence. Journal of Turbulence, 24(6-7), 330-348.](https://doi.org/10.1080/14685248.2023.2178654)
