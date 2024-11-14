# ZPATH

This git page has been created to gather "hands-on" material to be used in teaching. Se more details below.

To copy all the material on your PC, in the terminal do `mkdir zpath; git clone https://github.uio.no/zpath/software.git`

## Installation instructions 

To get all the software needed to run the programs the easiest is to install anaconda as described below. 

### First time setup (LINUX)

For other OSs see: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation

The following only needs to be done once:

```
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
bash Anaconda3-2020.11-Linux-x86_64.sh
```

Follow the instructions on the screen. 

**Recomandation**: When the installer asks `Do you wish the installer to initialize Anaconda3` answer `no` to avoid the conda environment to be enabled by default whenever you start a new shell. 

More details on installation can be found here:

**Linux**: https://docs.anaconda.com/anaconda/install/linux/

Activate anaconda in your terminal

`source <path-to-where-anaconda-is-installed>/etc/profile.d/conda.sh`

Build the environment (using the python3 yml file in the github repo) - this may take some time

`conda env create -f environment_fys5555_py3.yml`

Then load the environement, and your good to go!

`conda activate fys5555_py3`

### Whenever you start a new session

If you have installed conda as described above, when you start a new shell, all you have to do is 

```
source <path-to-where-anaconda-is-installed>/etc/profile.d/conda.sh
conda activate fys5555_py3
``` 
