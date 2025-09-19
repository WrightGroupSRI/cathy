# Hello `cathy`!

This is a command line application for manipulating catheter MR data recorded from the Wright Group HeartVista catheter tracking sequences. This has been tested under Python 3.6 and 3.8 on Linux, and Python 3.11 on Windows.

## Windows Setup
The development and testing of cathy has been done on Linux systems and most of the commands outside this section assume a bash terminal. If you are on Windows, the WSL (Windows Subsystem for Linux) can be used: this has not been tested with cathy, but should work similarly to a regular Linux system.

Alternatively, if you want to use the Windows versions of Python and Anaconda, this section describes how to set up a git-bash terminal to use these. If you haven't already, install git on Windows following the [official instructions](https://git-scm.com/download/win); this will install git-bash as well as git.

### Installing Python and Anaconda on Windows
This setup was tested using the default Windows downloads from the official Python and Anaconda sites:
- Python 3.11 was installed from the [python website](https://www.python.org/downloads/)
- Anaconda version 2023.07-1 was installed from the [Anaconda website](https://www.anaconda.com/download)

### Set up paths and aliases
In the git-bash terminal:
1. Get the path to anaconda and python: `CDIR=$(dirname $(find {HARDDRIVE_PATH} -name "_conda.exe" 2>/dev/null))`
    - NOTE: Replace {HARDDRIVE_PATH} with your hard drive path that has Anaconda installed. Ex. /c (if C: Drive)
2. echo $CDIR
    - NOTE: if there is more than one path listed then redefine variable CDIR with the correct path. 
    - `CDIR={PATH_TO_ANACONDA}` (replace {PATH_TO_ANACONDA} with the correct path)
3. Update the path for bash: `echo 'export PATH="$PATH:$CDIR:$CDIR/Scripts"' >> ~/.bashrc`
4. Create an alias for Python: `echo 'alias python="winpty python.exe"'>> ~/.bashrc`
5. Source the bash file: `source ~/.bashrc`
6. Test by checking the versions of python and conda: `python --version` and `conda --version`
    - if running python from git-bash does not work with an "Access is denied error", refer to the next section

#### Solving Python Access Issue
1. Launch Windows Settings and navigate to "App Execution Aliases"
2. Turn off "python.exe" and "python3.exe" App Installers

## Installing
Get this project from github using "git clone" or "gh repo clone", using the address from the "Code" button dropdown with the protocol of your choice.

### Option 1: Build a docker image
From the code directory:
`docker build -t cathy-mm .`

### Option 2: Use a Virtual Environment
You can use [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [venv](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) to isolate your Python development environments. Examples below will use conda 4.6+ (activating environments prior to conda 4.6 varies according to OS, please see the note in the intro of [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)).

You can start by creating a conda environment with the tested version of Python. We've called it cathyEnv below, but you can name it as you wish:
`conda create -n cathyEnv python=3.6`
   - Note: Minimum requirement for Windows is Python 3.8. `conda create -n cathyEnv python=3.8`

Then activate the environment:
`conda activate cathyEnv`

Subsequent package installations with pip and conda, as detailed in the next section, should be done with this new environment activated.

#### Installing special dependencies

cathy depends on the following packages:
- [catheter_utils](https://github.com/WrightGroupSRI/catheter_utils)
- [catheter_ukf](https://github.com/WrightGroupSRI/catheter_ukf)
- [dicom_utils](https://github.com/WrightGroupSRI/dicom_utils)
- [dicom_art](https://github.com/WrightGroupSRI/dicom_art)
- [get_gt](https://github.com/WrightGroupSRI/get_gt)

For each of the dependencies above **in the order shown**, from your terminal:
1. Perform the git clone. The git clone link is accessible from the "Code" button dropdown for each project page.
2. A new project directory will be created after the previous command. Run `pip install -e [project_directory]` to install the package.

##### ffmpeg
Additionally, ffmpeg needs to be installed for some features. This can be done in either of the following ways:

- To install ffmpeg in a system-level location, you can use your package manager on Linux or follow the instructions on the [official site](https://ffmpeg.org/download.html). Admin privileges on your machine will be required.
- To install it to a local user directory using conda, use:
`conda install conda-forge::ffmpeg`

Note: this has only been tested on Linux

#### Installing cathy
Now that you have installed the custom packages, you should be able to install cathy into the same environment. From the cathy project directory (the top-level directory containing this README.md file) run:
`pip install -e .`

## Commands
If you are using docker, precede the commands with "docker run --rm" and give permissions as necessary for file access (through bind mounts) and/or your display for graphical output.

Documentation about specific parameters and options is available for each specific command: for example, run `cathy info --help` to find out more about the info command.

*Reminder*: If you are using a virtual environment, you will need to activate the environment in which you installed cathy in order to run it. If you activate a different environment, exit the terminal, log out of or restart your computer, etc, you will no longer be in cathy's environment.

### `cathy info`:
Get information about a .projection file, or a directory containing related projections files.

### `cathy peek`:
Visualize catheter projections with optional coil location data.

### `cathy localize`:
Apply localization algorithms to process projections files into catheter coordinates files.

### `cathy apply-ukf`:
Smooth and fit catheter coordinate data using an unscented Kalman filter developed for this purpose.

### `cathy fft-projections`:

Extract 1D signals from tracking sequence projection files and calculates magnitude of absolute fourier transform and shift of signal. Data outputted as text file and will be saved to the fft subfolder.
Textfile name template is "fft_signals_coilX_recX_fovXXX_readoutsXXX.txt" (coil index, recording index, field of view (mm), # of readouts).

Example Format of text file:

PRJTXT (file format name)

1 (version)
    
SRI_April-2020-07-09T12_20_56.765 (tracking sequence folder name)
    
x y z x_timestamp average_timestamp readout_size (column headings: axes, timestamps, ... length of each readout)


Note: Data organized as columns, to be processed vertically for each column heading. 

usage: cathy fft-projections PATH_TO_TRACKING_SEQUENCE_FOLDER DESTINATION_PATH RECORDING_INDEXES [OPTIONAL: -d/--distal_index, -p/--proximal_index, -z/--dither_index]  

Note: For RECORDING_INDEXES can specify range or individual comma-separated values (ex. 0-2 or 0,1,2)

TODO: select and manage multiple dither indexes/folders

### `cathy coil-metrics`:
Calculate or plot tracking error results of localization algorithms (require output directory from cathy localize).

Usage: cathy coil-metrics calc [OPTIONAL: --dest, -d/--distal_index, -p/--proximal_index, -z/--dither_index, --expname] path_to_tracking_sequence_folder path_to_localization_results_folder path_to_groundtruth_folder/ localization_algorithms

Sample Data can be found on ircci-share drive. PATH=ircci-share/ircci/Data/MR_Tracking/trackSpoilers-09July2020

Example: cathy coil-metrics calc --dest PATH/SRI_April-2020-07-09T12_20_56.765/results -d 7 -p 6 -z 0 --expname VaringSpoilers PATH/SRI_April-2020-07-09T12_20_56.765 PATH/SRI_April-2020-07-09T12_20_56.765/results PATH/SRI_April-2020-07-09T12_20_56.765/groundtruth/ peak,centroid_around_peak,png

Usage: cathy coil-metrics plotter [OPTIONS: --expname, --xaxis] path_to_folder_with_TrackingErr.csv  

Note: *TrackingErr.csv created after running cathy coil-metrics calc and located at --dest from cathy coil-metrics calc. If --dest not used located at path_to_localization_results_folder

### `cathy scatter`:
Visualize catheter coordinates on a dicom image.

### `cathy gt-tool`:
Select object on dicom images to export ground truth positions.

### `cathy view-dicom`:
View dicom image.

## Resources

### Click
Click is the library used to generate the command line interface for cathy.
It generates a CLI from a tree of functions with appropriate decorators. Documentation and extensive examples can be found here:
https://click.palletsprojects.com/en/7.x/.
