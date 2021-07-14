# Hello `cathy`!

This is a command line application for manipulating catheter data from the Wright Group HeartVista catheter tracking sequences. This has been tested under Python 3.6.

## Installing
Get this project from the Wright Group gitlab using "git clone *address*", using the address listed on the project page. For example, if you are on the Sunnybrook network and have ssh keys set up:
`git clone git@panoptes.sri.utoronto.ca:wright-group/cathy.git`

### Virtual Environments
We recommend using [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [venv](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) to isolate your Python development environments. Examples below will use conda 4.6+ (activating environments prior to conda 4.6 varies according to OS, please see the note in the intro of [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)).

You can start by creating a conda environment with the tested version of Python. We've called it cathyEnv below, but you can name it as you wish:
`conda create -n cathyEnv python=3.6`

Then activate the environment:
`conda activate cathyEnv`

### Installing special dependencies
cathy depends on the following packages which are currently only available on the Wright Group gitlab on panoptes. If you are not on the internal network, you will need to [set up remote access](https://wrightgroup.sri.utoronto.ca/tiki-download_file.php?fileId=278) and replace "panoptes.sri.utoronto.ca" in the links below with "localhost":
- [catheter_utils](http://panoptes.sri.utoronto.ca:8088/wright-group/catheter_utils)
- [catheter_ukf](http://panoptes.sri.utoronto.ca:8088/wright-group/catheter_ukf)
- [dicom_utils](http://panoptes.sri.utoronto.ca:8088/wright-group/dicom_utils)
- [dicom_art](http://panoptes.sri.utoronto.ca:8088/wright-group/dicom_art)
- [get_gt](http://panoptes.sri.utoronto.ca:8088/wright-group/CLI_DataAugment) (from the CLI_DataAugment project)

Use "git clone" to get each of the above projects in their own directory.   Then, in each project directory, use "pip" to install it into the current environment. Install them in the order given above for dependencies to be satisfied:
`pip install -e .`

### Installing cathy
Now that you have installed the custom packages, you should be able to install cathy into the same environment. From the cathy project directory (the top-level directory containing this README.md file) run:
`pip install -e .`

## Commands

Documentation about specific parameters and options is available for each specific command: for example, run `cathy info --help` to find out more about the info command.

*Reminder*: You will need to activate the environment in which you installed cathy in order to run it. If you activate a different environment, exit the terminal, log out of or restart your computer, etc, you will no longer be in cathy's environment.

### `cathy info`:
Get information about a .projection file, or a directory containing related projections files.

### `cathy peek`:
Visualize catheter projections with optional coil location data.

### `cathy localize`:
Apply localization algorithms to process projections files into catheter coordinates files.

### `cathy apply-ukf`:
Smooth and fit catheter coordinate data using an unscented Kalman filter developed for this purpose.

### `cathy scatter`:
Visualize catheter coordinates on a dicom image.

### `cathy gt-tool`:
Select object on dicom images to export ground truth postions.

### `cathy view-dicom`:
View dicom image.

## Resources

### Click
Click is the library used to generate the command line interface for cathy.
It generates a CLI from a tree of functions with appropriate decorators. Documentation and extensive examples can be found here:
https://click.palletsprojects.com/en/7.x/.