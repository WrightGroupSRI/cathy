# Windows Setup
The development and testing of cathy has been done on Linux systems and most of the commands outside this section assume a bash terminal. If you are on Windows, the WSL (Windows Subsystem for Linux) can be used: this has not been tested with cathy, but should work similarly to a regular Linux system.

Alternatively, if you want to use the Windows versions of Python and Anaconda, this section describes how to set up a git-bash terminal to use these. If you haven't already, install git on Windows following the [official instructions](https://git-scm.com/download/win); this will install git-bash as well as git.

## Installing Python and Anaconda on Windows
This setup was tested using the default Windows downloads from the official Python and Anaconda sites:
- Python 3.11 was installed from the [python website](https://www.python.org/downloads/)
- Anaconda version 2023.07-1 was installed from the [Anaconda website](https://www.anaconda.com/download)

## Set up paths and aliases
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