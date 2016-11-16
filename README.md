# PyPES

PyPES is a Python module that will be develped into a Python-based independend softwere that will help making potential energy surface (PES).

This package is a product of Bowman Group and is curretly developed by Kee.

## Dependency requirement


Those are the dependencies that reqired to fully functionalize PyPES.

##Requirements without Version Specifiers ######
get-pip.py

numpy

matplotlib

## Requirements with Version Specifiers ######
For python download: https://www.python.org/downloads/
python >= 3             # Recomended version python 3.5.2

## Recommended user interface ######
ipython

# Module Installation Instruction ######
##All dependencies will be installed automatically without requiring root permission by installing Anaconda3. Strongly suggest using Anaconda to install all dependencies.

1. First install anaconda3 (meaning with Python 3 distribution)

        `bash Anaconda3-4.2.0-Linux-x86_64.sh`

2. Install PyPES module in the root directory

        `pip install -e .` (There is a dot in the end)

        `-e` here means editor mode. All your changes in the python source script will be reflected instantly.

        New packages are install in `/usr/local/lib/python3.X/site-packages/[PACKAGE].egg`


3. Then you can import configs anywhere from python or ipython

        `from pypes.configs import configs`

4. To uninstall:

        `pip uninstall pypes`

