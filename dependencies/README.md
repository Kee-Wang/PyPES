#Dependency requirement
####### Dependency requirement #######
#Those are the dependencies that reqired to fully functionalize PyPES.

###### Requirements without Version Specifiers ######
get-pip.py
numpy
matplotlib

###### Requirements with Version Specifiers ######
#For python download: https://www.python.org/downloads/
python >= 3             # Recomended version python 3.5.2

###### Recommended user interface ######
ipython

##### Module Installation Instruction ######

1. First install get-pip.py

	`python get-pip.py`

2. Install PyPES module in the root directory

	`pip install -e .` (There is a dot in the end)

	`-e` here means editor mode. All your changes in the python source script will be reflected instantly.

	New packages are install in `/usr/local/lib/python3.X/site-packages/[PACKAGE].egg`

3. To install dependencies:

	`pip install <dependency_name>`

4. Then you can import configs anywhere from python or ipython

	`from pypes.configs import configs`

