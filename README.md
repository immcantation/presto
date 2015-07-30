Version 0.5.0.beta:  July 29, 2015
===============================================================================

pRESTO is a toolkit for processing raw reads from high-throughput sequencing 
of lymphocyte repertoires. Dramatic improvements in high-throughput sequencing 
technologies now enable large-scale characterization of immunoglobulin 
repertoires, defined as the collection of trans-membrane antigen-receptor 
proteins located on the surface of T and B lymphocytes. The REpertoire 
Sequencing TOolkit (pRESTO) is composed of a suite of utilities to handle all 
stages of sequence processing prior to germline segment assignment. pRESTO is 
designed to handle either single reads or paired-end reads. It includes 
features for quality control, primer masking, annotation of reads with sequence 
embedded barcodes, generation of single-molecule consensus sequences, assembly 
of paired-end reads and identification of duplicate sequences. Numerous options 
for sequence sorting, sampling and conversion operations are also included.


Requirements
-------------------------------------------------------------------------------

Software        | Link
--------------- | -----------------------------
Python 3.4.0    | http://python.org
setuptools 2.0  | http://bitbucket.org/pypa/setuptools
NumPy 1.9       | http://numpy.org
SciPy 0.14      | http://scipy.org
pandas 0.16     | http://pandas.pydata.org
Biopython 1.65  | http://biopython.org
MUSCLE v3.8     | http://www.drive5.com/muscle
USEARCH v7.0    | http://www.drive5.com/usearch


Installation - Linux
-------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the full 
   SciPy stack using the instructions at http://scipy.org/install.html, then 
   install Biopython according to the 
   [instructions](http://biopython.org/DIST/docs/install/Installation.html).

2. Unzip the pRESTO bundle and run `python setup.py install --user`.


Installation - Windows
-------------------------------------------------------------------------------

1. Install Python 2.7.5+ from [Python](http://python.org/download).

2. Install pip and setuptools using these [installation instructions]
   (https://pip.pypa.io/en/latest/installing.html).

3. Install NumPy, SciPy, pandas and Biopython using the packages available 
   from the [Unofficial Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs)
   binary collection.

4. Unzip the pRESTO bundle, open a Command Prompt, and run
   `python setup.py install` from the pRESTO folder.
   
5. For a default installation of Python 2.7, the pRESTO scripts will be 
   installed into `C:\Python27\Scripts` and should be directly executable from 
   the Command Prompt. If this is not the case, then follow steps 5-6 below.
   
6. Add the installation directory to your `%Path%`. On Windows 7 the %Path% 
   setting is located under 'Control Panel' -> 'System and Security' -> 
   'System' -> 'Advanced System Settings' -> 'Environment variables' -> 
   'System variables' -> 'Path'.

7. Set the file association for Python ('.py') files by right-clicking on a 
   '.py' file, selecting 'Open with' -> 'Choose default program...', choosing the 
   'python.exe' executable from the Python 2.7 folder, and checking 
   'Always use the selected program'.
   

Installation - Mac OS X
-------------------------------------------------------------------------------

1. Install Xcode 3.2.6
   Available from the Apple store or 
   [developer downloads](http://developer.apple.com/downloads).
   If you have a newer version (eg, Xcode 4.6.3) that will work also,
   but Xcode 3 is free of charge.  If Xcode fails to install with an
   "Unknown Error", change the date on your system to some time in 2011,
   install Xcode, and then change the date back to the proper setting.

2. Install XQuartz 2.7.5
   Available from the [XQuartz project](http://xquartz.macosforge.org/landing).

3. Install Homebrew
   Follow the installation and post-installation [instructions](http://brew.sh).

4. Open a terminal and install gfortran (required for SciPy) using Homebrew
   (this can take an hour to install):  
   `> brew install gfortran`  
   If the above fails run this instead:  
   `> brew install --env=std gfortran`
   
5. Install Python 2.7.5+ and set it as the default python executable:  
   `> brew install python`  
   `> echo 'export PATH=/usr/local/bin:$PATH' >> ~/.profile`  
   Exit and reopen the terminal application so the PATH setting takes effect

6. Install NumPy, SciPy, pandas and Biopyton using the Python package manager:  
   `> pip install numpy`  
   `> pip install scipy`  
   `> pip install pandas`  
   `> pip install biopython`  
   
7. Unzip the pRESTO bundle, open a terminal window, and run
   `python setup.py install` from the pRESTO folder.
