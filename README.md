Version 0.4.6:  May 13, 2015
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
Python 2.7.5    | http://python.org
NumPy 1.7       | http://numpy.org
SciPy 0.12      | http://scipy.org
pandas 0.12     | http://pandas.pydata.org
Biopython 1.59  | http://biopython.org
MUSCLE v3.8     | http://www.drive5.com/muscle
USEARCH v7.0    | http://www.drive5.com/usearch


Installation - Linux
-------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the full 
   SciPy stack using the instructions at http://scipy.org/install.html, then 
   install Biopython according to the 
   [instructions](http://biopython.org/DIST/docs/install/Installation.html).

2. Unzip the pRESTO bundle into a directory of your choice and add that 
   directory to your `$PATH`.


Installation - Windows
-------------------------------------------------------------------------------

1. Install Python 2.7.5+ from [Python](http://python.org/download).

2. Install NumPy, SciPy, pandas and Biopython using the packages available 
   from the [Unofficial Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs)
   binary collection.

3. Unzip the pRESTO bundle into a directory of your choice and add that 
   directory to your `%Path%`.  On Windows 7 the `%Path%` setting is located under
   'Control Panel' -> 'System and Security' -> 'System' -> 
   'Advanced System Settings' -> 'Environment variables' -> 'System variables' 
   -> 'Path'

4. The pRESTO scripts should then be directly executable from the Command Prompt.
   If not, correct the file association for `.py` files by right-clicking on a 
   `.py` file, selecting 'Open with' -> 'Choose default program...', choosing the 
   `python.exe` from Python 2.7, and checking 'Always use the selected program'.
   

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
   
7. Add the pRESTO installation to your `PATH` setting. For example,
   if you copy the pRESTO scripts into `/Users/Username/presto`, then set:  
   `> echo 'export PATH=$HOME/presto:$PATH' >> ~/.profile`  
   Exit and reopen the terminal application so the PATH setting takes effect
