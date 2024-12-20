Installation
================================================================================

The simplest way to install the latest stable release of pRESTO is via pip::

    > pip3 install presto --user

The current development build can be installed using pip and git in a
similar fashion::

    > pip3 install git+https://github.com/immcantation/presto@master --user

If you currently have a development version installed, then you will likely
need to add the arguments ``--upgrade --no-deps --force-reinstall`` to the
pip3 command.

Requirements
--------------------------------------------------------------------------------

+  `Python 3.4.0 <https://python.org>`__
+  `setuptools 2.0 <https://github.com/pypa/setuptools>`__
+  `NumPy 1.8 <https://numpy.org/>`__
+  `SciPy 0.14 <https://scipy.org/>`__
+  `pandas 0.24 <https://pandas.pydata.org>`__
+  `Biopython 1.77 <https://biopython.org>`__
+  AlignSets requires `MUSCLE v3.8 <https://www.drive5.com/muscle>`__
+  ClusterSets requires `USEARCH v7.0 <https://www.drive5.com/usearch>`__,
   `vsearch v2.3.2 <https://github.com/torognes/vsearch>`__, or
   `CD-HIT v4.6.8 <https://weizhongli-lab.org/cd-hit/>`__
+  AssemblePairs-reference requires `USEARCH v7.0 <https://www.drive5.com/usearch>`__
   or `BLAST+ 2.5 <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>`__

Linux
--------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the
   full SciPy stack using the
   `instructions <https://scipy.org/install.html>`__, then install
   Biopython according to its
   `instructions <https://biopython.org/DIST/docs/install/Installation.html>`__.

2. Download the `pRESTO bundle <https://github.com/immcantation/presto/tags>`__ and run::

   > pip3 install presto-x.y.z.tar.gz --user

Mac OS X
--------------------------------------------------------------------------------

1. Install Xcode, which is available from the Apple store or
   `developer downloads <https://developer.apple.com/downloads>`__.

2. Older versions Mac OS X will require you to install XQuartz 2.7.5, which is available
   from the `XQuartz project <https://xquartz.macosforge.org/landing>`__.

3. Install Homebrew following the installation and post-installation
   `instructions <https://brew.sh>`__.

4. Install Python 3.4.0+ and set the path to the python3 executable::

   > brew install python3
   > echo 'export PATH=/usr/local/bin:$PATH' >> ~/.profile

5. Exit and reopen the terminal application so the PATH setting takes effect.

6. You may, or may not, need to install gfortran (required for SciPy). Try
   without first, as this can take an hour to install and is not needed on
   newer releases. If you do need gfortran to install SciPy, you can install it
   using Homebrew::

   > brew install gfortran

   If the above fails run this instead::

   > brew install --env=std gfortran

7. Install NumPy, SciPy, pandas and Biopython using the Python package
   manager::

   > pip3 install numpy scipy pandas biopython

8. Download the `pRESTO bundle <https://github.com/immcantation/presto/tags>`__,
   open a terminal window, change directories
   to download location, and run::

   > pip3 install presto-x.y.z.tar.gz

Windows
--------------------------------------------------------------------------------

1. Install Python 3.4.0+ from `Python <https://python.org/downloads>`__,
   selecting both the options 'pip' and 'Add python.exe to Path'.

2. Install NumPy, SciPy, pandas and Biopython.

3. Download the `pRESTO bundle <https://github.com/immcantation/presto/tags>`__,
   open a Command Prompt, change directories to
   the download folder, and run::

   > pip install presto-x.y.z.tar.gz

4. For a default installation of Python 3.4, the pRESTO scripts will be
   installed into ``C:\Python34\Scripts`` and should be directly
   executable from the Command Prompt. If this is not the case, then
   follow step 5 below.

5. Add both the ``C:\Python34`` and ``C:\Python34\Scripts`` directories
   to your ``%Path%``. On Windows 7 the ``%Path%`` setting is located
   under Control Panel -> System and Security -> System -> Advanced
   System Settings -> Environment variables -> System variables -> Path.

6. If you have trouble with the ``.py`` file associations, try adding ``.PY``
   to your ``PATHEXT`` environment variable. Also, opening a
   command prompt as Administrator and run::

    > assoc .py=Python.File
    > ftype Python.File="C:\Python34\python.exe" "%1" %*
