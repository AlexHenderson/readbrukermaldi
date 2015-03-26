# readbrukermaldi

**MATLAB code to read the Bruker Flex file format used to encode MALDI files.**

## Version

	Version 1.1
	March 2015

## Usage

	[mass,spectra,foldernames,filenames] = readbrukermaldi(foldernames);

**or**

	[mass,spectra,foldernames,filenames] = readbrukermaldi();

(The second version prompts for one or more folder names.)

## Input

	Zero, one or more folder names.

## Output

	'mass' a mass vector
    'spectra' a matrix of spectra in rows
    'foldernames' a matrix of folder names used in the order the spectra appear
    'filenames' a matrix of fid file names used in the order the spectra appear

## Requirements

The Bruker file format comprises a folder with a collection of files in subfolders. We're interested in two files: fid and acqus. The problem with reading multiple Bruker spectra is that they can't appear in the same folder since they'd overwrite one another. Therefore we need to have
a method of selecting multiple folders rather than files. MATLAB doesn't have this built in, so we're using a component from the MATLAB File Exchange: **uipickfiles** [http://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles--uigetfile-on-steroids
](http://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles--uigetfile-on-steroids)
## Notes

There is NO binning into fixed mass steps. Where the spectra are misaligned, the data is interpolated (linearly). The maximum and minimum mass limits are determined by the spectrum with the smallest mass range, such that the spectra matrix only contains the mass range that overlaps all input spectra. The data are aligned such that each column of the spectra matrix corresponds to the same mass.

Some of the ideas were taken from the following projects:

- Sebastian Gibb's **readBrukerFlexData**, implemented in R, for interpretation of the calibration constants.  
[https://github.com/sgibb/readBrukerFlexData](https://github.com/sgibb/readBrukerFlexData), and,
- Implementation of the mass calculation taken from:
[https://github.com/sgibb/readBrukerFlexData/commits/master/R/tof2mass-functions.R](https://github.com/sgibb/readBrukerFlexData/commits/master/R/tof2mass-functions.R) 
Commit number 75a05247e631b2b5c9ba9698f6acff83df57e79e (initial import)
- Martin Strohalm and BH Clowers' **toolz** project, implemented in Python, for interpretation of the tof spacing and offsets. 
[https://code.google.com/p/toolz/](https://code.google.com/p/toolz/)
Information taken from: 
[http://toolz.googlecode.com/svn-history/r43/trunk/SimpleView/flexReader.py](http://toolz.googlecode.com/svn-history/r43/trunk/SimpleView/flexReader.py)

If you need R or Python code to read the Bruker Flex format, I recommend you check out these resources.

##Updates
A more updated version of this code may be found at: [https://github.com/AlexHenderson/readbrukermaldi](https://github.com/AlexHenderson/readbrukermaldi "MATLAB FileExchange")

## Legal bit

Copyright (c) Alex Henderson, 2015

Licensed under the GNU Lesser General Public License (LGPL) version 3 [https://www.gnu.org/copyleft/lesser.html](https://www.gnu.org/copyleft/lesser.html) 

Other licensing options are available, please contact Alex for details. If you use this file in your work, please acknowledge the author(s) in your publications. 

## Contact

alex.henderson @ manchester.ac.uk

