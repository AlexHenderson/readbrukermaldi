function [mass,spectra,foldernames,filenames] = readbrukermaldi(foldernames)
% READBRUKERMALDI Reads the Bruker Flex spectrum file format
% Version 1.3
%
% usage: 
% [mass,spectra,filenames] = readbrukermaldi(foldernames);
% or
% [mass,spectra,filenames] = readbrukermaldi();
%  (The second version prompts for one or more folder names.)
%
% Takes zero, one or more file names. 
% Returns:  'mass' a mass vector
%           'spectra' a matrix of spectra in rows
%           'foldernames' a matrix of folder names used in the order the 
%                         spectra appear
% There is NO binning into fixed mass steps. Where the spectra are
% misaligned the data is interpolated (linearly). The maximum and minimum
% mass limits are determined by the spectrum with the smallest mass range,
% such that the spectra matrix only contains the mass range that overlaps
% all input spectra. The data are aligned such that each column of the
% spectra matrix corresponds to the same mass.
%
%   Copyright (c) 2015-2017, Alex Henderson.
%   Contact email: alex.henderson@manchester.ac.uk
%   Licenced under the GNU Lesser General Public License (LGPL) version 3
%   https://www.gnu.org/copyleft/lesser.html
%   Other licensing options are available, please contact Alex for details
%   If you use this file in your work, please acknowledge the author(s) in
%   your publications. 
%
% Version 1.3, Alex Henderson April 2017

% Version 1.3, Alex Henderson April 2017
%   Incorporated find_value2 function. 
% Version 1.2, Alex Henderson April 2015
%   Fixed bug when passing in a folder name. 
% Version 1.1, Alex Henderson March 2015
%   Added Linux compatibility. 
% Version 1.0  Alex Henderson, February 2015
%   Based on biotofspectrum version 1.1

% Code to read the file format is at the end of this file

% The Bruker file format comprises a folder with a collection of files in
% subfolders. We're interested in two files: fid and acqus. The problem
% with reading multiple Bruker spectra is that they can't appear in the
% same folder since they'd overwrite one another. Therefore we need to have
% a method of selecting multiple folders rather than files. MATLAB doesn't
% have this built in, so we're using a component from the MATLAB File
% Exchange: uipickfiles
% http://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles--uigetfile-on-steroids

%% Check whether uipickfiles function is available
if ~exist('uipickfiles','file')
    error('This function relies ''uipickfiles'' available from MATLAB File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles--uigetfile-on-steroids');
end

%% Check whether the user has supplied a list of folders
if (~exist('foldernames', 'var'))
    foldernames = uipickfiles();
    if (isfloat(foldernames) && (foldernames==0))
        % Nothing chosen
        return;
    end
end

% The uipickfiles code returns a cell array of filenames of dimension 
% 1 x many. We need a list of dimension many x 1. However, we may have been
% passed a column of cells containing foldernames. Another possibility is
% that we have a char array passed in as a row vector.
if (iscell(foldernames))
     if (size(foldernames,2) > 1)
        foldernames = foldernames';
     end
end
numberoffolders = size(foldernames,1);

% With the foldernames (either supplied or chosen using the GUI) we look
% for files called 'fid'

% We need to traverse a number of folders since the Windows DIR command
% can't look for a specific file name in a remote folder. Therefore we need
% to record where we start from so we can return back there later. 
wherewewere = pwd;

% Somewhere to put the filenames
filenames = [];

try
    % We're using a try/catch mechanism to make sure we return to the
    % original folder if we get an error. 

    for i=1:numberoffolders
        folder = char(foldernames(i,:));
        % Move into the folder
        cd(folder);
        % Look for files called 'fid' in all subfolders
        status = 1;
        if (ispc())
            [status,result] = system('dir fid /S /B');
        else
            if (isunix())
                [status,result] = system(['find "', folder, '" -name "fid"']);
            end
        end
            
        if (status == 0)
            % We found something
            try 
                % strsplit introduced after R2012a
                filesinfolder = strsplit(result,'\n');
            catch exception
                [filesinfolder, matches] = regexp(result, '\n', 'split', 'match'); 
            end
            filenames = vertcat(filenames,filesinfolder(1:end-1)');
        end
    end

catch exception
    
    % Go back to the original folder if there is an error
    cd(wherewewere);
    error(exception.identifier, exception.message);
end

% Go back to the original folder anyway
cd(wherewewere);

numberoffiles = size(filenames,1);

%% Now we have a list of 'fid' files so we can open them and extract the spectra

for i=1:numberoffiles

    [mass_i, spectrum_i] = brukerflex(char(filenames(i,:)));
    if (i==1)
        % First time through we initialise the data array
        spectra=zeros(numberoffiles,length(spectrum_i));
        spectra(1,:)=spectrum_i;
        mass=mass_i;
    end
    
    needtointerpolate=0;
    if(length(mass_i) ~= length(mass))
        % Different number of data points, so we need to interpolate the
        % data. This is examined separately otherwise, if the two vectors
        % are of different length, MATLAB raises an error. 
        needtointerpolate=1;
    else
        if(mass_i ~= mass)
            % Different values of mass, so we need to interpolate the data
            needtointerpolate=1;
        end
    end
    
    if(needtointerpolate)
        % First determine the range over which the mismatched spectra
        % overlap. Then truncate both the mass vector and data matrix for
        % the data already processed and the new data.
        
       lowmass=max(mass(1), mass_i(1));
       highmass=min(mass(end), mass_i(end));
       
       idx=find_value2(mass,[lowmass,highmass]);
       mass=mass(idx(1):idx(2));
       spectra=spectra(:,idx(1):idx(2));
       
       idx=find_value2(mass_i,[lowmass,highmass]);
       mass_i=mass_i(idx(1):idx(2));
       spectrum_i=spectrum_i(idx(1):idx(2));
       
       % Now interpolate the new spectrum vector to match the existing
       % data.       
       spectrum_i=interp1(mass_i,spectrum_i,mass,'linear');
    end
    
    spectra(i,:)=spectrum_i;
end

% Sometimes the interpolation turns up a NaN in either the first or last
% channel. Possibly both. Here we truncate the data to remove them. 
if(find(isnan(spectra(:,1))))
    mass=mass(2:end);
    spectra=spectra(:,2:end);
end
if(find(isnan(spectra(:,end))))
    mass=mass(1:end-1);
    spectra=spectra(:,1:end-1);
end

end % function readbrukermaldi

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [mass, spectrum] = brukerflex(fidfilename)

% BRUKERFLEX Reads the Bruker Flex spectrum file format
% Version 1.0
%
% syntax: [mass, spectrum] = brukerflex(filename);
% or
% syntax: [mass, spectrum] = brukerflex();
%
% The second version prompts for a file name.
%
% Version 1.0
% Reads the Bruker Flex spectrum file format
%
% Takes a file name
% Returns two equal length column vectors - mass and spectrum
% There is NO binning into fixed mass steps 
%
%   Copyright (c) Alex Henderson, February 2015
%   Contact email: alex.henderson@manchester.ac.uk
%   Licenced under the GNU Lesser General Public License (LGPL) version 3
%   https://www.gnu.org/copyleft/lesser.html
%   Other licensing options are available, please contact Alex for details
%   If you use this file in your work, please acknowledge the author(s) in
%   your publications. 
%
% Version 1.0, February 2015

% Version 1.0  Alex Henderson, February 2015
%   Based on biotofspectrum version 1.1

% Credits to...
% Sebastian Gibb for interpretation of the calibration contants. 
% https://github.com/sgibb/readBrukerFlexData
% Implementation of the mass calculation taken from:
% https://github.com/sgibb/readBrukerFlexData/commits/master/R/tof2mass-functions.R
% Commit number 75a05247e631b2b5c9ba9698f6acff83df57e79e (initial import)
%
% Martin Strohalm and BH Clowers for interpretation of the tof spacing and
% offsets. 
% https://code.google.com/p/toolz/
% Implementation taken from: 
% http://toolz.googlecode.com/svn-history/r43/trunk/SimpleView/flexReader.py
%
% If you need R or Python code to read the Bruker Flex format, check out
% these resources.

try
    [filehandle,message] = fopen(fidfilename);
    if(filehandle == -1)
        error(message); 
    end;

    spectrum = fread(filehandle,inf,'uint');
    fclose(filehandle);
catch exception
    fclose(filehandle);
    error(exception.identifier, exception.message);
end
    

[pathstr] = fileparts(fidfilename);
acqusfilename = fullfile(pathstr,'acqus');

acqus = fileread(acqusfilename);

matchStr = regexp(acqus,'##\$ML1= (.+?) ','tokens');
constant1 = str2double(matchStr{1});

matchStr = regexp(acqus,'##\$ML2= (.+?) ','tokens');
constant2 = str2double(matchStr{1});

matchStr = regexp(acqus,'##\$ML3= (.+?) ','tokens');
constant3 = str2double(matchStr{1});

matchStr = regexp(acqus,'##\$DELAY= (.+?) ','tokens');
delay = str2double(matchStr{1});

matchStr = regexp(acqus,'##\$DW= (.+?) ','tokens');
step = str2double(matchStr{1});

numpoints = length(spectrum);

tof = transpose(delay : step : (delay+step*(numpoints-1)));

picoseconds = 1e12;

a = constant3;
b = sqrt(picoseconds / constant1);
c = constant2 - tof;

mass = (-b + sqrt((b * b) - (4 * a * c))) / (2 * a);
sign = 1;
if (mass  < 0)
    sign = -1;
end
mass = (mass .^2) * sign;

end % function brukerflex

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [vec_index]=find_value2(vector,value)

% This function finds the indices for a value in vector
%
% Normal vector, not saisir structure
%
% Usage [vec_index]=find_value2(vector,value)

% Modified version of find_value to allow more than one value to be
% calculated. Alex Henderson, January 2013


if (numel(value) > 1)
    vec_index = zeros(size(value));
    for i = 1:length(value)
        [output_value,vec_index(i)] = min(abs(vector-value(i)));
    end
else
    [output_value,vec_index] = min(abs(vector-value));
end

end
