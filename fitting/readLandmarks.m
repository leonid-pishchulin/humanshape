%{
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Author: Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
%}

function [landmarks] = readLandmarks(filename)

%lndFileLandmarkPoints loads the points from the file 'filename' into the 'LandmarkPoints' cell array

% check number and type of arguments
if nargin < 1
    error('Function requires one input argument');
elseif ~isstr(filename)
    error('Input must be a string representing a filename');
end


% Open the file.  If this returns a -1, we did not open the file
% successfully.
fid = fopen(filename);
if fid==-1
    error('File not found or permission denied');
end


% Initializing flags to 0
dataFlag = 0; % 1 if the file pointer is in the data section of the file
eofFlag = 0; % 1 if the file pointer is in the last line of the file

% Initializing the temp variables
lineNr = 0; % indicates the nr of lines read
char CurrentLine[] = ''; % variable to store the current line contents

landmarks = [];

while (eofFlag == 0)
    
    CurrentLine = fgetl(fid);
    lineNr = lineNr+1;
    
    if strcmp(CurrentLine,'END =')
        dataFlag = 0;
        eofFlag =1;
    end
    
    if (dataFlag == 1)
        ids = cell2mat(textscan(CurrentLine,'%*u %d %d %*f %*f %*f %*f %*s %*s','delimiter','/t'));
        if (ids(1) == -999 || ids(2) == -999)
            landmarks = [landmarks; {[NaN] [NaN] [NaN]}];
        else
            landmarks = [landmarks;textscan(CurrentLine,'%*u %*u %*u %*f %f %f %f %*s %*s','delimiter','/t')];
        end
        
    end
    
    if strcmp(CurrentLine,'AUX =')
        dataFlag = 1;
    end
    
end

fclose(fid);

landmarks = cell2mat(landmarks);
end

