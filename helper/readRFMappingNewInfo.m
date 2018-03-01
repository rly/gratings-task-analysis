function rfMappingNewInfo = readRFMappingNewInfo(filename)

%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\rfMappingNewInfo.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2018/02/25 18:47:05

%% Initialize variables.
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%	column3: double (%f)
%   column4: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
sessionName = dataArray{:, 1};
blockInd = num2cell(dataArray{:, 2});
mode = num2cell(dataArray{:, 3});
resultsFileName = dataArray{:, 4};

%% Place vars into struct
info = [sessionName, blockInd, mode, resultsFileName];
headers = {'sessionName', 'blockInd', 'mode', 'resultsFileName'};
rfMappingNewInfo = cell2struct(info, headers, 2);
