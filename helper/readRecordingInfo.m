function recordingInfo = readRecordingInfo()

%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\recordingInfo.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/10/21 15:14:57

%% Initialize variables.
filename = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\recordingInfo.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: text (%q)
%	column2: text (%q)
%   column3: text (%q)
%	column4: text (%q)
%   column5: text (%q)
%	column6: text (%q)
%   column7: text (%q)
%	column8: text (%q)
%   column9: text (%q)
%	column10: text (%q)
%   column11: text (%q)
%   column12: text (%q)
%   column13: text (%q)
%   column14: text (%q)
%   column15: text (%q)
%   column16: text (%q)
%   column17: text (%q)
%   column18: text (%q)
%   column19: text (%q)
%   column20: text (%q)
%   column21: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

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
pl2FileName = dataArray{:, 1};
sessionName = dataArray{:, 2};
areaName = dataArray{:, 3};
blockNames = cellfun(@(x) strsplit(x, ', '), dataArray{:, 4}, 'UniformOutput', false);
gratingsTask3DIndices = cellfun(@(x) cellfun(@(y) str2double(y), strsplit(x, ', ')), dataArray{:, 5}, 'UniformOutput', false);
gratingsTask3DLogIndices = cellfun(@(x) cellfun(@(y) str2double(y), strsplit(x, ', ')), dataArray{:, 6}, 'UniformOutput', false);
gratingsTask0DIndices = cellfun(@(x) cellfun(@(y) str2double(y), strsplit(x, ', ')), dataArray{:, 7}, 'UniformOutput', false);
gratingsTask0DLogIndices = cellfun(@(x) cellfun(@(y) str2double(y), strsplit(x, ', ')), dataArray{:, 8}, 'UniformOutput', false);
vepmIndices = cellfun(@(x) cellfun(@(y) str2double(y), strsplit(x, ', ')), dataArray{:, 9}, 'UniformOutput', false);
aepmIndices = cellfun(@(x) cellfun(@(y) str2double(y), strsplit(x, ', ')), dataArray{:, 10}, 'UniformOutput', false);
spikeChannelPrefix = dataArray{:, 11};
spikeChannelsToLoad = cellfun(@(x) eval(x), dataArray{:, 12}, 'UniformOutput', false); % warning: matlab injection security risk
muaChannelsToLoad = cellfun(@(x) eval(x), dataArray{:, 13}, 'UniformOutput', false); % warning: matlab injection security risk
lfpChannelsToLoad = cellfun(@(x) eval(x), dataArray{:, 14}, 'UniformOutput', false); % warning: matlab injection security risk
spkcChannelsToLoad = cellfun(@(x) eval(x), dataArray{:, 15}, 'UniformOutput', false); % warning: matlab injection security risk
directChannelsToLoad = cellfun(@(x) eval(x), dataArray{:, 16}, 'UniformOutput', false); % warning: matlab injection security risk
pldChannels = cellfun(@(x) eval(x), dataArray{:, 17}, 'UniformOutput', false); % warning: matlab injection security risk
plvChannels = cellfun(@(x) eval(x), dataArray{:, 18}, 'UniformOutput', false); % warning: matlab injection security risk
pmChannels = cellfun(@(x) eval(x), dataArray{:, 19}, 'UniformOutput', false); % warning: matlab injection security risk
piChannels = cellfun(@(x) eval(x), dataArray{:, 20}, 'UniformOutput', false); % warning: matlab injection security risk
notes = dataArray{:, 21};


%% Place vars into struct
info = [pl2FileName, sessionName, areaName, blockNames, gratingsTask3DIndices, ...
        gratingsTask3DLogIndices, vepmIndices, aepmIndices, spikeChannelPrefix, ...
        spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad, ...
        pldChannels, plvChannels, pmChannels, piChannels];
headers = {'pl2FileName', 'sessionName', 'areaName', 'blockNames', 'gratingsTask3DIndices', ...
        'gratingsTask3DLogIndices', 'vepmIndices', 'aepmIndices', 'spikeChannelPrefix', ...
        'spikeChannelsToLoad', 'muaChannelsToLoad', 'lfpChannelsToLoad', 'spkcChannelsToLoad', 'directChannelsToLoad', ...
        'pldChannels', 'plvChannels', 'pmChannels', 'piChannels'};
recordingInfo = cell2struct(info, headers, 2);

