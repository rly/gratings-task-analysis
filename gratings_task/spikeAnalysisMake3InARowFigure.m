clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';

filePathEnds = {...
%         'M20170211/M20170211_PUL_15d-g2-g3-g4-g5-g6-evokedSpiking-v4.mat' % evoked, no delay PM
        'M20170211/M20170211_PUL_6d-g2-g3-g4-g5-g6-evokedSpiking-v5.mat' % no delay 1, delay2 PLd
        'M20170608/M20170608_PUL_43a-g2-g3-g4-evokedSpiking-v5.mat' % delay1, no delay2 PLv
%         'M20170211/M20170211_PUL_8b-g2-g3-g4-g5-g6-evokedSpiking-v2.mat' % delay1 and weak delay2, no evoked (also response) PM
        }; 


processedDataDir = sprintf('%s/posterFigs', processedDataRootDir);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end

nLoc = 4;

for i = 1:size(filePathEnds)
    evokedSpikingFilePath = sprintf('%s/%s', processedDataRootDir, filePathEnds{i});
    indexSlash = strfind(filePathEnds{i}, '/');
    indexES = strfind(filePathEnds{i}, '-evokedSpiking');
    filePathEndRoot = filePathEnds{i}(indexSlash(1)+1:indexES-1);
    
%     plotFileName = sprintf('%s/%s-posterFig3.png', processedDataDir, filePathEndRoot);
%     fprintf('\tSaving figure to file %s...\n', plotFileName);
%     quickSpdfAllEvents3InARow(evokedSpikingFilePath, nLoc, plotFileName);
    
    plotFileName = sprintf('%s/%s-posterFig5.png', processedDataDir, filePathEndRoot);
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    quickSpdfAllEvents5InARow(evokedSpikingFilePath, nLoc, plotFileName);
    
    print(sprintf('%s/%s-posterFig5.svg', processedDataDir, filePathEndRoot), '-dsvg');
end