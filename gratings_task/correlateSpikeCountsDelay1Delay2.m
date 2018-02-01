clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';

filePathEnds = {...
        'M20170211/M20170211_PUL_20a-g2-g3-g4-g5-g6-evokedSpiking.mat' % evoked, no delay PM
        'M20170211/M20170211_PUL_25a-g2-g3-g4-g5-g6-evokedSpiking.mat' % delay1, no delay2 PI
        'M20170211/M20170211_PUL_6d-g2-g3-g4-g5-g6-evokedSpiking.mat' % delay2, no delay1 PLd
        'M20170211/M20170211_PUL_8b-g2-g3-g4-g5-g6-evokedSpiking.mat' % delay1 and delay2, no evoked (also response) PM
        }; 


processedDataDir = sprintf('%s/posterFigs', processedDataRootDir);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

nLoc = 4;

for i = 1:size(filePathEnds)
    evokedSpikingFilePath = sprintf('%s/%s', processedDataRootDir, filePathEnds{i});
    indexSlash = strfind(filePathEnds{i}, '/');
    indexES = strfind(filePathEnds{i}, '-evokedSpiking.mat');
    filePathEndRoot = filePathEnds{i}(indexSlash(1)+1:indexES-1);
    
    plotFileName = sprintf('%s/%s-posterFig.png', processedDataDir, filePathEndRoot);
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    
    ES = load(evokedSpikingFilePath);
    cueTargetDelayInRFRate = ES.averageFiringRatesByCount.cueTargetDelayHold.trialRateByLoc{ES.inRFLoc};
    targetDimDelayInRFRate = ES.averageFiringRatesByCount.targetDimDelay.trialRateByLoc{ES.inRFLoc};
    
    figure;
    hold on;
    plot(cueTargetDelayInRFRate, targetDimDelayInRFRate, ...
            '.', 'MarkerSize', 10);
    
    [sortCueTargetDelayInRFRate,sortCTDelayInd] = sort(cueTargetDelayInRFRate);
    targetDimDelayInRFRateSortedByCTDelay = targetDimDelayInRFRate(sortCTDelayInd);
    
    nTrials = numel(cueTargetDelayInRFRate);
    topThirdIndices = round(nTrials * 2/3)+1:nTrials;
    bottomThirdIndices = 1:round(nTrials * 1/3);
    
    plot(sortCueTargetDelayInRFRate(topThirdIndices), ...
            targetDimDelayInRFRateSortedByCTDelay(topThirdIndices), ...
            '.', 'MarkerSize', 10);
        
    [h,p] = ttest2(targetDimDelayInRFRateSortedByCTDelay(topThirdIndices), ...
            targetDimDelayInRFRateSortedByCTDelay(bottomThirdIndices))
    % ideally shuffle test?
    % how to deal with correlated increases in firing in a trial?
    
    
end