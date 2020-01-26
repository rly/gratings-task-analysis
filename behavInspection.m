clear
close all
outputDir = '/Users/labmanager/Documents/MATLAB/BurstSep4/';
v = 14;

fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigMcCartney.csv.bak');
% fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigFerdy.csv.bak');
sessionInfo = textscan(fid, '%d8%s%s', 'Delimiter', ',', 'HeaderLines' ,1-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fid);
counterFR = zeros(1,numel(sessionInfo{1}));
nUnitsPerSession = zeros(1,numel(sessionInfo{1}));

for sessioni = 1:numel(sessionInfo{1})
    clearvars -except sessionInfo sessioni outputDir counterFR nUnitsPerSession v
    sessionInd = sessionInfo{1}(sessioni);
    processedDataRootDir = '/Volumes/scratch/rly/gratings-task-analysis/processed_data/';
    dataDirRoot = '/Volumes/kastner/ryanly/McCartney/merged';
    %     dataDirRoot = '/Volumes/kastner/ryanly/Ferdy/merged';
    muaDataDirRoot = '/Volumes/scratch/rly/simple-mua-detection/processed_data/';
    suaMuaDataDirRoot = muaDataDirRoot;
    recordingInfoFileName = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
    channels = strsplit(sessionInfo{3}{sessioni},'-');
    channelsToLoad = 0;%str2double(channels{1}):str2double(channels{2});
    muaChannelsToLoad = 0;%channelsToLoad;
    lfpChannelsToLoad = 0;%str2double(channels{1}):str2double(channels{2});
    lfpChannels = 0;%lfpChannelsToLoad;
    isZeroDistractors = 0;
    numRandomizations = 0;
    isLoadSortedSua = 0;
    isLoadMua = 0;
    isLoadAllSpikes = 0;
    %%
    if isZeroDistractors
        scriptName = 'SUA_MUA_GRATINGS_0D';
        taskName = 'GRATINGS_0D';
    else
        scriptName = 'SUA_MUA_GRATINGS';
        taskName = 'GRATINGS';
    end
    %% Plot raster, firing rate, pre-post ISI, waveform/hist & BRI
        [R, D, processedDataDir, blockName] = loadRecordingData(processedDataRootDir, ...
            dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, ...
            taskName, scriptName, isLoadSortedSua, isLoadMua, 0, 0, ...
            [],[],isLoadAllSpikes);
    
    sessionName = R.sessionName;
    fprintf('Processing %s...\n', sessionName);
    dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
    gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));
    UE = getUsefulEvents2(gratingsTaskLogDir, R.gratingsTaskLogIndices, 4, D, blockName);
    dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
    gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));
    % totalTimeOverall = sum(D.blockStopTimes(R.blockIndices) - D.blockStartTimes(R.blockIndices));
    minFiringRateOverall = 5; % Hz
    nUnits = numel(D.allUnitStructs);
    % kCluster = 4;
    % centroid = zeros(nUnits,kCluster,2);
    
    FigH = figure; histogram2(UE.targetDimDelayDur(UE.targetDimDelayDur>0),UE.rt(UE.targetDimDelayDur>0),'BinWidth',[25 .025],'DisplayStyle','tile')
    xlim([450 1100]);ylim([0.250 0.750])
    plotFileName = sprintf('%s/Behavior-sessionInd%d-sessioni%d-v%d.png', outputDir,sessionInd, sessioni, v);
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    saveas(FigH, plotFileName);
    close all
end




