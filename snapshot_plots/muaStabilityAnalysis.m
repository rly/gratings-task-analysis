function muaStabilityAnalysis(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)

v = 10;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('MUA Unit Stability Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('MUA Channel to Load: %d\n', muaChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

%% load recording information
[R, D] = loadRecordingData(processedDataRootDir, dataDirRoot, muaDataDirRoot, ...
        recordingInfoFileName, sessionInd, muaChannelsToLoad, 'Gratings', 'Gratings', 1, 0);
assert(numel(R.blockNames) == numel(D.blockStartTimes));

nUnits = numel(D.allMUAStructs);
fprintf('Processing %d multi-units...\n', nUnits);

%% test for significant spike rate drift over the gratings task blocks
firstGTaskBlockStartTime = D.blockStartTimes(R.gratingsTask3DIndices(1));
lastGTaskBlockStopTime = D.blockStopTimes(R.gratingsTask3DIndices(end));
windowSize = 300; % seconds
movingWindowStep = floor(windowSize/2);

isUnitStable = false(nUnits, 1);
unitNames = cell(nUnits, 1);
for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitNames{j} = spikeStruct.name;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitNames{j}, j, ...
            nUnits, round(j/nUnits*100));
    
    % get all spike times 
    spikeIndices = spikeStruct.ts >= firstGTaskBlockStartTime & ...
                spikeStruct.ts <= lastGTaskBlockStopTime;
    spikeTimes = spikeStruct.ts(spikeIndices);
    
    % adjust spike times to account for gaps between blocks
    adjustFactor = zeros(size(spikeStruct.ts));
    for i = 2:numel(R.gratingsTask3DIndices)
        spikeIndices2 = spikeStruct.ts >= D.blockStartTimes(R.gratingsTask3DIndices(i)) & ...
                spikeStruct.ts <= D.blockStopTimes(R.gratingsTask3DIndices(i));
        % difference between start(i) and stop(i-1)
        adjustFactor(spikeIndices2) = D.blockStartTimes(R.gratingsTask3DIndices(i)) - D.blockStopTimes(R.gratingsTask3DIndices(i - 1));
    end
    spikeTimesAdj = spikeTimes - adjustFactor;    
    
    % bin spikes in moving, overlapping windows
    binStartTimes = firstGTaskBlockStartTime:movingWindowStep:lastGTaskBlockStopTime-windowSize;
    numWindows = numel(binStartTimes);
    spikesPerBin = nan(numWindows, 1);
    for k = 1:numWindows
        spikesPerBin(k) = sum(spikeTimesAdj >= binStartTimes(k) & spikeTimesAdj < binStartTimes(k) + windowSize);
    end
    
    % test for time series stationarity
    % one issue: this method preferentially labels decreasing time series
    % as stationary. reversing the order of the data yields different
    % results.
    [hIsSpikeRateStationaryFwd,pIsSpikeRateStationaryFwd,stat,cValue,reg] = adftest(spikesPerBin, 'model', 'TS', 'alpha', 0.1);
    [hIsSpikeRateStationaryRev,pIsSpikeRateStationaryRev,stat,cValue,reg] = adftest(spikesPerBin(end:-1:1), 'model', 'TS', 'alpha', 0.1);
    
    % test for monotonic trend with strong alpha
    % REALLY SHOULD be using the reverse test where H0 = there is a trend
    % or the above stationarity test, perhaps with a bigger window size
    [hIsSpikeRateTrend,pIsSpikeRateTrend] = Mann_Kendall(spikesPerBin, 0.005);
    
    minMeanSpikeRate = 0.2;
    if mean(spikesPerBin) / windowSize >= minMeanSpikeRate
        if ~hIsSpikeRateTrend
            fprintf('Hooray! Null hypothesis not rejected...\n');
            isUnitStable(j) = 1;
        end

        if hIsSpikeRateStationaryFwd && hIsSpikeRateStationaryRev
            fprintf('Hooray! Spike rate is stationary.\n');
        end
    end
    
    % alternatively, test whether bins are drawn from a uniform
    % distribution
    
%     figure;
%     plot(spikesPerBin);

%     if hIsSpikeRateStationary
%         title('Hooray! Spike rate is stationary.');
%         fprintf('Hooray! Spike rate is stationary.\n');
%         stop
%     else
%         title('Spike rate has a unit root. It is not stationary.');
%     end
end

fprintf('\nStable units:\n');
cellfun(@(x) fprintf('%s\n', x), unitNames(isUnitStable));
