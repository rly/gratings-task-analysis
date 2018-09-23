function D = loadPL2(fileName, suaMuaDataDirRoot, sessionName, areaName, isLoadSortedSua, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad)

fprintf('----------- %s -----------\n', fileName);
fprintf('----------- %s, %s -----------\n', sessionName, areaName);

if nargin < 15
    error('Missing arguments');
end
if nargin == 15
    spikeChannelsToLoad = makeRowVector(spikeChannelsToLoad);
    muaChannelsToLoad = makeRowVector(muaChannelsToLoad);
    lfpChannelsToLoad = makeRowVector(lfpChannelsToLoad);
    spkcChannelsToLoad = makeRowVector(spkcChannelsToLoad);
    directChannelsToLoad = makeRowVector(directChannelsToLoad);
end

%% load basic info
dataInfo = PL2GetFileIndex(fileName);

D = var2struct(fileName, sessionName, areaName, isLoadSortedSua, isLoadLfp, isLoadSpkc);

D.totalRecSec = dataInfo.DurationOfRecordingSec;
D.totalTicks = dataInfo.DurationOfRecordingTicks;
D.timestampFrequency = dataInfo.TimestampFrequency;
fprintf('Recording time: %d m, %d s. Last tick: %d. Timestamp freq: %d Hz\n', ...
        floor(D.totalRecSec/60), floor(mod(D.totalRecSec, 60)), ...
        D.totalTicks, D.timestampFrequency);

fprintf('Data file has events: \n\t');
nEventCh = 0;
D.events = cell(1);
for i = 1:numel(dataInfo.EventChannels)
    if dataInfo.EventChannels{i}.NumEvents > 0
        nEventCh = nEventCh + 1;
        fprintf('%s(%d), ', dataInfo.EventChannels{i}.Name, dataInfo.EventChannels{i}.NumEvents);
        ts = PL2EventTs(fileName, dataInfo.EventChannels{i}.Name);
        assert(strcmp(dataInfo.EventChannels{i}.Name, sprintf('EVT%02d', nEventCh)));
        D.events{nEventCh} = ts.Ts;
    end
end
fprintf('\n');
fprintf('Data file has %d event channels.\n', nEventCh);

D.blockStartTimes = PL2StartStopTs(fileName, 'start');
D.blockStopTimes = PL2StartStopTs(fileName, 'stop');
assert(numel(D.blockStartTimes) == numel(D.blockStopTimes));

fprintf('Data file has %d blocks.\n', numel(D.blockStartTimes));
    
D.allUnitStructs = {};


% %% process spike times per channel into one cell of structs
% 
% if isLoadSpikes
%     % count units and get channel indices
%     D.spikeChannelIndices = zeros(1, 125); % mapping from channel number to index in dataInfo.SpikeChannels{i}
%     D.nUnitsByCh = zeros(1, 125);
%     D.nActiveSpikeChannels = 0;
%     % go in order through the spike channels and use the sorting with the
%     % largest source number for each spike channel
%     for i = 1:numel(dataInfo.SpikeChannels)
%         if dataInfo.SpikeChannels{i}.Enabled && ...
%                 strcmp(dataInfo.SpikeChannels{i}.SourceName, spikeChannelPrefix) && ...% ONLY LOAD UNITS THAT START WITH SPK_SPKC
%                 dataInfo.SpikeChannels{i}.Channel <= 125 % SKIP extra analog signal data
%             D.spikeChannelIndices(dataInfo.SpikeChannels{i}.Channel) = i;
%             D.nUnitsByCh(dataInfo.SpikeChannels{i}.Channel) = dataInfo.SpikeChannels{i}.NumberOfUnits;
%             D.nActiveSpikeChannels = D.nActiveSpikeChannels + 1;
%         end
%     end
% 
%     fprintf('Data file has %d active spike channels.\n', D.nActiveSpikeChannels);
%     
%     D.nSpikeCh = numel(spikeChannelsToLoad);
%     D.nUnits = sum(D.nUnitsByCh(spikeChannelsToLoad));
% 
%     % save each unit into D
%     fprintf('\tProcessing %d units from %d channels into D.allSpikeStructs...\n', D.nUnits, D.nSpikeCh);
%     D.allSpikeStructs = cell(D.nUnits, 1);
%     unitInd = 0;
%     spikeFs = D.timestampFrequency;
%     for i = spikeChannelsToLoad
%         % consider unsorted waveforms too - id 0
%         % PL2Ts(spikesFileName, origSpikeVarName, 0);
%         channelInd = D.spikeChannelIndices(i);
%         if channelInd == 0
%             continue;
%         end
%         spikeChannel = dataInfo.SpikeChannels{channelInd};
%         for j = 1:spikeChannel.NumberOfUnits
%             unitInd = unitInd + 1;
%             spikeStruct = struct();
%             spikeStruct.name = sprintf('%s_%s_%d%c', ...
%                     sessionName, areaName, ...
%                     spikeChannel.Channel, 'a'+j-1);
%             spikeStruct.sessionName = sessionName;
%             spikeStruct.areaName = areaName;
%             spikeStruct.channelID = spikeChannel.Channel;
%             spikeStruct.unitID = j;
%             spikeStruct.unitIDChar = sprintf('%d%c', spikeChannel.Channel, 'a'+j-1);
%             spikeStruct.unitIndInSession = unitInd;
%             spikeStruct.threshold = spikeChannel.Threshold * spikeChannel.CoeffToConvertToUnits * 1000; % mV
%             spikeStruct.thresholdTime = spikeChannel.PreThresholdSamples / spikeChannel.SamplesPerSecond;
%             spikeStruct.Fs = spikeChannel.SamplesPerSecond;
%             origSpikeVarName = spikeChannel.Name;
%             waveInfo = PL2Waves(fileName, origSpikeVarName, j);
%             spikeStruct.wf = waveInfo.Waves;
%             spikeStruct.ts = waveInfo.Ts; % timestamps in seconds
%             if isempty(spikeStruct.ts) % don't add more more information if there are no waveforms
%                 D.allSpikeStructs{unitInd} = spikeStruct;
%                 continue;
%             end
%             spikeStruct.meanWf = mean(waveInfo.Waves);
%             spikeStruct.sdWf = std(waveInfo.Waves);
%             [~,troughInd] = min(spikeStruct.meanWf); % not necessarily the trough associated with thresh crossing
%             [postTroughPeak,relPeakInd] = max(spikeStruct.meanWf(troughInd:end)); % peak must be after trough
%             spikeStruct.troughToPeakTime = (relPeakInd-1)/spikeFs;
%             
%             smoothedWf = movmean(spikeStruct.meanWf, 5); % moving average over 5 data points
%             d1 = diff(spikeStruct.meanWf) > 0; % 1 if moving up, 0 if moving down
%             d1Smooth = diff(smoothedWf) > 0; % 1 if moving up, 0 if moving down
%             d2Smooth = diff(d1Smooth); % 0->1 = 1 = trough, 1->0 = -1 = peak
%             spikeStruct.firstTroughIndex = find(d2Smooth == 1, 1, 'first') + 1;
%             
% %             figure_tr_inch(6, 6);
% %             hold on;
% %             plot(1:numel(spikeStruct.meanWf), spikeStruct.meanWf, 'LineWidth', 2);
% %             plot(troughInd, trough, 'k.', 'MarkerSize', 40);
% %             plot(troughInd + relPeakInd - 1, postTroughPeak, 'k.', 'MarkerSize', 40);
% %             x1 = spikeStruct.meanWf(1:end-1);
% %             x1(d1 == 1) = NaN;
% %             x2 = spikeStruct.meanWf(1:end-1);
% %             x2(d1 == 0) = NaN;
% %             plot(1:numel(spikeStruct.meanWf), spikeStruct.meanWf, 'b.', 'MarkerSize', 20);
% %             plot(1:numel(x1), x1, 'r.', 'MarkerSize', 20);
% %             plot(1:numel(x2), x2, 'g.', 'MarkerSize', 20);
%             
%             spikeStruct.fullWidthHalfMax = NaN;
%             % find point closest to half max before and after max
%             if postTroughPeak > 0
%                 halfMax = 1/2 * postTroughPeak;
% %                 plot(1:numel(spikeStruct.meanWf), halfMax*ones(size(spikeStruct.meanWf)), 'm');
%                 firstDownPrePeakInd = find(d1(troughInd:troughInd+relPeakInd-2) == 0 & spikeStruct.meanWf(troughInd:troughInd+relPeakInd-2) <= halfMax, 1, 'last') + troughInd - 1;
%                 firstUpPostPeakInd = find(d1(troughInd+relPeakInd:end) == 1 & spikeStruct.meanWf(troughInd+relPeakInd:end-1) <= halfMax, 1, 'first') + troughInd + relPeakInd - 1;
%                 if isempty(firstDownPrePeakInd)
%                     firstDownPrePeakInd = troughInd;
%                 end
%                 if isempty(firstUpPostPeakInd)
%                     firstUpPostPeakInd = numel(spikeStruct.meanWf);
%                 end
%                 if halfMax >= spikeStruct.meanWf(firstUpPostPeakInd)
%                     % both series are monotonic, can reverse interp1 args
%                     x = firstDownPrePeakInd:troughInd+relPeakInd-1;
%                     y = spikeStruct.meanWf(x);
%                     prePeakHalfMaxInd = interp1(y, x, halfMax, 'linear');
%                     x = troughInd+relPeakInd:firstUpPostPeakInd;
%                     y = spikeStruct.meanWf(x);
%                     postPeakHalfMaxInd = interp1(y, x, halfMax, 'linear');
%                     spikeStruct.fullWidthHalfMax = (postPeakHalfMaxInd - prePeakHalfMaxInd)/spikeFs;
% 
% %                     plot(firstDownPrePeakInd, spikeStruct.meanWf(firstDownPrePeakInd), 'c.', 'MarkerSize', 30);
% %                     plot(firstUpPostPeakInd, spikeStruct.meanWf(firstUpPostPeakInd), 'c.', 'MarkerSize', 30);
% %                     plot(prePeakHalfMaxInd, halfMax, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
% %                     plot(postPeakHalfMaxInd, halfMax, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
%                 end
%             end
%             % TODO compute full-width half-minimum around trough (sharpness
%             % of trough)
%             
%             peakSmoothedIndices = find(d2Smooth == -1) + 1;
%             troughSmoothedIndices = find(d2Smooth == 1) + 1;
%             % if it starts with a peak that wasn't picked up by diff
%             if smoothedWf(1) > smoothedWf(2) && ...
%                     (isempty(peakSmoothedIndices) || ...
%                     peakSmoothedIndices(1) ~= 1) 
%                 peakSmoothedIndices = [1 peakSmoothedIndices];
%             end
%             % if it ends with a peak that wasn't picked up by diff
%             % note: the end tends to be noisy and may require more
%             % smoothing
%             if smoothedWf(end-1) < smoothedWf(end) && ...
%                     (isempty(peakSmoothedIndices) || ...
%                     peakSmoothedIndices(1) ~= numel(smoothedWf))
%                 peakSmoothedIndices = [peakSmoothedIndices numel(smoothedWf)];
%             end
%             % highly unlikely to start/end with a trough. ignore that case
% %             figure;
% %             hold on;
% %             plot(1:numel(spikeStruct.meanWf), spikeStruct.meanWf, 'LineWidth', 2);
% %             plot(1:numel(spikeStruct.meanWf), smoothedWf, 'LineWidth', 2);
% %             plot(1:numel(spikeStruct.meanWf), smoothedWf - 1/4 * spikeStruct.sdWf, ':', 'LineWidth', 2);
% %             plot(1:numel(spikeStruct.meanWf), smoothedWf + 1/4 * spikeStruct.sdWf, ':', 'LineWidth', 2);
% %             plot(1:numel(spikeStruct.meanWf), zeros(size(spikeStruct.meanWf)), 'k');
% %             title(spikeStruct.name, 'Interpreter', 'none');
%             % only consider peaks > 1/4 SD from 0, troughs < 1/4 SD from 0
%             minSDFactorFromZero = 1;
%             goodPeaks = false(size(peakSmoothedIndices));
%             for k = 1:numel(peakSmoothedIndices)
%                 goodPeaks(k) = smoothedWf(peakSmoothedIndices(k)) - minSDFactorFromZero * spikeStruct.sdWf(peakSmoothedIndices(k)) > 0;
%             end
%             goodTroughs = false(size(troughSmoothedIndices));
%             for k = 1:numel(troughSmoothedIndices)
%                 goodTroughs(k) = smoothedWf(troughSmoothedIndices(k)) + minSDFactorFromZero * spikeStruct.sdWf(troughSmoothedIndices(k)) < 0;
%             end
%             d2s = zeros(size(d2Smooth));
%             d2s(troughSmoothedIndices(goodTroughs)) = 't'; % trough
%             d2s(peakSmoothedIndices(goodPeaks)) = 'p'; % peak
%             d2s(d2s == 0) = [];
%             spikeStruct.peakSmoothedIndices = peakSmoothedIndices(goodPeaks);
%             spikeStruct.troughSmoothedIndices = troughSmoothedIndices(goodTroughs);
%             spikeStruct.peakSmoothedAmps = smoothedWf(spikeStruct.peakSmoothedIndices);
%             spikeStruct.troughSmoothedAmps = smoothedWf(spikeStruct.troughSmoothedIndices);
%             spikeStruct.inflectionPattern = char(d2s);
%             spikeStruct.numInflections = numel(d2s);
%             spikeStruct = classifyCellClass(spikeStruct);
%             spikeStruct.isMUA = false;
%             
% %             text(0.01, 0.05, char(d2s), 'Units', 'normalized', 'FontSize', 16);
% %             text(0.01, 0.12, sprintf('%0.3f ms', spikeStruct.fullWidthHalfMax*1000), 'Units', 'normalized', 'FontSize', 16);
% 
%             D.allSpikeStructs{unitInd} = spikeStruct;
% %             close;
%         end
%     end
%     
%     muaToRemove = zeros(numel(D.allSpikeStructs), 1);
%     for j = 1:numel(D.allSpikeStructs)
%         if isempty(D.allSpikeStructs{j}.ts)
%             muaToRemove(j) = 1;
%         end
%     end
%     fprintf('\tRemoving %d units due to lack of spikes.\n', sum(muaToRemove));
%     D.allSpikeStructs(muaToRemove == 1) = [];
%     
%     clear spikeStruct;
%     fprintf('\tdone.\n');
% end

%% process spike times per channel into one cell of structs

if isLoadSortedSua
    % count units and get channel indices
%     D.spikeChannelIndices
%     D.nUnitsByCh
%     D.nActiveSpikeChannels
%     
%     D.nSpikeCh
%     D.nUnits

    %% read spike sorting quality metrics
    sortQualityNotesFileName = 'spike sorting notes.xlsx';
    xlsSheet = 1;
    xlRange = 'A2:H3000';

    [sortQualityNotesNums,sortQualityNotesText] = xlsread(sortQualityNotesFileName, xlsSheet, xlRange);
    assert(size(sortQualityNotesNums, 2) == 6);
    assert(size(sortQualityNotesText, 2) == 8);
    sortQualityNotesText(:,2:7) = [];
    assert(size(sortQualityNotesNums, 1) == size(sortQualityNotesText, 1));
    [sortQualityNotesNums,i] = trimNanRows(sortQualityNotesNums);
    sortQualityNotesText(i,:) = [];

    % sortQualityNotesNums(:,1) = channel number
    % sortQualityNotesNums(:,2) = unit id (0 = unsorted)
    % sortQualityNotesNums(:,3) = 0 or 1 if waveform has typical shape
    % sortQualityNotesNums(:,4) = 0-5 quality of unit separation
    % sortQualityNotesNums(:,5) = start time of unit in seconds
    % sortQualityNotesNums(:,6) = end time of unit in seconds
    
    % sortQualityNotesText(:,1) = session name
    % sortQualityNotesText(:,2) = comments

    % save each unit into D
    D.allSpikeStructs = [];
    
    unitInd = 0; % index in session
    spikeFs = D.timestampFrequency;
    D.spikeChannelsToLoad = spikeChannelsToLoad;
    for i = spikeChannelsToLoad
        % load sorted SUA data
        % indexing of offline sorter starts at 0
        suaFilePath = sprintf('%s/%s-SUA_%03d.mat', suaMuaDataDirRoot, sessionName, i-1);
        suaData = load(suaFilePath, sprintf('wfData%d', i-1));
        suaData = suaData.(sprintf('wfData%d', i-1));
        % suaData should have waveform x data where
        % suaData(:,1) = channel number
        % suaData(:,2) = unit id (0 = unsorted)
        % suaData(:,3) = timestamp in ms (time of wf minimum after threshold crossing)
        % suaData(:,4:nSamples+3) = waveform in microvolts
        assert(all(suaData(:,1) == i));
        
        muaFilePath = sprintf('%s/%s-SPKC%03d-MUA.mat', suaMuaDataDirRoot, sessionName, i);
        muaData = load(muaFilePath); % for getting threshold data
        
        nUnitsThisCh = max(suaData(:,2));
        for j = 1:nUnitsThisCh
            qualityNotesMatchInd = strcmp(sortQualityNotesText(:,1), sessionName) & ...
                    sortQualityNotesNums(:,1) == i & ...
                    sortQualityNotesNums(:,2) == j;
            if ~any(qualityNotesMatchInd)
                warning('No matches in SUA quality notes for session %s, channel %d, unit %d\n', sessionName, i, j);
            elseif sum(qualityNotesMatchInd) > 1
                warning('Too many matches in SUA quality notes for session %s, channel %d, unit %d\n', sessionName, i, j);
            end
            assert(sum(qualityNotesMatchInd) == 1);
            hasTypicalWaveformShape = (sortQualityNotesNums(qualityNotesMatchInd,3) == 1);
            separationQuality = sortQualityNotesNums(qualityNotesMatchInd,4);
            % INCLUDE ONLY UNITS WITH TYPICAL WAVEFORM SHAPE AND GOOD
            % SEPARATION QUALITY FOR NOW
            if ~(hasTypicalWaveformShape && separationQuality >= 3)
                continue;
            end
            
            unitStartTime = sortQualityNotesNums(qualityNotesMatchInd,5);
            unitEndTime = sortQualityNotesNums(qualityNotesMatchInd,6);
            if isempty(unitStartTime)
                unitStartTime = D.blockStartTimes(1);
            end
            if isempty(unitEndTime)
                unitEndTime = D.blockStopTimes(end);
            end
            sortComments = sortQualityNotesText(:,2);
            
            unitInd = unitInd + 1;
            unitMatch = suaData(:,2) == j;
            spikeStruct = struct();
            spikeStruct.name = sprintf('%s_%s_%d%c', ...
                    sessionName, areaName, i, 'a'+j-1);
            spikeStruct.sessionName = sessionName;
            spikeStruct.areaName = areaName;
            spikeStruct.channelID = i;
            spikeStruct.unitID = j;
            spikeStruct.unitIDChar = sprintf('%d%c',i, 'a'+j-1);
            spikeStruct.unitIndInSession = unitInd;
            spikeStruct.Fs = spikeFs;
            spikeStruct.hasTypicalWaveformShape = hasTypicalWaveformShape;
            spikeStruct.separationQuality = separationQuality;
            spikeStruct.unitStartTime = unitStartTime;
            spikeStruct.unitEndTime = unitEndTime;
            spikeStruct.sortComments = sortComments;
            % get threshold data from muaData mat file
            spikeStruct.threshold = nanmean(muaData.thresholdParams.thresholds); % mV
            spikeStruct.thresholdTime = muaData.thresholdParams.nPreThresholdSamples / spikeFs;
            spikeStruct.thresholdParams = muaData.thresholdParams;
            spikeStruct.wf = suaData(unitMatch,4:end) / 1000; % now in millivolts
            spikeStruct.ts = suaData(unitMatch,3); % seconds
            if isempty(spikeStruct.ts) % don't add more information if there are no waveforms
                D.allSpikeStructs{unitInd} = spikeStruct;
                continue;
            end
            spikeStruct.meanWf = mean(spikeStruct.wf);
            spikeStruct.sdWf = std(spikeStruct.wf);
            [~,troughInd] = min(spikeStruct.meanWf); % not necessarily the trough associated with thresh crossing
            [postTroughPeak,relPeakInd] = max(spikeStruct.meanWf(troughInd:end)); % peak must be after trough
            spikeStruct.troughToPeakTime = (relPeakInd-1)/spikeFs;
            
            smoothedWf = movmean(spikeStruct.meanWf, 5); % moving average over 5 data points
            d1 = diff(spikeStruct.meanWf) > 0; % 1 if moving up, 0 if moving down
            d1Smooth = diff(smoothedWf) > 0; % 1 if moving up, 0 if moving down
            d2Smooth = diff(d1Smooth); % 0->1 = 1 = trough, 1->0 = -1 = peak
            spikeStruct.firstTroughIndex = find(d2Smooth == 1, 1, 'first') + 1;
            
%             figure_tr_inch(6, 6);
%             hold on;
%             plot(1:numel(spikeStruct.meanWf), spikeStruct.meanWf, 'LineWidth', 2);
%             plot(troughInd, trough, 'k.', 'MarkerSize', 40);
%             plot(troughInd + relPeakInd - 1, postTroughPeak, 'k.', 'MarkerSize', 40);
%             x1 = spikeStruct.meanWf(1:end-1);
%             x1(d1 == 1) = NaN;
%             x2 = spikeStruct.meanWf(1:end-1);
%             x2(d1 == 0) = NaN;
%             plot(1:numel(spikeStruct.meanWf), spikeStruct.meanWf, 'b.', 'MarkerSize', 20);
%             plot(1:numel(x1), x1, 'r.', 'MarkerSize', 20);
%             plot(1:numel(x2), x2, 'g.', 'MarkerSize', 20);
            
            spikeStruct.fullWidthHalfMax = NaN;
            % find point closest to half max before and after max
            if postTroughPeak > 0
                halfMax = 1/2 * postTroughPeak;
%                 plot(1:numel(spikeStruct.meanWf), halfMax*ones(size(spikeStruct.meanWf)), 'm');
                firstDownPrePeakInd = find(d1(troughInd:troughInd+relPeakInd-2) == 0 & spikeStruct.meanWf(troughInd:troughInd+relPeakInd-2) <= halfMax, 1, 'last') + troughInd - 1;
                firstUpPostPeakInd = find(d1(troughInd+relPeakInd:end) == 1 & spikeStruct.meanWf(troughInd+relPeakInd:end-1) <= halfMax, 1, 'first') + troughInd + relPeakInd - 1;
                if isempty(firstDownPrePeakInd)
                    firstDownPrePeakInd = troughInd;
                end
                if isempty(firstUpPostPeakInd)
                    firstUpPostPeakInd = numel(spikeStruct.meanWf);
                end
                if halfMax >= spikeStruct.meanWf(firstUpPostPeakInd)
                    % both series are monotonic, can reverse interp1 args
                    x = firstDownPrePeakInd:troughInd+relPeakInd-1;
                    y = spikeStruct.meanWf(x);
                    prePeakHalfMaxInd = interp1(y, x, halfMax, 'linear');
                    x = troughInd+relPeakInd:firstUpPostPeakInd;
                    y = spikeStruct.meanWf(x);
                    postPeakHalfMaxInd = interp1(y, x, halfMax, 'linear');
                    spikeStruct.fullWidthHalfMax = (postPeakHalfMaxInd - prePeakHalfMaxInd)/spikeFs;

%                     plot(firstDownPrePeakInd, spikeStruct.meanWf(firstDownPrePeakInd), 'c.', 'MarkerSize', 30);
%                     plot(firstUpPostPeakInd, spikeStruct.meanWf(firstUpPostPeakInd), 'c.', 'MarkerSize', 30);
%                     plot(prePeakHalfMaxInd, halfMax, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
%                     plot(postPeakHalfMaxInd, halfMax, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
                end
            end
            % TODO compute full-width half-minimum around trough (sharpness
            % of trough)
            
            peakSmoothedIndices = find(d2Smooth == -1) + 1;
            troughSmoothedIndices = find(d2Smooth == 1) + 1;
            % if it starts with a peak that wasn't picked up by diff
            if smoothedWf(1) > smoothedWf(2) && ...
                    (isempty(peakSmoothedIndices) || ...
                    peakSmoothedIndices(1) ~= 1) 
                peakSmoothedIndices = [1 peakSmoothedIndices];
            end
            % if it ends with a peak that wasn't picked up by diff
            % note: the end tends to be noisy and may require more
            % smoothing
            if smoothedWf(end-1) < smoothedWf(end) && ...
                    (isempty(peakSmoothedIndices) || ...
                    peakSmoothedIndices(1) ~= numel(smoothedWf))
                peakSmoothedIndices = [peakSmoothedIndices numel(smoothedWf)];
            end
            % highly unlikely to start/end with a trough. ignore that case
%             figure;
%             hold on;
%             plot(1:numel(spikeStruct.meanWf), spikeStruct.meanWf, 'LineWidth', 2);
%             plot(1:numel(spikeStruct.meanWf), smoothedWf, 'LineWidth', 2);
%             plot(1:numel(spikeStruct.meanWf), smoothedWf - 1/4 * spikeStruct.sdWf, ':', 'LineWidth', 2);
%             plot(1:numel(spikeStruct.meanWf), smoothedWf + 1/4 * spikeStruct.sdWf, ':', 'LineWidth', 2);
%             plot(1:numel(spikeStruct.meanWf), zeros(size(spikeStruct.meanWf)), 'k');
%             title(spikeStruct.name, 'Interpreter', 'none');
            % only consider peaks > 1/4 SD from 0, troughs < 1/4 SD from 0
            minSDFactorFromZero = 1;
            goodPeaks = false(size(peakSmoothedIndices));
            for k = 1:numel(peakSmoothedIndices)
                goodPeaks(k) = smoothedWf(peakSmoothedIndices(k)) - minSDFactorFromZero * spikeStruct.sdWf(peakSmoothedIndices(k)) > 0;
            end
            goodTroughs = false(size(troughSmoothedIndices));
            for k = 1:numel(troughSmoothedIndices)
                goodTroughs(k) = smoothedWf(troughSmoothedIndices(k)) + minSDFactorFromZero * spikeStruct.sdWf(troughSmoothedIndices(k)) < 0;
            end
            d2s = zeros(size(d2Smooth));
            d2s(troughSmoothedIndices(goodTroughs)) = 't'; % trough
            d2s(peakSmoothedIndices(goodPeaks)) = 'p'; % peak
            d2s(d2s == 0) = [];
            spikeStruct.peakSmoothedIndices = peakSmoothedIndices(goodPeaks);
            spikeStruct.troughSmoothedIndices = troughSmoothedIndices(goodTroughs);
            spikeStruct.peakSmoothedAmps = smoothedWf(spikeStruct.peakSmoothedIndices);
            spikeStruct.troughSmoothedAmps = smoothedWf(spikeStruct.troughSmoothedIndices);
            spikeStruct.inflectionPattern = char(d2s);
            spikeStruct.numInflections = numel(d2s);
            spikeStruct = classifyCellClass(spikeStruct);
            spikeStruct.isMUA = false;
            
%             text(0.01, 0.05, char(d2s), 'Units', 'normalized', 'FontSize', 16);
%             text(0.01, 0.12, sprintf('%0.3f ms', spikeStruct.fullWidthHalfMax*1000), 'Units', 'normalized', 'FontSize', 16);

            D.allSpikeStructs{unitInd} = spikeStruct;
%             close;
        end
    end
    
    toRemove = zeros(numel(D.allSpikeStructs), 1);
    for j = 1:numel(D.allSpikeStructs)
        if isempty(D.allSpikeStructs{j}.ts)
            toRemove(j) = 1;
        end
    end
    fprintf('\tRemoving %d units due to lack of spikes.\n', sum(toRemove));
    D.allSpikeStructs(toRemove == 1) = [];
    
    D.allUnitStructs(1:numel(D.allSpikeStructs)) = D.allSpikeStructs;
    
    clear spikeStruct;
    fprintf('\tdone.\n');
end

%% TODO
if isLoadMua
    
    muaDirString = sprintf('%s/%s-*-MUA.mat', suaMuaDataDirRoot, sessionName);
    muaDirContents = dir(muaDirString);
    muaDirContents = muaDirContents(~cellfun('isempty', {muaDirContents.date}) & cell2mat({muaDirContents.bytes}) > 0); % rm invalid entries
    D.nActiveMUAChannels = numel(muaDirContents);
    fprintf('Found data for %d MUA channels.\n', D.nActiveMUAChannels);
        
    D.nMUACh = numel(muaChannelsToLoad);

    % save each unit into D
    fprintf('\tProcessing %d MUAs into D.allMUAStructs...\n', D.nMUACh);
    D.allMUAStructs = cell(D.nMUACh, 1);
    spikeFs = D.timestampFrequency;
    muaInd = 0;
    for i = muaChannelsToLoad
        muaFilePath = sprintf('%s/%s-SPKC%03d-MUA.mat', suaMuaDataDirRoot, sessionName, i);
        muaData = load(muaFilePath);
        
        muaInd = muaInd + 1;
        muaStruct = struct();
        muaStruct.name = sprintf('%s_%s_%d%c', ...
                sessionName, areaName, ...
                i, 'M');
        muaStruct.sessionName = sessionName;
        muaStruct.areaName = areaName;
        muaStruct.channelID = i;
        muaStruct.unitID = i;
        muaStruct.unitIDChar = sprintf('%d%c', i, 'M');
        muaStruct.unitIndInSession = muaInd;       
        muaStruct.Fs = spikeFs;
        muaStruct.threshold = nanmean(muaData.thresholdParams.thresholds); % mV
        muaStruct.thresholdTime = (muaData.thresholdParams.nPreThresholdSamples + 1) / spikeFs;
        muaStruct.thresholdParams = muaData.thresholdParams;
        muaStruct.wf = muaData.wf;
        muaStruct.ts = muaData.ts; % timestamps in seconds
        if isempty(muaStruct.ts) % don't add more more information if there are no waveforms
            D.allMUAStructs{muaInd} = muaStruct;
            continue;
        end
        muaStruct.meanWf = mean(muaStruct.wf);
        muaStruct.sdWf = std(muaStruct.wf);

        % fill in same fields as spikeStruct
        muaStruct.troughToPeakTime = [];
        muaStruct.firstTroughIndex = [];
        muaStruct.fullWidthHalfMax = [];
        muaStruct.peakSmoothedIndices = [];
        muaStruct.troughSmoothedIndices = [];
        muaStruct.peakSmoothedAmps = [];
        muaStruct.troughSmoothedAmps = [];
        muaStruct.inflectionPattern = [];
        muaStruct.numInflections = [];
        muaStruct.physClass = 'MUA';
        muaStruct.isMUA = true;

        D.allMUAStructs{muaInd} = muaStruct;
    end
    
    muaToRemove = zeros(numel(D.allMUAStructs), 1);
    for j = 1:numel(D.allMUAStructs)
        if isempty(D.allMUAStructs{j}.ts)
            muaToRemove(j) = 1;
        end
    end
    fprintf('\tRemoving %d units due to lack of spikes.\n', sum(muaToRemove));
    D.allMUAStructs(muaToRemove == 1) = [];
    
    D.allUnitStructs(numel(D.allUnitStructs)+1:numel(D.allUnitStructs)+numel(D.allMUAStructs)) = D.allMUAStructs;
    
    clear muaStruct;
    fprintf('\tdone.\n');
end

%% process lfps into one matrix
if isLoadLfp
    % get fragment sizes and channel indices
    D.lfpChannelIndices = 0; % mapping from channel number to index in dataInfo.AnalogChannels{i}
    D.nLfpTime = 0;
    D.lfpFragTs = 0;
    D.lfpFragCounts = 0;
    D.lfpFs = 0;
    for i = 1:numel(dataInfo.AnalogChannels)
        if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'FP') && ...
                dataInfo.AnalogChannels{i}.NumValues > 0 && dataInfo.AnalogChannels{i}.Channel <= 125
            D.lfpChannelIndices(dataInfo.AnalogChannels{i}.Channel) = i;
            if ~D.nLfpTime
                D.nLfpTime = dataInfo.AnalogChannels{i}.NumValues;
                adInfo = PL2Ad(fileName, dataInfo.AnalogChannels{i}.Name);
                D.lfpFragTs = adInfo.FragTs;
                D.lfpFragCounts = adInfo.FragCounts;
                D.lfpFs = adInfo.ADFreq;
                clear adInfo;
            else
                assert(D.nLfpTime == dataInfo.AnalogChannels{i}.NumValues);
            end
        end
    end

    D.nLfpChAll = sum(D.lfpChannelIndices > 0);
    
    fprintf('Data file has %d LFP channels.\n', D.nLfpChAll);

    D.nLfpCh = numel(lfpChannelsToLoad);

    % save the data into D
    fprintf('\tProcessing %d LFPs...\n', D.nLfpCh);
%     D.lfps = nan(D.nLfpCh, D.nLfpTime);
    D.adjLfps = nan(D.nLfpCh, D.nLfpTime);
    D.lfpNames = cell(D.nLfpCh, 1);
    lfpChCounter = 0;
    for i = lfpChannelsToLoad
        channelInd = D.lfpChannelIndices(i);
        lfpChCounter = lfpChCounter + 1;
        adInfo = PL2Ad(fileName, dataInfo.AnalogChannels{channelInd}.Name);
%             D.lfps(lfpChCounter,:) = adInfo.Values;
        adjLfp = padNaNsToAccountForDropsPL2(adInfo);
        if lfpChCounter == 1
            % adjust size of matrix after padding
            D.adjLfps = nan(D.nLfpCh, numel(adjLfp));
        end
        D.adjLfps(lfpChCounter,:) = adjLfp;

        D.lfpNames{lfpChCounter} = sprintf('%s_%s_%s', ...
                    sessionName, areaName, ...
                    dataInfo.AnalogChannels{channelInd}.Name);
    end
    clear adInfo;
    fprintf('\tdone.\n');
end

%% process SPKC into one matrix
startTime = 1; % from 1 second in
endTime = 11; % to 11 seconds in
maxTimeToSave = (endTime - startTime)*40000;

if isLoadSpkc
    % get fragment sizes and channel indices
    D.spkcChannelIndices = 0; % mapping from channel number to index in dataInfo.AnalogChannels{i}
    D.nSpkcTime = 0;
    D.spkcFragTs = 0;
    D.spkcFragCounts = 0;
    D.spkcFs = 0;
    for i = 1:numel(dataInfo.AnalogChannels)
        if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'SPKC') && ...
                dataInfo.AnalogChannels{i}.NumValues > 0 && ...
                dataInfo.AnalogChannels{i}.Channel <= 125 % SKIP extra analog signal data
            D.spkcChannelIndices(dataInfo.AnalogChannels{i}.Channel) = i;
            if ~D.nSpkcTime
                D.nSpkcTime = dataInfo.AnalogChannels{i}.NumValues;
                adInfo = PL2Ad(fileName, dataInfo.AnalogChannels{i}.Name);
                D.spkcFragTs = adInfo.FragTs;
                D.spkcFragCounts = adInfo.FragCounts;
                D.spkcFs = adInfo.ADFreq;
                clear adInfo;
            else
                assert(D.nSpkcTime == dataInfo.AnalogChannels{i}.NumValues);
            end
        end
    end
    
    D.nSpkcChAll = sum(D.spkcChannelIndices > 0);
    fprintf('Data file has %d SPKC channels.\n', D.nSpkcChAll);

    D.nSpkcCh = numel(spkcChannelsToLoad);

    % save the data into D
    fprintf('\tProcessing %d SPKCs into D.spkcs...\n', D.nSpkcCh);
    D.spkcs = nan(D.nSpkcCh, maxTimeToSave);
    D.spkcNames = cell(D.nSpkcCh, 1);
    spkcChCounter = 0;
    for i = spkcChannelsToLoad
        channelInd = D.spkcChannelIndices(i);
        spkcChCounter = spkcChCounter + 1;
%             adInfo = PL2Ad(fileName, dataInfo.AnalogChannels{channelInd}.Name);
        adInfo = PL2AdTimeSpan(fileName, dataInfo.AnalogChannels{channelInd}.Name, startTime, endTime);
        % index time (1:maxTimeToSave) because there might be an extra timestamp
        D.spkcs(spkcChCounter,:) = adInfo.Values(1:maxTimeToSave);
        % TODO pad with NaNs around each fragment
        D.spkcNames{spkcChCounter} = sprintf('%s_%s_%s', ...
                    sessionName, areaName, ...
                    dataInfo.AnalogChannels{channelInd}.Name);
    end
    clear adInfo;
    fprintf('\tdone.\n');
end

%% process direct analog channels into one matrix
if isLoadDirect
    % get fragment sizes and channel indices
    D.directChannelIndices = 0; % mapping from channel number to index in dataInfo.AnalogChannels{i}
    D.nDirectTime = 0;
    D.directFragTs = 0;
    D.directFragCounts = 0;
    D.directFs = 0;
    for i = 1:numel(dataInfo.AnalogChannels)
        if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'Direct')
            D.directChannelIndices(dataInfo.AnalogChannels{i}.Channel) = i;
            if ~D.nDirectTime
                D.nDirectTime = dataInfo.AnalogChannels{i}.NumValues;
                adInfo = PL2Ad(fileName, dataInfo.AnalogChannels{i}.Name);
                D.directFragTs = adInfo.FragTs;
                D.directFragCounts = adInfo.FragCounts;
                D.directFs = adInfo.ADFreq;
                clear adInfo;
            else
                assert(D.nDirectTime == dataInfo.AnalogChannels{i}.NumValues);
            end
        end
    end

    D.nDirectChAll = sum(D.directChannelIndices > 0);
    
    fprintf('Data file has %d Direct channels.\n', D.nDirectChAll);

    D.nDirectCh = numel(directChannelsToLoad);

    % save the data into D
    fprintf('\tProcessing %d Direct channels...\n', D.nDirectCh);
%     D.directs = nan(D.nDirectCh, D.nDirectTime);
    D.adjDirects = nan(D.nDirectCh, D.nDirectTime);
    D.directNames = cell(D.nDirectCh, 1);
    directChCounter = 0;
    for i = directChannelsToLoad
        channelInd = D.directChannelIndices(i);
        directChCounter = directChCounter + 1;
        adInfo = PL2Ad(fileName, dataInfo.AnalogChannels{channelInd}.Name);
%             D.directs(directChCounter,:) = adInfo.Values;
        adjDirect = padNaNsToAccountForDropsPL2(adInfo);
        if directChCounter == 1
            % adjust size of matrix after padding
            D.adjDirects = nan(D.nDirectCh, numel(adjDirect));
        end
        D.adjDirects(directChCounter,:) = adjDirect;

        D.directNames{directChCounter} = sprintf('%s_%s_%s', ...
                    sessionName, areaName, ...
                    dataInfo.AnalogChannels{channelInd}.Name);
    end
    clear adInfo;
    fprintf('\tdone.\n');
end