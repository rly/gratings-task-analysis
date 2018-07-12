%% read spike sorting quality metrics
sortQualityNotesFileName = 'spike sorting notes.xlsx';
xlsSheet = 1;
xlRange = 'A2:E3000';

[sortQualityNotesNums,sortQualityNotesSessions] = xlsread(sortQualityNotesFileName, xlsSheet, xlRange);
assert(size(sortQualityNotesNums, 2) == 4);
assert(size(sortQualityNotesNums, 1) == size(sortQualityNotesSessions, 1));
[sortQualityNotesNums,i] = trimNanRows(sortQualityNotesNums);
sortQualityNotesSessions(i) = [];

% sortQualityNotesNums(:,1) = channel number
% sortQualityNotesNums(:,2) = unit id (0 = unsorted)
% sortQualityNotesNums(:,3) = 0 or 1 if waveform has typical shape
% sortQualityNotesNums(:,4) = 0-5 quality of unit separation

%% load SUA
fileNames = dir('C:\Users\Ryan\Documents\MATLAB\gratings-task-data\MUA\*-SUA_*.mat');
nFiles = numel(fileNames);

wfs = cell(nFiles, 1); % prelalocate with underestimate
sessionNames = cell(nFiles, 1);
channelNum = nan(nFiles, 1);
unitIndInChannel = nan(nFiles, 1);
hasTypicalWaveformShape = false(nFiles, 1);
separationQuality = nan(nFiles, 1);
unitCount = 0;

for i = 1:nFiles
    if rem(i, 5) == 0 % print status
        fprintf('Processing files: %d/%d = %d%%\n', i, nFiles, round(i / nFiles * 100))
    end    
    
    filePath = sprintf('%s/%s', fileNames(i).folder, fileNames(i).name);
    L = load(filePath);
    
    sessionName = fileNames(i).name(1:strfind(fileNames(i).name, '-') - 1);
    channelInd = str2double(fileNames(i).name(strfind(fileNames(i).name, '_') + 1:end-4)); % indices start at 0
    suaData = load(filePath, sprintf('wfData%d', channelInd));
    suaData = suaData.(sprintf('wfData%d', channelInd));
    % suaData should have waveform x data where
    % suaData(:,1) = channel number
    % suaData(:,2) = unit id (0 = unsorted)
    % suaData(:,3) = timestamp in ms (time of wf minimum after threshold crossing)
    % suaData(:,4:nSamples+3) = waveform in microvolts
    
    nUnitsThisCh = max(suaData(:,2));
    for j = 1:nUnitsThisCh
        unitCount = unitCount + 1;
        unitMatch = suaData(:,2) == j;
        
        sessionNames{unitCount} = sessionName;
        channelNum(unitCount) = channelInd + 1;
        unitIndInChannel(unitCount) = j;
        wfs{unitCount} = suaData(unitMatch,4:end);
                
        qualityNotesMatchInd = strcmp(sortQualityNotesSessions, sessionName) & ...
                sortQualityNotesNums(:,1) == channelNum(unitCount) & ...
                sortQualityNotesNums(:,2) == unitIndInChannel(unitCount);
        if ~any(qualityNotesMatchInd)
            warning('No matches in SUA quality notes for session %s, channel %d, unit %d, file #%d\n', sessionName, channelNum(unitCount), j, i);
        elseif sum(qualityNotesMatchInd) > 1
            warning('Too many matches in SUA quality notes for session %s, channel %d, unit %d, file #%d\n', sessionName, channelNum(unitCount), j, i);
        end
        assert(sum(qualityNotesMatchInd) == 1);
        hasTypicalWaveformShape(unitCount) = (sortQualityNotesNums(qualityNotesMatchInd,3) == 1);
        separationQuality(unitCount) = sortQualityNotesNums(qualityNotesMatchInd,4);
    end
end
fprintf('Done processing files!\n');

%% load recording information
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
recordingInfo = readRecordingInfo(recordingInfoFileName);

%% save
% wfs variable is huge
saveFileName = 'allSuaData.mat';
save(saveFileName, 'recordingInfo', 'wfs', 'sessionNames', ...
        'channelNum', 'unitIndInChannel', 'hasTypicalWaveformShape', ...
        'separationQuality', 'unitCount', '-v7.3');
    
%% load
clear;
saveFileName = 'allSuaData.mat';
load(saveFileName);

%% compute mean wfs
meanWfs = cellfun(@mean, wfs, 'UniformOutput', false);
nUnits = unitCount;

%% determine whether each unit is in pulvinar
isInPulvinar = false(nUnits, 1);
for i = 1:nUnits
    sessionMatch = find(strcmp({recordingInfo.sessionName}, sessionNames{i}));
    for j = 1:numel(sessionMatch)
        R = recordingInfo(sessionMatch(j));
        if any(channelNum(i) == [R.dPulChannels R.vPulChannels])
            isInPulvinar(i) = 1;
        end
    end
end

%% summarize
nGoodSingleUnits = sum(isInPulvinar & hasTypicalWaveformShape & separationQuality >= 3);
fprintf('Found %d pulvinar units with typical waveform shape and good separation quality\n', nGoodSingleUnits);

nGreatSingleUnits = sum(isInPulvinar & hasTypicalWaveformShape & separationQuality >= 4);
fprintf('Found %d pulvinar units with typical waveform shape and great separation quality\n', nGreatSingleUnits);

nAwesomeSingleUnits = sum(isInPulvinar & hasTypicalWaveformShape & separationQuality == 5);
fprintf('Found %d pulvinar units with typical waveform shape and awesome separation quality\n', nAwesomeSingleUnits);

%% compute trough to peak time and plot histogram
peakInd = nan(nUnits, 1);
troughToPeakTime = nan(nUnits, 1);
spikeFs = 40000;
tToPSplitTime = 0.325;
nsCount = 0;
for i = 1:nUnits
    if isInPulvinar(i) && hasTypicalWaveformShape(i) && separationQuality(i) >= 3
        meanWf = meanWfs{i} / 1000; % mV
        [~,troughInd] = min(meanWf); % not necessarily the trough associated with thresh crossing
        [~,relPeakInd] = max(meanWf(troughInd:end)); % peak must be after trough
        peakInd(i) = relPeakInd + troughInd - 1;
        troughToPeakTime(i) = (relPeakInd-1)/spikeFs * 1000; % ms
        if troughToPeakTime(i) < tToPSplitTime
            nsCount = nsCount + 1;
            fprintf('NS #%d: Session %s, Channel %d, Unit %d\n', nsCount, sessionNames{i}, channelNum(i), unitIndInChannel(i));
        end
    end
end
troughToPeakTime(isnan(troughToPeakTime)) = [];

figure;
histogram(troughToPeakTime, 0:0.05:1);

%% plot all units considered as single units
t = (-16:40) / 40;
cols = lines(6);

figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'ML', 0.12, 'MB', 0.12);
hold on;
for i = 1:nUnits
    if isInPulvinar(i) && hasTypicalWaveformShape(i) && separationQuality(i) >= 3
        plot(t, meanWfs{i} / 1000, 'LineWidth', 1, 'Color', cols(rem(i, 6)+1,:));
        plot(t(peakInd(i)), meanWfs{i}(peakInd(i)) / 1000, '.', 'MarkerSize', 15, 'Color', cols(rem(i, 6)+1,:));
    end
end
plot(tToPSplitTime * [1 1], [-1 1], 'k--');
xlim(t([1 end]));
ylim([-0.25 0.17]);
box off;
grid on;
xlabel('Time from Trough (ms)');
ylabel('Voltage (mV)');
set(gca, 'XTick', -0.4:0.2:1);
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);
set(gcf, 'Color', 'w');

plotFileName = 'singleUnitWaveforms.png';
export_fig(plotFileName, '-nocrop');

%% plot PC space of mean waveforms
isGoodSingleUnit = isInPulvinar & hasTypicalWaveformShape & separationQuality >= 3;
meanWfsGood = meanWfs(isGoodSingleUnit);
meanWfsGood = cell2mat(meanWfsGood);
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca(meanWfsGood);
fprintf('PCA: %d variables, %d observations\n', size(meanWfsGood, 2), size(meanWfsGood, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC2 explains %0.1f%% of the variance.\n', pcaPctExplained(2));
fprintf('\tPC3 explains %0.1f%% of the variance.\n', pcaPctExplained(3));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', sum(pcaPctExplained(1:2)));
fprintf('\tPC1 + PC2 + PC3 explain %0.1f%% of the variance.\n', sum(pcaPctExplained(1:3)));

figure;
plot(-pcaScore(:,1), pcaScore(:,2), '.', 'MarkerSize', 20);

figure;
hold on;
plot(-pcaCoeff(:,1), 'LineWidth', 2);
plot(pcaCoeff(:,2), 'LineWidth', 2);
plot(pcaCoeff(:,3), 'LineWidth', 2);

% how to capture faster return to baseline on NS units?