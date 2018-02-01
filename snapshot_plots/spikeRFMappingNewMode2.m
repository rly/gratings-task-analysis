
clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';

%%
pl2FileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170529\d1-33.00mm-d2-32.00mm-all_merged_noWB-sort.pl2';
sessionName = 'M20170529';
areaName = 'PUL';
blockNames = {'vepm4', 'vepm5', 'aepm1', 'aepm2', 'aepm3', 'rfm1', 'rfm2', 'rfm3', 'rfm4', 'g1', 'g2', 'g3', 'g4', 'rfm5', 'rest1'};
rfmResultsRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170529\20170529\';
rfmResultsFiles = {'multi_flash_rf_mapping_results_20170529-173737.json', ...
    'multi_flash_rf_mapping_results_20170529-174741.json', ...
    'multi_flash_rf_mapping_results_20170529-175055.json', ...
    'multi_flash_rf_mapping_results_20170529-180028.json'};
blockInds = 7;

%%
pl2FileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170608\d1_34.50mm_d2_37.00mm_all_merged_noWB-sort.pl2';
sessionName = 'M20170608';
areaName = 'PUL';
blockNames = {'rest3', 'vepm1', 'aepm1', 'rfm1', 'rfm2', 'rfm3', 'rfm4', 'g1', 'g2', 'g3', 'g4', 'rfm5', 'vepm2', 'rest4'};
rfmResultsRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170608\20170608\';
rfmResultsFiles = {'multi_flash_rf_mapping_results_20170608-155719.json', ...
    'multi_flash_rf_mapping_results_20170608-160603.json', ...
    'multi_flash_rf_mapping_results_20170608-160844.json', ...
    'multi_flash_rf_mapping_results_20170608-161651.json'};
blockInds = 5;

%% setup and load data
fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Analysis - Spikes\n');
fprintf('Loading %s...\n', pl2FileName);

tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
D = loadPL2(pl2FileName, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc);

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(blockNames(blockInds), '-');

%% remove spike and event times not during RFM task to save memory
D = trimSpikeTimesAndEvents(D, blockInds);

%%
% task as of 1/10/17 (or earlier)
% 325ms fixation
% 25ms of pre-flashes period (continued fixation)
% up to N=4 flashes, each 100ms duration, separated by 200ms ISI
% 250-500ms after 5th flash until FP dims
% 260-775ms after FP dims is correct response window
% 4 distances to fixation 150:60:330
% 11 polar angles (-5:5)*pi/8
% 4 grating orientations (0:3)*pi/4
% 3 reps each

% EVT05 - begin pre-cue period
% EVT06 - flash onset
% EVT07 - FP dim
% EVT08 - juice onset

origFlashEvents = D.events{6};
preFlashesEvents = D.events{5};
nFlashes = numel(origFlashEvents);

% % update based on the PCL code
% % 4 possible, coded in events 9-10: 01, 01, 10, 11, in the order below
% distsToFix = [180 240 320];
% % 7 possible, coded in events 11-14: 0001, 0010, ..., 1011, in the order below
% polarAngles = (-5:5)*pi/8;
% % 4 possible, coded in events 15-16: 00, 01, 10, 11, in the order below
% gratingAngles = (0:3)*pi/4;
% 
% [flashOnsets,flashStats] = decodeFlashParams(origFlashEvents, ...
%         D.events, distsToFix, polarAngles, gratingAngles);
    
%%
trialStructs = loadjson([rfmResultsRoot rfmResultsFiles{2}]);

clear flashParams;
flashCount = 0;
for i = 1:numel(trialStructs)
    for j = 2:numel(trialStructs{i}.events)
        flashCount = flashCount + 1;
        if flashCount == 1
            flashParams = trialStructs{i}.events{j}.flashParams;
        else
            flashParams(flashCount) = trialStructs{i}.events{j}.flashParams;
        end
    end
end
distsToFix = cell2mat({flashParams.distToFix});
polarAngles = cell2mat({flashParams.polarAngle});
diameters = cell2mat({flashParams.diameter});
gratingAngles = cell2mat({flashParams.gratingAngle});

%         stimId: 34
%        distToFix: 180
%       polarAngle: 0.4490
%         diameter: 59
%     gratingAngle: 45
%         contrast: 0.5000
%        driftRate: 0

% flashStats(i,1:2) = distance(i) * [cos(theta(i)) sin(theta(i))];
% flashStats(i,3) = gratingAngles(gratingAngleBinaryCode + 1);
% flashStats(i,4:6) = [distBinaryCode + 1 polarAngleBinaryCode gratingAngleBinaryCode + 1];


%%
simultSignalTol = 0.001; % tolerance level for simultaneous signals --
% if two events are within this number of  seconds from each other, they 
% represent the same signal. 
flashDuration = 0.1;
flashDurationTol = 0.04;

allEventTimes = D.events;
flashEventTimes = D.events{6};
flashOnsets = nan(size(flashEventTimes));
stimIDs = nan(size(flashEventTimes));
for i = 1:numel(flashEventTimes)  
    isFoundEvt = 0;
    foundEvts = cell(16, 1);
    for j = 9:16
        foundEvts{j} = find(abs(flashEventTimes(i) - allEventTimes{j}) < simultSignalTol);
        % there should be at most 1 event9, event10, etc within the
        % simultaneous signal tolerance of each flashEventTime
        assert(numel(foundEvts{j}) <= 1)
        if ~isempty(foundEvts{j})
            isFoundEvt = 1;
        end
    end
    if ~isFoundEvt
        fprintf('huh...\n');
    end

    allFoundEvtTimes = [allEventTimes{9}(foundEvts{9}); 
            allEventTimes{10}(foundEvts{10}); 
            allEventTimes{11}(foundEvts{11}); 
            allEventTimes{12}(foundEvts{12}); 
            allEventTimes{13}(foundEvts{13});
            allEventTimes{14}(foundEvts{14}); 
            allEventTimes{15}(foundEvts{15}); 
            allEventTimes{16}(foundEvts{16})];
    if (isempty(allFoundEvtTimes))
        continue;
    end;
    % TODO: do pairwise comparisons instead of consecutive comparisons
    % but this should be good enough
    assert(all(diff(allFoundEvtTimes) < simultSignalTol));
    flashOnset = allFoundEvtTimes(1);
    
    % there should be a second event (in events 9-16) marking the flash
    % offset, indicating the trial was not aborted mid-flash
    isFoundEndEvt = 0;
    foundEndEvts = cell(16, 1);
    for j = 9:16
        foundEndEvts{j} = find(abs(allEventTimes{j} - flashOnset - flashDuration) < flashDurationTol);
        % there should be at most 1 event9, event10, etc within the
        % simultaneous signal tolerance of each flashEventTime
        assert(numel(foundEndEvts{j}) <= 1)
        if ~isempty(foundEndEvts{j})
            isFoundEndEvt = 1;
        end
    end
    if ~isFoundEndEvt
        fprintf('huhhuh...\n');
    end
    allFoundEndEvtTimes = [allEventTimes{9}(foundEndEvts{9}); 
            allEventTimes{10}(foundEndEvts{10}); 
            allEventTimes{11}(foundEndEvts{11}); 
            allEventTimes{12}(foundEndEvts{12}); 
            allEventTimes{13}(foundEndEvts{13});
            allEventTimes{14}(foundEndEvts{14}); 
            allEventTimes{15}(foundEndEvts{15}); 
            allEventTimes{16}(foundEndEvts{16})];
    if (isempty(allFoundEndEvtTimes))
        continue;
    end;
    % TODO: do pairwise comparisons instead of consecutive comparisons
    % but this should be good enough
    assert(all(diff(allFoundEndEvtTimes) < simultSignalTol));
    flashOffset = allFoundEndEvtTimes(1);
    
    flashOnsets(i) = flashOnset;
    stimIDs(i) = bi2de(cellfun(@isempty, foundEvts(9:16))' == 0);
end

flashOnsets(isnan(stimIDs)) = [];
stimIDs(isnan(stimIDs)) = [];

assert(all(stimIDs' == cell2mat({flashParams.stimId})));

%%
nUnits = numel(D.allSpikeStructs);
fprintf('Processing %d units...\n', nUnits);

for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));
        
    if isempty(spikeTimes)
        fprintf('\tNo spikes. Skipping...\n');
        continue;
    end

%% compute psth for each flash condition
% this does not subtract out per-condition baseline activity
psthWindow = [0.25 0.25]; % secs before, secs after
analysisWindowOffset = [0.02 0.22]; % secs offset from event onset
analysisWindow = analysisWindowOffset + psthWindow(1);
analysisWindowDur = diff(analysisWindowOffset);

baselineWindowOffset = [-0.2 0];
baselineWindow = baselineWindowOffset + psthWindow(1);
baselineWindowDur = diff(baselineWindowOffset);

kernelSigma = 0.01;
nTime = fix(5*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
psthT = linspace(0, sum(psthWindow), nTime);
t = psthT - psthWindow(1);

% for raster
allAlignedSpikeTimes = createdatamatpt(spikeTimes, flashOnsets, psthWindow);

% TODO consider removing quiet trials

distsToFixUnique = sort(unique(distsToFix));
polarAnglesUnique = sort(unique(polarAngles));

coordinatesIndex = zeros(size(flashOnsets));
coordinatesUnique = cell(1, 1);
coordinatesCount = 0;
for k = 1:numel(distsToFixUnique)
    for l = 1:numel(polarAnglesUnique)
        theseFlashOnsets = flashOnsets(distsToFix == distsToFixUnique(k) & polarAngles == polarAnglesUnique(l));
        if isempty(theseFlashOnsets)
            continue;
        end
        coordinatesCount = coordinatesCount + 1;
        coordinatesUnique{coordinatesCount} = distsToFixUnique(k) * [cos(polarAnglesUnique(l)) sin(polarAnglesUnique(l))];
        coordinatesIndex(distsToFix == distsToFixUnique(k) & polarAngles == polarAnglesUnique(l)) = coordinatesCount;
    end
end

gratingAnglesUnique = sort(unique(gratingAngles)); % degrees

allPsths = nan(numel(coordinatesUnique), numel(gratingAnglesUnique), nTime);
trialCounts = zeros(numel(coordinatesUnique), numel(gratingAnglesUnique));
responseMetric = zeros(numel(coordinatesUnique), numel(gratingAnglesUnique));
baselineMetric = zeros(numel(coordinatesUnique), numel(gratingAnglesUnique));

for k = 1:numel(coordinatesUnique)
    for l = 1:numel(gratingAnglesUnique)
        theseFlashOnsets = flashOnsets(coordinatesIndex' == k & gratingAngles == gratingAnglesUnique(l));
        trialCounts(k,l) = numel(theseFlashOnsets);
        if isempty(theseFlashOnsets)
            continue;
        end

        dataStruct = createdatamatpt(spikeTimes, theseFlashOnsets, psthWindow);

        theseNSpikesPostFlash = zeros(numel(dataStruct), 1);
        theseNSpikesBaseline = zeros(numel(dataStruct), 1);
        % maybe could just collapse the struct and then do the time check
        for n = 1:numel(dataStruct)
            theseNSpikesPostFlash(n) = sum(analysisWindow(1) <= dataStruct(n).times & ...
                dataStruct(n).times <= analysisWindow(2));
            theseNSpikesBaseline(n) = sum(baselineWindow(1) <= dataStruct(n).times & ...
                dataStruct(n).times <= baselineWindow(2));
        end
        responseMetric(k,l) = sum(theseNSpikesPostFlash) / numel(dataStruct); % average num spikes per trial
        baselineMetric(k,l) = sum(theseNSpikesBaseline) / numel(dataStruct); % average num spikes per trial

        dataStruct = createdatamatpt(spikeTimes, theseFlashOnsets, psthWindow);
        R = fixedPsth(dataStruct, kernelSigma, 2, psthT);
        if ~isempty(R)
            allPsths(k,l,:) = R;
        end
    end
end

%%
% figure TODOs
% add separate colorbars for big and for mini, or use same color axis?
% merge x axis for raw data
% label groups of raw data
% y tick off
f = figure_tr_inch(14, 8); clf
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

modTitle = sprintf('Spike RF Mapping: %s (%d flashes)', unitName, nFlashes);
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
waveformW = 0.1;
rasterByTimeW = 0.1;
rawW = 0.11;
distW = 0.63;
waveformH = 0.15;
rasterByTimeH = 0.32;
rawH = 0.375;
distH = rawH;

waveformLeft1 = 0.05;
rasterByTimeLeft1 = waveformLeft1;
rawLeft1 = rasterByTimeLeft1 + rasterByTimeW + 0.03;
rawLeft2 = rawLeft1;
distLeft1 = rawLeft1 + rawW + 0.05;
distLeft2 = distLeft1;
btm = 0.07;
loc1Btm = btm + rawH + 0.09;
infoTextTop = btm + 0.18;
rasterByTimeBtm = infoTextTop + 0.09;
waveformBtm = rasterByTimeBtm + rasterByTimeH + 0.1;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
plotSpikeWaveform(D, j);

%% plot raster
axes('Position', [rasterByTimeLeft1 rasterByTimeBtm rasterByTimeW rasterByTimeH]); 
hold on;

window = psthWindow;
data = allAlignedSpikeTimes;
rasterY = 0;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
lineHeight = 5;
for k = 1:numel(data)
    rasterY = rasterY + 1;
    if ~isempty(data(k).times)
        adjTimes = data(k).times - window(1);
        for l = 1:numel(data(k).times)
            plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
        end
    end
end

yBounds = [0 rasterY+1];
plot([0 0], yBounds, 'k');
plot([analysisWindowOffset(1) analysisWindowOffset(1)], yBounds, '-', 'Color', [0 0.7 0]);
plot([analysisWindowOffset(2) analysisWindowOffset(2)], yBounds, '-', 'Color', [0 0.7 0]);
xlim([-0.1 0.25]);
ylim(yBounds);
set(gca, 'YDir', 'reverse');
xlabel('Time from flash onset (s)'); 
ylabel('Flash number');
titleH = title('Raster by time', 'Interpreter', 'None');
set(gca, 'Layer', 'top');

%% info
writeUnitInfo(spikeStruct, axBig, -0.03, infoTextTop);

%% plot VES by grating orientation
for k = 1:numel(coordinatesUnique)
    if k == 1
        axes('Position', [rawLeft1 loc1Btm rawW rawH]);
    else
        axes('Position', [rawLeft1 btm rawW rawH]);
    end
    hold on;
    channelSep = 1;
    count1 = 1;
    yScale = 0.07;
    cols = lines(numel(gratingAnglesUnique));
    for l = 1:numel(gratingAnglesUnique)
        count1 = count1 + channelSep;
        plot(t, squeeze(allPsths(k,l,:))*yScale + count1, 'Color', cols(l,:));
    end
    count1 = count1 + 3;
    plot([0 0], [0 count1], '-', 'Color', [0.3 0.3 0.3]);
    fillH = jbfill(analysisWindowOffset, [-1 -1], [count1 count1], [0.3 1 0.3], ones(3, 1), 1, 0.1);
    uistack(fillH, 'bottom');
    xlim([-0.1 0.25]);
    ylim([1 count1]);
    set(gca, 'YTickLabel', {});
    if k == numel(coordinatesUnique)
        xlabel('Time from flash onset (s)');
    end
    title(sprintf('VES by orientation - Loc %d', k));
end

%% plot distribution of responses per orientation

% compute spike rate stats, not using psth
avgSpikeRatePostFlash = responseMetric / analysisWindowDur;
avgSpikeRateBaseline = baselineMetric / baselineWindowDur;
meanSpikeRate = mean(avgSpikeRatePostFlash(:));
sdSpikeRate = std(avgSpikeRatePostFlash(:));
maxSpikeRate = max(avgSpikeRatePostFlash(:));
baselineSpikeRate = mean(avgSpikeRateBaseline(:));

for k = 1:numel(coordinatesUnique)
    if k == 1
        distAx = axes('Position', [distLeft1 loc1Btm distW distH]);
    else
        distAx = axes('Position', [distLeft1 btm distW distH]);
    end
    
    bar(gratingAnglesUnique, avgSpikeRatePostFlash(k,:) - avgSpikeRateBaseline(k,:), 'FaceColor', lines(1));
    xlabel('Grating Orientation (deg)');
    ylabel('Spike Rate (Hz; baseline-corr)');
    title(sprintf('Orientation Tuning Curve - Location %d', k));
    set(gca, 'XTick', gratingAnglesUnique);
    diffGratingAngles = diff(gratingAnglesUnique);
    xlim(gratingAnglesUnique([1 end]) + [-0.75 0.75]*diffGratingAngles(1));
    box off;
    if k ~= numel(coordinatesUnique)
        xlabel('');
    end
end

plotFileName = sprintf('%s/%s-%s.png', processedDataDir, unitName, blockName);
export_fig(plotFileName, '-nocrop');
close;

end % for each unit