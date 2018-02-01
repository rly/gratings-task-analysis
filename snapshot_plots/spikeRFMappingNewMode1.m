
clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

% sessionInd = 4; blockInds = x; % until there is a better way
% rfmResultsFiles = {'multi_flash_rf_mapping_results_20170529-173737.json', ...
%     'multi_flash_rf_mapping_results_20170529-174741.json', ...
%     'multi_flash_rf_mapping_results_20170529-175055.json', ...
%     'multi_flash_rf_mapping_results_20170529-180028.json'};

sessionInd = 9; blockInds = 4; % until there is a better way
rfmResultsFiles = {'multi_flash_rf_mapping_results_20170608-155719.json'};

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% setup and load data
fprintf('\n-------------------------------------------------------\n');
fprintf('RF Mapping Analysis - Spikes\n');
fprintf('Loading %s...\n', pl2FileName);

tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 0;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(blockNames(blockInds), '-');

rfmResultsRoot = sprintf('%s/%s/%s/', dataDirRoot, sessionName, sessionName(2:end));

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
trialStructs = loadjson([rfmResultsRoot rfmResultsFiles{1}]);

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
allEventTimes{15} = -1;
allEventTimes{16} = -1; % TODO grating angle doesn't change. there are no EVT15, EVT16
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
        fprintf('huhhuh... %0.1f\n', flashEventTimes(i));
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



%% quick hack
fpStimIDs = cell2mat({flashParams.stimId});
fpCount = numel(fpStimIDs);
% for i = numel(stimIDs):-1:1
%     if stimIDs(i) == fpStimIDs(fpCount)
%         fpCount = fpCount - 1;
%         continue;
%     end
%     if stimIDs(i) == fpStimIDs(fpCount - 1)
%         flashParams(fpCount) = [];
%         fpCount = fpCount - 2;
%     end
% end

fprintf('Found %d flash event times and %d flashParam stimIDs.\n', ...
        numel(flashEventTimes), fpCount);

%%
% numel(stimIDs)
% numel(cell2mat({flashParams.stimId}))
% a = [stimIDs; 0];
% b = cell2mat({flashParams.stimId});
%% temp
% flashParams(1) = [];

%%
distsToFix = cell2mat({flashParams.distToFix});
polarAngles = cell2mat({flashParams.polarAngle});
diameters = cell2mat({flashParams.diameter});
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

allPsths = nan(numel(distsToFixUnique), numel(polarAnglesUnique), nTime);
trialCounts = zeros(numel(distsToFixUnique), numel(polarAnglesUnique));
responseMetric = zeros(numel(distsToFixUnique), numel(polarAnglesUnique));
baselineMetric = zeros(numel(distsToFixUnique), numel(polarAnglesUnique));

for k = 1:numel(distsToFixUnique)
    for l = 1:numel(polarAnglesUnique)
        theseFlashOnsets = flashOnsets(distsToFix == distsToFixUnique(k) & polarAngles == polarAnglesUnique(l));
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
f = figure_tr_inch(13, 7.5); clf
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
heatW = 0.6;
waveformH = 0.15;
rasterByTimeH = 0.32;
rawH = 0.84;
heatH = 0.84;

waveformLeft1 = 0.05;
rasterByTimeLeft1 = waveformLeft1;
rawLeft1 = rasterByTimeLeft1 + rasterByTimeW + 0.03;
rawLeft2 = rawLeft1 + rawW + 0.02;
heatLeft = rawLeft2 + rawW + 0.01;
btm = 0.07;
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

%% plot VES by eccentricity
axes('Position', [rawLeft1 btm rawW rawH]); 
hold on;
channelSep = 1;
groupSep = numel(polarAnglesUnique) / numel(distsToFixUnique) * 2;
count1 = 1;
yScale = 0.07;
cols = lines(numel(polarAnglesUnique));
for k = 1:numel(distsToFixUnique)
    for l = 1:numel(polarAnglesUnique)
        count1 = count1 + channelSep;
        plot(t, squeeze(allPsths(k,l,:))*yScale + count1, 'Color', cols(l,:));
    end
    count1 = count1 + groupSep;
end
count1 = count1 + 6;
plot([0 0], [0 count1], '-', 'Color', [0.3 0.3 0.3]);
fillH = jbfill(analysisWindowOffset, [-1 -1], [count1 count1], [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
xlim([-0.1 0.25]);
ylim([0 count1]);
set(gca, 'YTickLabel', {});
xlabel('Time from flash onset (s)');
title('VES by eccentricity');

%% plot VES by polar angle
axes('Position', [rawLeft2 btm rawW rawH]); 
hold on;
channelSep = 1;
groupSep = 2;
count2 = 1;
yScale = 0.07;
cols = lines(numel(polarAnglesUnique));
for l = 1:numel(polarAnglesUnique)
    for k = 1:numel(distsToFixUnique)
        count2 = count2 + channelSep;
        plot(t, squeeze(allPsths(k,l,:))*yScale + count2, 'Color', cols(l,:));
    end
    count2 = count2 + groupSep;
end
count2 = count2 + 6;
plot([0 0], [0 count2], '-', 'Color', [0.3 0.3 0.3]);
fillH = jbfill(analysisWindowOffset, [-1 -1], [count2 count2], [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
xlim([-0.1 0.25]);
ylim([0 count2]);
set(gca, 'YTickLabel', {});
% xlabel('Time from flash onset (s)');
title('VES by polar angle');
assert(count1 == count2); % i.e. the y axes are the same for the two raw plots

%% make heatmap of average number of spikes (doesn't use psth)
mapScale = 1/10;
mapDim = round(800 * mapScale) + 1; % x and y: -400 to 400 
rfmapByOri = zeros(numel(distsToFix), mapDim, mapDim); % extra dim uses all trials
rfmapByOriCount = zeros(numel(distsToFix), mapDim, mapDim);
rfmapSmoothByOri = zeros(size(rfmapByOri));
mapXOffset = (mapDim + 1)/2;
mapYOffset = (mapDim + 1)/2;
numPixelsPerDegree = 28;

meanAllResponses = mean(responseMetric(:));
baselineResponse = mean(baselineMetric(:));
threshold = baselineResponse;
thresholdName = 'baseline';

gaussianFilterExtent = mapDim;
for k = 1:numel(distsToFixUnique)
    % each distToFix has a unique diameter
    matchingDiameters = diameters(distsToFix == distsToFixUnique(k));
    assert(all(matchingDiameters == matchingDiameters(1)))
    diskDiameter = matchingDiameters(1) * mapScale;
    
    uniformFilter = fspecial('disk', diskDiameter/2);
    uniformFilter = uniformFilter/max(uniformFilter(:));
    
    gaussianFilterSigma = diskDiameter / 2;
    gaussianFilter = fspecial('gaussian', [gaussianFilterExtent gaussianFilterExtent], ...
            gaussianFilterSigma);
    gaussianFilter = gaussianFilter / max(gaussianFilter(:));
    
    % better way: gaussian filter after uniform filter
    
%     filterExtent = diskDiameter + floor(68 * mapScale); 
%     noisyEdgeUniformFilter = zeros(filterExtent, filterExtent);
%     startUniformFilterInd = round((filterExtent-size(uniformFilter, 1))/2) + 1;
%     uniformFilterInds = startUniformFilterInd:startUniformFilterInd+size(uniformFilter, 1)-1;
%     noisyEdgeUniformFilter(uniformFilterInds,uniformFilterInds) = uniformFilter;
    
%     gaussianFilter = fspecial('gaussian', [34 34]*2, 34);
%     noisyEdgeUniformFilter2 = conv2(, gaussianFilter, 'same');
%     noisyEdgeUniformFilter2(noisyEdgeUniformFilter == 1) = 0;
%     noisyEdgeUniformFilter2 = noisyEdgeUniformFilter2 / max(noisyEdgeUniformFilter2(:));
%     noisyEdgeUniformFilter3 = noisyEdgeUniformFilter + noisyEdgeUniformFilter2;
%     noisyEdgeUniformFilter3(noisyEdgeUniformFilter3 > 1) = 1;
%     
%     figure; imagesc(noisyEdgeUniformFilter3); colorbar;
%     figure; plot(noisyEdgeUniformFilter3(floor(size(noisyEdgeUniformFilter3, 1)/2)+1,:))
    
    for l = 1:numel(polarAnglesUnique)
        flashX = distsToFixUnique(k) * cos(polarAnglesUnique(l));
        flashY = distsToFixUnique(k) * sin(polarAnglesUnique(l));
        mapX = round(flashX * mapScale) + mapXOffset;
        mapY = round(flashY * mapScale) + mapYOffset;
        if trialCounts(k,l) > 0
            % weighted mean across orientations
            % threshold at baseline
            rfmapByOri(k,mapX,mapY) = max(0, responseMetric(k,l,:) - threshold);
            rfmapByOriCount(k,mapX,mapY) = trialCounts(k,l); 
        end
    end
    rfmapSmoothByOri(k,:,:) = imfilter(squeeze(rfmapByOri(k,:,:)), gaussianFilter, 'replicate');
end

rfmapAllOri = squeeze(sum(rfmapByOri, 1));
rfmapAllOriCount = squeeze(sum(rfmapByOriCount, 1));
rfmapAllSmoothByOriAll = squeeze(sum(rfmapSmoothByOri(:,:,:), 1));


%% plot smoothed heatmap for all and individual orientations

% compute spike rate stats, not using psth
avgSpikeRatePostFlash = responseMetric / analysisWindowDur;
avgSpikeRateBaseline = baselineMetric / baselineWindowDur;
meanSpikeRate = mean(avgSpikeRatePostFlash(:));
sdSpikeRate = std(avgSpikeRatePostFlash(:));
maxSpikeRate = max(avgSpikeRatePostFlash(:));
baselineSpikeRate = mean(avgSpikeRateBaseline(:));

% start with all orientations together
heatAx = axes('Position', [heatLeft btm heatW heatH]); 
hold on;
plotRfMapSmooth(rfmapAllSmoothByOriAll, rfmapAllOriCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset);
title(sprintf('VES Heatmap (%d-%d ms after flash onset; all orientations)', round(analysisWindowOffset * 1000)), 'Interpreter', 'none');
textParams = {'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]};
text(0.02, 0.08, {sprintf('Threshold: %s', thresholdName), ...
        sprintf('Baseline rate: %0.2f Hz', baselineSpikeRate), ...
        sprintf('Mean rate: %0.2f Hz', meanSpikeRate), ...
        sprintf('Std rate: %0.2f Hz', sdSpikeRate), ...
        sprintf('Max rate: %0.2f Hz (%0.1f SDs)', ...
        maxSpikeRate, (maxSpikeRate - meanSpikeRate) / sdSpikeRate)}, ...
        textParams{:});
axis(heatAx, 'square'); % shouldn't do much if i set the dims properly
colorbar;
caxis([0 5]); 

plotFileName = sprintf('%s/%s-%s.png', processedDataDir, unitName, blockName);
export_fig(plotFileName, '-nocrop');
close;

end % for each unit