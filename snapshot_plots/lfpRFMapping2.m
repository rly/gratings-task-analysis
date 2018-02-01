% NEW RF mapping M20170529+

%%
clear;

lfpsFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170529\d1-33.00mm-d2-32.00mm-all_merged_noWB.pl2';
animalName = 'M';
sessionName = '20170529';
areaName = 'PUL';
blockNames = {'vepm4', 'vepm5', 'aepm1', 'aepm2', 'aepm3', 'rfm1', 'rfm2', 'rfm3', 'rfm4', 'g1', 'g2', 'g3', 'g4', 'rfm5', 'rest1'};
blockInds = 6;
rfmMode = 1;

logDir = sprintf('%s\\%s\\', fileparts(lfpsFileName), sessionName);
logFiles = {...
        'multi_flash_rf_mapping_params_20170529-173737.txt', ... 
        'multi_flash_rf_mapping_params_20170529-174741.txt', ...
        'multi_flash_rf_mapping_params_20170529-175055.txt', ...
        'multi_flash_rf_mapping_params_20170529-180028.txt', ...
        'multi_flash_rf_mapping_params_20170529-194452.txt'};

%%
fprintf('\n-------------------------------------------------------\n');
fprintf('RF Mapping Analysis - LFPs\n');
fprintf('Loading %s...\n', lfpsFileName);
tic;
isLoadSpikes = 1;
isLoadLfp = 1;
isLoadSpkc = 0;
D = loadPL2(lfpsFileName, sessionName, animalName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc); 

processedDataDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/%s%s', animalName, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(blockNames(blockInds), '-');

%% remove spike and event times not during RFM task to save memory
D = trimSpikeTimesAndEvents(D, blockInds);

%%
paramsLogFileName = sprintf('%s\\%s', logDir, logFiles{1});
fid = fopen(paramsLogFileName, 'r');

while 1
    tline = fgetl(fid);
    if startsWith(tline, 'stimuli_filename')
        break;
    end
end
% stimuli[1]=[55,180,2.468,59,45,0.5,0]
% {1} = stimulus order
% {2} = stimulus ID
% {3} = eccentricity
% {4} = polar angle (radians)
% {5} = diameter
% {6} = gratings orientation (degrees)
% {7} = contrast
% {8} = drift rate
stimuliData = textscan(fid, 'stimuli[%d]=[%d,%f,%f,%f,%f,%f,%f]');
fclose(fid);

if rfmMode == 1
    numStimuli = 684;
end
assert(numel(stimuliData{1}) == numStimuli);

distsToFix = unique(stimuliData{3});
polarAngles = unique(stimuliData{4});
gratingAngles = unique(stimuliData{6});

%%
% EVT05 - begin pre-cue period
% EVT06 - flash onset
% EVT07 - FP dim
% EVT08 - juice onset

doOutlierCheckPlot = 0;

periFlashWindowOffset = [-0.25 0.3]; % ms around flash
baselineWindowOffset = [-0.2 0]; % ms - 200ms minimum between flashes and before first flash
expandedPlotWindowOffset = [-0.3 1.6]; % ms - 150ms minimum between flashes and before first flash

origFlashEvents = D.events{6};
preFlashesEvents = D.events{5};
nFlashes = numel(origFlashEvents);

% currently no outlier detection
maxAbsYNonOutlier = 1000;%0.25;
outlierMaxSD = 12;
outlierCheckWindowOffset = periFlashWindowOffset; 

%% check codes
[flashOnsets,decimalCodes] = decodeFlashParams2(origFlashEvents, D.events);

for i = 1:numel(flashOnsets)
    assert(decimalCodes(i) == stimuliData{2}(mod(i-1, numStimuli) + 1));
end

%%
% assert(numel(D.fragTs) == 1);
% origFlashEvents = D.events{6} - D.fragTs(1); % adjust so that fragTs is not needed - just index
% preFlashesEvents = D.events{5} - D.fragTs(1);
% flashOnsets = flashOnsets - D.fragTs(1);
% nFlashes = numel(flashOnsets);
Fs = D.lfpFs;
nChannels = D.nLfpCh;

plotFileNamePrefix = sprintf('%s%s_%s_FPall', animalName, sessionName, areaName);

%% low-pass filter data even further
% ideally don't low-pass at 300 and then again at another freq
hiCutoffFreqCustom = 100; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
D.adjLfps(isnan(D.adjLfps)) = 0;
D.adjLfpsLP = filtfilt(bFirLowPassCustom, 1, D.adjLfps');
D.adjLfpsLP = D.adjLfpsLP';
D.adjLfps = [];

%% calculate and plot per-trial baseline
fPlot = figure_tr_inch(6, 5);
[tBaseline,averageBaselineResponses] = computeResponsesInWindow(D.adjLfpsLP, preFlashesEvents, ...
        baselineWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaseline, averageBaselineResponses, 'yScale', 1000);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baseline.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot baseline responses overlaid
figure_tr_inch(6, 5, 6, 0);
hold on;
plot(tBaseline, averageBaselineResponses - repmat(mean(averageBaselineResponses, 2), 1, size(averageBaselineResponses, 2)));
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baselineOverlaid.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot expanded range around trial onset
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesExp] = computeResponsesInWindow(D.adjLfpsLP, preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesExp, 'yScale', 100);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot expanded range around trial onset after cross-channel CAR
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesCARExp] = computeResponsesInWindow(D.adjLfpsLP - nanmean(D.adjLfpsLP, 1), preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesCARExp, 'yScale', 500);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel (CAR)');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineCARExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% detect events with outlier activity (e.g. amp saturated)
outlierCheckStartIndices = round((origFlashEvents + outlierCheckWindowOffset(1)) * Fs); 
outlierCheckEndIndices = outlierCheckStartIndices + round(diff(outlierCheckWindowOffset) * Fs) - 1;
outlierCheckT = outlierCheckWindowOffset(1):1/Fs:outlierCheckWindowOffset(2)-1/Fs;
nOutlierCheckTime = numel(outlierCheckT);
assert(all(nOutlierCheckTime == (outlierCheckEndIndices - outlierCheckStartIndices + 1)));

isFlashEventOutlier = zeros(nFlashes, 1);
fprintf('Checking for flashes with outlier activity (e.g. amp saturated)...\n');
for j = 1:nChannels
    fprintf('Processing %s...\n', D.lfpNames{j});
    
    % normalize each channel by its mean and standard deviation
    % units are now standard deviations away from the mean
    channelDataCAR = D.adjLfpsLP(j,:) - nanmean(D.adjLfpsLP, 1);
    channelDataNorm = (channelDataCAR - nanmean(channelDataCAR)) / nanstd(channelDataCAR);
    
    % plot per channel responses before outlier removal
    if doOutlierCheckPlot
        fHandles(j) = figure_tr_inch(12,6);
        suptitle(sprintf('Outlier Check on Channel %d', j));
        subaxis(1,2,1);
        hold on;
    end
    
    for i = 1:nFlashes
        dataInOutlierCheckPeriod = channelDataNorm(outlierCheckStartIndices(i):outlierCheckEndIndices(i));
        % mark outlier flashes
        if any(isnan(dataInOutlierCheckPeriod))
            isFlashEventOutlier(i) = 1;
            fprintf('\tDetected %d outlier flashes because of NaN in trial\n', ...
                    sum(isnan(dataInOutlierCheckPeriod)));
        elseif any(abs(dataInOutlierCheckPeriod) > maxAbsYNonOutlier)
            isFlashEventOutlier(i) = 1;
            fprintf('\tDetected %d outlier flashes because absolute value of voltage in trial is too high (> %0.1f)\n', ...
                    sum(abs(dataInOutlierCheckPeriod) > maxAbsYNonOutlier), maxAbsYNonOutlier);
        elseif any(abs(dataInOutlierCheckPeriod) > outlierMaxSD)
            isFlashEventOutlier(i) = 1;
            fprintf('\tDetected %d outlier flashes because voltage in trial is too extreme (> %d SDs)\n', ...
                    sum(abs(dataInOutlierCheckPeriod) > outlierMaxSD), outlierMaxSD);
        end
    
        if doOutlierCheckPlot
            plot(outlierCheckT, dataInOutlierCheckPeriod);
            xlim(outlierCheckWindowOffset);
        end
    end
    
    if doOutlierCheckPlot
        title('Original Responses');
        xlabel('Time from Flash Onset');
        ylabel('SDs from Overall Mean');
        origYLim = ylim();
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        ylim(origYLim);
    end
end
clear channelDataNorm;

% remove outlier flashes
fprintf('\tDetected %d outlier flashes. Removing...\n', sum(isFlashEventOutlier));
flashEvents = origFlashEvents;
flashEvents(isFlashEventOutlier == 1) = []; % TODO look into flashOnsets vs flashEvents
nFlashes = numel(flashOnsets);

%% plot per channel responses after outlier removal

% recompute baseline (?)

startIndices = round((flashOnsets + periFlashWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periFlashWindowOffset) * Fs) - 1;
t = periFlashWindowOffset(1):1/Fs:periFlashWindowOffset(2)-1/Fs;
nTime = numel(t);
assert(all(nTime == (endIndices - startIndices + 1)));

responses = nan(nChannels, nTime, nFlashes);
channelDataNorm = nan(size(D.adjLfpsLP));
for j = 1:nChannels
    % normalize each channel by mean baseline across trials and time and
    % std across time of mean baseline across trials
    % TODO: bar plot of average baseline response and std by channel
    fprintf('Processing %s...\n', D.lfpNames{j});
    
%     channelDataBaselined{j} = (channelData{j} - nanmean(averageBaselineResponses(j,:))) / nanstd(averageBaselineResponses(j,:));
    channelDataCAR = D.adjLfpsLP(j,:) - nanmean(D.adjLfpsLP, 1);
    channelDataNorm(j,:) = (channelDataCAR - nanmean(channelDataCAR)) / nanstd(channelDataCAR);
    
    % can vectorize???
    for i = 1:nFlashes
        responses(j,:,i) = channelDataNorm(j,startIndices(i):endIndices(i));
    end
end

% plot for sanity check
if doOutlierCheckPlot
    for j = 1:nChannels
        figure(fHandles(j));
        subaxis(1,2,2);
        hold on;
        plot(t, squeeze(responses(j,:,:)));
        xlim(outlierCheckWindowOffset);
        title('After Outlier Removal');
        xlabel('Time from Flash Onset');
        origYLim = ylim();
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        ylim(origYLim);
    end
end

%% process RF for each channel
for j = 1:nChannels

channelName = D.lfpNames{j};
fprintf('Processing %s (%d/%d = %d%%)... \n', channelName, j, ...
        nChannels, round(j/nChannels*100));
lfpAdj = channelDataNorm(j,:);
responsesThisCh = squeeze(responses(j,:,:));

%%
periEventWindowOffset = [-0.25 0.25];
t = periEventWindowOffset(1):1/Fs:periEventWindowOffset(2)-1/Fs;

analysisWindowOffset = [0.025 0.2];
analysisWindowLogical = t >= analysisWindowOffset(1) & t < analysisWindowOffset(2);

baselineWindowOffset = [-0.2 0];
baselineWindowLogical = t >= baselineWindowOffset(1) & t < baselineWindowOffset(2);

averageFlashResponse = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles)+1, numel(t));
responseMetric = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles)+1);
baselineMetric = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles)+1);
trialCounts = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles));
m2 = numel(gratingAngles) + 1; % for indexing the no orientation case

for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles)
        for m = 1:numel(gratingAngles)
            theseFlashOnsets = flashOnsets(stimuliData{3} == distsToFix(k) & stimuliData{4} == polarAngles(l) & stimuliData{6} == gratingAngles(m));
            trialCounts(k,l,m) = numel(theseFlashOnsets);
            if isempty(theseFlashOnsets)
                continue;
            end
                        
            rawSignals = nan(numel(theseFlashOnsets), numel(t));
            % convert flash time to lfp variable index
            % slightly more accurate than createdatamatc b/c of rounding after offset
            startIndices = round((theseFlashOnsets + periEventWindowOffset(1)) * Fs); 
            endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;
            
            for n = 1:numel(theseFlashOnsets)
                rawSignals(n,:) = lfpAdj(startIndices(n):endIndices(n));
            end
            averageFlashResponse(k,l,m,:) = mean(rawSignals, 1); % mean over flashes
%             averageFlashResponseMeanOverTime(k,l,m) = mean(abs(averageFlashResponse(k,l,m,:)));
            responseMetric(k,l,m) = max(averageFlashResponse(k,l,m,analysisWindowLogical)) - ...
                    min(averageFlashResponse(k,l,m,analysisWindowLogical));
            baselineMetric(k,l,m) = max(averageFlashResponse(k,l,m,baselineWindowLogical)) - ...
                    min(averageFlashResponse(k,l,m,baselineWindowLogical));
        end
        
        % case: regardless of orientation 
        theseFlashOnsets = flashOnsets(stimuliData{3} == distsToFix(k) & stimuliData{4} == polarAngles(l));
        if isempty(theseFlashOnsets)
            continue;
        end
        
        rawSignals = nan(numel(theseFlashOnsets), numel(t));
        % convert flash time to lfp variable index
        % slightly more accurate than createdatamatc b/c of rounding after offset
        startIndices = round((theseFlashOnsets + periEventWindowOffset(1)) * Fs); 
        endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;
        
        for n = 1:numel(theseFlashOnsets)
            rawSignals(n,:) = lfpAdj(startIndices(n):endIndices(n));
        end
        averageFlashResponse(k,l,m2,:) = mean(rawSignals, 1); % mean over flashes
%         averageFlashResponseMeanOverTime(k,l,m2) = mean(abs(averageFlashResponse(k,l,m2,:)));
        responseMetric(k,l,m2) = max(averageFlashResponse(k,l,m2,analysisWindowLogical)) - ...
                min(averageFlashResponse(k,l,m2,analysisWindowLogical));
        baselineMetric(k,l,m2) = max(averageFlashResponse(k,l,m2,baselineWindowLogical)) - ...
                    min(averageFlashResponse(k,l,m2,baselineWindowLogical));
    end
end

rawSignals = nan(numel(flashOnsets), numel(t));
% convert flash time to lfp variable index
% slightly more accurate than createdatamatc b/c of rounding after offset
startIndices = round((flashOnsets + periEventWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;

for n = 1:numel(flashOnsets)
    rawSignals(n,:) = lfpAdj(startIndices(n):endIndices(n));
end

%%
% figure TODOs
% add separate colorbars for big and for mini, or use same color axis?
% merge x axis for raw data
% label groups of raw data
% y tick off
f = figure_tr_inch(14, 8); clf
set(gcf,'Color','white');
set(gcf,'renderer','painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig,'Title'), 'Visible', 'on')

modTitle = sprintf('%s (%d flashes)', channelName, nFlashes);
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
waveformW = 0.1;
rasterByTimeW = 0.1;
rawW = 0.11;
heatW = 0.6;
% miniHeatW = 0.18;
waveformH = 0.15;
rasterByTimeH = 0.35;
rawH = 0.85;
heatH = 0.85;
% miniHeatH = 0.385;

waveformLeft1 = 0.05;
rasterByTimeLeft1 = waveformLeft1;
rawLeft1 = rasterByTimeLeft1 + rasterByTimeW + 0.03;
rawLeft2 = rawLeft1 + rawW + 0.02;
heatLeft = rawLeft2 + rawW;
% miniHeatLeft1 = heatLeft + heatW + 0.005;
% miniHeatLeft2 = miniHeatLeft1 + miniHeatW + 0.01;
btm = 0.07;
rasterByTimeBtm = btm + 0.25;
waveformBtm = rasterByTimeBtm + rasterByTimeH + 0.1;
% miniHeatBtm2 = btm + miniHeatH + 0.08;

%%
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
hold on;


%%
axes('Position', [rasterByTimeLeft1 rasterByTimeBtm rasterByTimeW rasterByTimeH]); 
hold on;

nRawToShow = 20;
indRawToShow = round(linspace(1, size(rawSignals, 1), nRawToShow));
rasterY = 0;
yScale = 0.5;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
for k = 1:numel(indRawToShow)
    rasterY = rasterY + 1;
    plot(t, rawSignals(indRawToShow(k),:)*yScale + rasterY, lineParams{:});
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
titleH = title('Example LFPs by time', 'Interpreter', 'None');
set(gca, 'Layer', 'top');


%% info
textParams = {'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'none'};
text(axBig, -0.03, 0.11, {''}, ...
        textParams{:});

%%
axes('Position', [rawLeft1 btm rawW rawH]); 
hold on;
channelSep = 1;
groupSep = (numel(polarAngles) - 1) / (numel(distsToFix) - 1) * 2;
count1 = 1;
yScale = 1;
cols = lines(numel(polarAngles));
for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles)
        count1 = count1 + channelSep;
        plot(t, squeeze(averageFlashResponse(k,l,m2,:))*yScale + count1, 'Color', cols(l,:));
    end
    if k < numel(distsToFix)
        count1 = count1 + groupSep;
    end
end
count1 = count1 + 2;
plot([0 0], [0 count1], '-', 'Color', [0.3 0.3 0.3]);
plot([analysisWindowOffset(1) analysisWindowOffset(1)], [0 count1], '-', 'Color', [0 0.7 0]);
plot([analysisWindowOffset(2) analysisWindowOffset(2)], [0 count1], '-', 'Color', [0 0.7 0]);
xlim([-0.1 0.25]);
ylim([0 count1]);
set(gca, 'YTickLabel', {});
xlabel('Time from flash onset (s)');
title('VEP by eccentricity');

%%
axes('Position', [rawLeft2 btm rawW rawH]); 
hold on;
channelSep = 1;
groupSep = 2;
count2 = 1;
yScale = 1;
cols = lines(numel(polarAngles));
for l = 1:numel(polarAngles)
    for k = 1:numel(distsToFix)
        count2 = count2 + channelSep;
        plot(t, squeeze(averageFlashResponse(k,l,m2,:))*yScale + count2, 'Color', cols(l,:));
    end
    if l < numel(polarAngles)
        count2 = count2 + groupSep;
    end
end
count2 = count2 + 2;
plot([0 0], [0 count2], '-', 'Color', [0.3 0.3 0.3]);
plot([analysisWindowOffset(1) analysisWindowOffset(1)], [0 count2], '-', 'Color', [0 0.7 0]);
plot([analysisWindowOffset(2) analysisWindowOffset(2)], [0 count2], '-', 'Color', [0 0.7 0]);
xlim([-0.1 0.25]);
ylim([0 count2]);
set(gca, 'YTickLabel', {});
% xlabel('Time from flash onset (s)');
title('VEP by polar angle');
assert(count1 == count2); % i.e. the y axes are the same for the two raw plots

%% make heatmap
mapScale = 1/10;
mapDim = round(800 * mapScale) + 1; % x and y: -400 to 400 
rfmapByOri = zeros(numel(distsToFix), numel(gratingAngles)+1, mapDim, mapDim); % extra dim uses all trials
rfmapByOriCount = zeros(numel(distsToFix), numel(gratingAngles)+1, mapDim, mapDim);
rfmapSmoothByOri = zeros(size(rfmapByOri));
mapXOffset = (mapDim + 1)/2;
mapYOffset = (mapDim + 1)/2;
numPixelsPerDegree = 28;

meanAllResponses = mean(responseMetric(:));
baselineResponse = mean(baselineMetric(:));
threshold = baselineResponse;
thresholdName = 'baseline';

gaussianFilterExtent = mapDim;
for k = 1:numel(distsToFix)
    % slight adjustment for eve movement, which is radius 68 but here let's treat it as diameter 68
    diskDiameter = floor((distsToFix(k) / 150 * 3 * 28 * 2 - 68) * mapScale); % pixels, floor to nearest even
    uniformFilter = fspecial('disk', diskDiameter/2);
    uniformFilter = uniformFilter/max(uniformFilter(:));
    
    gaussianFilterSigma = diskDiameter / 3;
    gaussianFilter = fspecial('gaussian', [gaussianFilterExtent gaussianFilterExtent], ...
            gaussianFilterSigma);
    gaussianFilter = gaussianFilter / max(gaussianFilter(:));
    
    for l = 1:numel(polarAngles)
        flashX = distsToFix(k) * cos(polarAngles(l));
        flashY = distsToFix(k) * sin(polarAngles(l));
        mapX = round(flashX * mapScale) + mapXOffset;
        mapY = round(flashY * mapScale) + mapYOffset;
        for m = 1:numel(gratingAngles)
            % threshold at baseline
            rfmapByOri(k,m,mapX,mapY) = max(0, responseMetric(k,l,m) - threshold);
            rfmapByOriCount(k,m,mapX,mapY) = trialCounts(k,l,m); 
        end
%         if sum(trialCounts(k,l,:)) > 0
%             % weighted mean across orientations
%             % threshold at baseline
%             rfmapByOri(k,m2,mapX,mapY) = max(0, sum(trialCounts(k,l,:).*responseMetric(k,l,:)) ...
%                     / sum(trialCounts(k,l,:)) - threshold);
%             rfmapByOriCount(k,m2,mapX,mapY) = sum(trialCounts(k,l,:)); 
%         end
        rfmapByOri(k,m2,mapX,mapY) = max(0, responseMetric(k,l,m2) - threshold);
        rfmapByOriCount(k,m2,mapX,mapY) = sum(trialCounts(k,l,:)); 
    end
    % maintain separate smooth maps for each distsToFix
    for m = 1:numel(gratingAngles) + 1
        rfmapSmoothByOri(k,m,:,:) = imfilter(squeeze(rfmapByOri(k,m,:,:)), gaussianFilter, 'replicate');
    end
end

rfmapAllOri = squeeze(sum(rfmapByOri(:,numel(gratingAngles)+1,:,:), 1));
rfmapAllOriCount = squeeze(sum(rfmapByOriCount(:,numel(gratingAngles)+1,:,:), 1));

rfmapSmoothByOriAll = squeeze(sum(rfmapSmoothByOri, 1)); % sum over the different smooth maps for each distsToFix
rfmapByOriCountAll = squeeze(sum(rfmapByOriCount, 1)); % sum over the different smooth maps for each distsToFix
rfmapAllSmoothByOriAll = squeeze(sum(rfmapSmoothByOri(:,numel(gratingAngles)+1,:,:), 1));

%%
% flatFilterSize = 4;
% flatFilter = ones(flatFilterSize, flatFilterSize);
% figure;
% hold on;
% rfmapSmooth = conv2(rfmapAllOri, flatFilter, 'same');
% imagesc(((1:mapDim) - mapXOffset) / numPixelsPerDegree, ...
%         ((1:mapDim) - mapYOffset) / numPixelsPerDegree, rfmapSmooth');
% xlim(([1 mapDim] - mapXOffset) / numPixelsPerDegree);
% ylim(([1 mapDim] - mapYOffset) / numPixelsPerDegree);
% xlabel('Degrees horizontal');
% ylabel('Degrees vertical');
% colorbar;
% 
% title(sprintf('%s - Response Heatmap (20-220ms after flash onset)', spkVarName), 'Interpreter', 'none');

%% plot smoothed heatmap for all and individual orientations

% compute lfp stats
sdAllResponses = std(responseMetric(:));
maxAllResponses = max(responseMetric(:));

% start with all orientations together
heatAx = axes('Position', [heatLeft btm heatW heatH]); 
hold on;
plotRfMapSmooth(rfmapAllSmoothByOriAll, rfmapAllOriCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset);
title(sprintf('VEP Heatmap (%d-%d ms after flash onset; all orientations)', round(analysisWindowOffset * 1000)), 'Interpreter', 'none');
textParams = {'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]};
text(0.02, 0.1, sprintf('Threshold: %s', thresholdName), 'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]);
text(0.02, 0.08, sprintf('Baseline amp: %0.2f', baselineResponse), 'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]);
text(0.02, 0.06, sprintf('Mean amp: %0.2f', meanAllResponses), 'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]);
text(0.02, 0.04, sprintf('Std amp: %0.2f', sdAllResponses), 'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]);
text(0.02, 0.02, sprintf('Max amp: %0.2f (%0.1f SDs)', ...
        maxAllResponses, (maxAllResponses - meanAllResponses) / sdAllResponses), ...
        textParams{:});
axis(heatAx, 'square'); % shouldn't do much if i set the dims properly
colorbar;
caxis([0 6*sdAllResponses]); 

% % iterate through all orientations
% maxRfmapSmooth = 0;
% miniHeatAx = nan(numel(gratingAngles), 1);
% for m = 1:numel(gratingAngles)
%     if m == 1
%         miniHeatAx(m) = axes('Position', [miniHeatLeft1 miniHeatBtm2 miniHeatW miniHeatH]); 
%         oriText = '-';
%     elseif m == 2
%         miniHeatAx(m) = axes('Position', [miniHeatLeft2 miniHeatBtm2 miniHeatW miniHeatH]); 
%         oriText = '/';
%     elseif m == 3
%         miniHeatAx(m) = axes('Position', [miniHeatLeft1 btm miniHeatW miniHeatH]); 
%         oriText = '|';
%     elseif m == 4
%         miniHeatAx(m) = axes('Position', [miniHeatLeft2 btm miniHeatW miniHeatH]); 
%         oriText = '\';
%     else
%         error('unknown m');
%     end
%     hold on;
%     
%     rfmapSmoothThisOri = squeeze(rfmapSmoothByOriAll(m,:,:));
%     rfmapThisOriCount = squeeze(rfmapByOriCountAll(m,:,:));
% 
%     plotRfMapSmooth(rfmapSmoothThisOri, rfmapThisOriCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset);
% 
%     title(sprintf('Orientation: %0.2f %s', gratingAngles(m), oriText), 'Interpreter', 'none');
%     axis(miniHeatAx(m), 'square'); % shouldn't do much if i set the dims properly
%     xlabel('');
%     ylabel('');
%     colorbar off;
%     
%     if max(rfmapSmoothThisOri(:)) > maxRfmapSmooth
%         maxRfmapSmooth = max(rfmapSmoothThisOri(:));
%     end
% end
% 
% % set colorbar axis to be the same for all individual ori plots
% for m = 1:numel(gratingAngles)
%     caxis(miniHeatAx(m), [0 maxRfmapSmooth]);
% end
% 
plotFileName = sprintf('%s/%s-%s.png', processedDataDir, channelName, blockName);
export_fig(plotFileName, '-nocrop');
close;

end


