%% setup
clear;

addpath(genpath('matlab_sdk'));
addpath(genpath('C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\util\'));
addpath(genpath('C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\'));
addpath(genpath('C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\bursting'));
% the data file has online sorted spikes -- ignore those; the offline
% sorted spikes in the spikes file are much better

%%
dataFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\20170308\20170308-mcc-pul-d1_34.0mm_d2_34.45mm_rest2.pl2';
spikesFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\20170308\20170308-mcc-pul-d1_34.0mm_d2_34.45mm_rest2-sort2.pl2';

animalName = 'M';
sessionName = '20170308';
areaName = 'PUL';

%%
dataFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170127\20170127-all_merged.pl2';
spikesFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170127\20170127-all_merged-sort1.pl2';

animalName = 'M';
sessionName = '20170127';
areaName = 'PUL';

%%
dataFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170327\20170327-mcc-pul-d1_34.00mm_d2_34.00mm_allExceptRest_merged_noWB-sort1.pl2';
spikesFileName = dataFileName;

animalName = 'M';
sessionName = '20170327';
areaName = 'PUL';

%% load pl2 data
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
D = loadPL2(spikesFileName, sessionName, animalName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc); 

%% plot all spike waveforms
nWfTime = numel(D.allSpikeStructs{1}.meanWf);
spikeFs = 40000;
t = (1:nWfTime)/(spikeFs/1000);
textParams = {'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'none'};
yBounds = [-0.1 0.08];
fsClassThresh = 0.35/1000; 
for i = 1:numel(D.allSpikeStructs)
    figure;
    hold on;
    plot(t, D.allSpikeStructs{i}.meanWf, 'LineWidth', 3);
    plot([t(1) t(end)], [0 0], 'Color', 0.5*ones(3, 1));
    sd = D.allSpikeStructs{i}.sdWf;
    jbfill(t, D.allSpikeStructs{i}.meanWf + sd, D.allSpikeStructs{i}.meanWf - sd, lines(1), ones(3, 1), 1, 0.5);
    xlabel('Time (ms)');
    ylim(yBounds);
    title(sprintf('Unit %d Waveform (Mean +/- SD)', i));
    if startsWith(D.allSpikeStructs{i}.inflectionPattern, 'tp') || ...
            (startsWith(D.allSpikeStructs{i}.inflectionPattern, 'ptp') && D.allSpikeStructs{i}.firstPeakSmoothedAmp > abs(D.allSpikeStructs{i}.firstTroughSmoothedAmp))
        if D.allSpikeStructs{i}.troughToPeakTime <= fsClassThresh
            fsClass = 'Narrow';
        else
            fsClass = 'Broad';
        end
    else
        fsClass = '';
    end
    
    text(0.63, 0.27, {D.allSpikeStructs{i}.name, ...
            sprintf('Session Name: %s', D.allSpikeStructs{i}.sessionName), ...
            sprintf('Channel ID: %d', D.allSpikeStructs{i}.channelID), ...
            sprintf('Num Spikes: %d', numel(D.allSpikeStructs{i}.ts)), ...
            sprintf('Mean SD: %0.3f', mean(std(D.allSpikeStructs{i}.wf))), ...
            sprintf('Trough to Peak Time: %0.2f ms', D.allSpikeStructs{i}.troughToPeakTime*1000), ...
            sprintf('Num Inflections: %d', D.allSpikeStructs{i}.numInflections), ...
            sprintf('Inflection Pattern: %s', D.allSpikeStructs{i}.inflectionPattern), ...
            sprintf('Putative classification: %s', fsClass), ...
            sprintf('SPKC High-pass Filter: 300 Hz')}, ...
            textParams{:});
    
    plotFileName = sprintf('%s-waveform-unit%d.png', ...
            sessionName, i);
    export_fig(plotFileName, '-nocrop');
    close;
end

%% for trough-first waveforms, plot distribution of trough to peak time
troughToPeakTimeAll = nan(numel(D.allSpikeStructs), 1);
for i = 1:numel(D.allSpikeStructs)
    spikeStruct = D.allSpikeStructs{i};
    if all(spikeStruct.peakSmoothedAmps > 0) && ...
            (startsWith(spikeStruct.inflectionPattern, 'tp') || ...
            (startsWith(spikeStruct.inflectionPattern, 'ptp') && ...
            spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
            spikeStruct.peakSmoothedAmps(1) < 1/4 * max(spikeStruct.peakSmoothedAmps)))
        troughToPeakTimeAll(i) = D.allSpikeStructs{i}.troughToPeakTime;
    end
end
troughToPeakTimeAll(isnan(troughToPeakTimeAll)) = [];
troughToPeakTimeAll = troughToPeakTimeAll * 1000; % convert s to ms
figure;
hold on;
histogram(troughToPeakTimeAll, 0:0.05:0.8);
origYLim = ylim();
plot([0.35 0.35], [0 origYLim(2)], 'k-', 'LineWidth', 3);
ylim(origYLim);
xlabel('Trough to Peak Time (ms)');
ylabel('Number of Units');
title(sprintf('Trough to Peak Times for Trough-First Waveforms (N=%d/%d)', ...
        numel(troughToPeakTimeAll), numel(D.allSpikeStructs)));

plotFileName = sprintf('%s-troughToPeakTimes.png', sessionName);
export_fig(plotFileName, '-nocrop');
close;

figure;
hold on;
grid on;
fsClassThresh = 0.35/1000;
uncClassThresh = 0.35/1000;
cols = lines(2);
nWfTime = numel(D.allSpikeStructs{i}.meanWf);
spikeFs = 40000;
t = (1:nWfTime)/(spikeFs/1000);
for i = 1:numel(D.allSpikeStructs)
    spikeStruct = D.allSpikeStructs{i};
    if all(spikeStruct.peakSmoothedAmps > 0) && ...
            (startsWith(spikeStruct.inflectionPattern, 'tp') || ...
            (startsWith(spikeStruct.inflectionPattern, 'ptp') && ...
            spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
            spikeStruct.peakSmoothedAmps(1) < 1/4 * max(spikeStruct.peakSmoothedAmps)))
        tShift = spikeStruct.firstTroughIndex/(spikeFs/1000);
        if spikeStruct.troughToPeakTime <= fsClassThresh
             plot(t - tShift, spikeStruct.meanWf / max(spikeStruct.meanWf), 'Color', cols(1,:));
             fprintf('%s - NS\n', spikeStruct.name);
        elseif spikeStruct.troughToPeakTime <= uncClassThresh
             plot(t - tShift, spikeStruct.meanWf / max(spikeStruct.meanWf), 'Color', 0.3*ones(3, 1));
             fprintf('%s - UNC\n', spikeStruct.name);
        else
             plot(t - tShift, spikeStruct.meanWf / max(spikeStruct.meanWf), 'Color', cols(2,:));
             fprintf('%s - BS\n', spikeStruct.name);
        end
    end
end
xlim([-0.2 0.6]);
ylim([-3 1]);
title(sprintf('Mean Waveforms (N=%d/%d)', ...
        numel(troughToPeakTimeAll), numel(D.allSpikeStructs)));
xlabel('Time from Trough (ms)');
ylabel('Normalized Amplitude');
plotFileName = sprintf('%s-troughFirstWfsOverlaid.png', sessionName);
export_fig(plotFileName, '-nocrop');
close;

[dip,p_value] = HartigansDipSignifTest(sort(troughToPeakTimeAll), 1000)

% TODO
% TODO
% compute trough to peak using first peak not max peak

% save('temp-grant-troughToPeakTimeAll-temp1.mat', 'troughToPeakTimeAll');
save('temp-grant-troughToPeakTimeAll-temp2.mat', 'troughToPeakTimeAll');

%%
a1 = load('temp-grant-troughToPeakTimeAll-temp1.mat');
a2 = load('temp-grant-troughToPeakTimeAll-temp2.mat');

%%
troughToPeakTimeAll = [a1.troughToPeakTimeAll; a2.troughToPeakTimeAll];

%% temp for grant figure (20170127 and 20170327)

figure_tr_inch(2.25, 2);
set(gcf, 'Color', 'white');
subaxis(1, 1, 1, 'ML', 0.2, 'MT', 0.05, 'MR', 0.05, 'MB', 0.25);
hold on;
histogram(troughToPeakTimeAll, 0:0.05:0.8);
xlim([0.1 0.75]);
box off;
origYLim = [0 30];%ylim();
plot([0.35 0.35], [0 origYLim(2)], 'k-', 'LineWidth', 3);
ylim(origYLim);
xlabel('Trough to Peak Time (ms)');
ylabel('Number of Units');
% title(sprintf('Spike Waveform Trough to Peak Times (N=%d Units)', ...
%         numel(troughToPeakTimeAll)));
% set(gca, 'FontSize', 12);
    
plotFileName = sprintf('grantFig-troughToPeakTimes.png');
export_fig(plotFileName, '-nocrop');

%% plot all waveforms overlaid
figure;
hold on;
grid on;
nWfTime = numel(D.allSpikeStructs{i}.meanWf);
spikeFs = 40000;
t = (1:nWfTime)/(spikeFs/1000);
for i = 1:numel(D.allSpikeStructs)
    plot(t, D.allSpikeStructs{i}.meanWf);
end
title(sprintf('Mean Waveforms (N=%d)', numel(D.allSpikeStructs)));
xlabel('Time (ms)');
plotFileName = sprintf('%s-allWfsOverlaid.png', sessionName);
export_fig(plotFileName, '-nocrop');
close;

%% histograms for spike rate for BS vs NS cell classes

% count up spikes in the time between 25th percentile of spike times and
% 75th percentile, then divide by that duration to get spike rate
% do this to handle spikes that appear or fade over time
fsClassThresh = 0.35/1000; 

frNSCells = nan(numel(D.allSpikeStructs), 1);
frBSCells = nan(numel(D.allSpikeStructs), 1);
durations = nan(numel(D.allSpikeStructs), 1);
for i = 1:numel(D.allSpikeStructs)
    if (startsWith(spikeStruct.inflectionPattern, 'tp') || ...
            (startsWith(spikeStruct.inflectionPattern, 'ptp') && ...
            spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
            spikeStruct.peakSmoothedAmps(1) < spikeStruct.peakSmoothedAmps(2))) && ...
            all(spikeStruct.peakSmoothedAmps > 0)
        ts = D.allSpikeStructs{i}.ts;
        spikeTime25Pctile = quantile(ts, 0.25);
        spikeTime75Pctile = quantile(ts, 0.75);
        duration = spikeTime75Pctile - spikeTime25Pctile;
        numSpikes = numel(ts) / 2; % 50% of spikes
        if duration >= 5*60 % cell should be stable at least 5 minutes
            if D.allSpikeStructs{i}.troughToPeakTime <= fsClassThresh
                % narrow spiking
                frNSCells(i) = numSpikes / duration;
            else
                % broad spiking
                frBSCells(i) = numSpikes / duration;
            end
        end
    end
end
frNSCells(isnan(frNSCells)) = [];
frBSCells(isnan(frBSCells)) = [];
frNSCells(frNSCells < 0.05) = []; % min firing rate to be considered
frBSCells(frBSCells < 0.05) = [];

median(frNSCells) - median(frBSCells)

figure;
histogram(frNSCells, 0:2:16);

figure;
histogram(frBSCells, 0:2:16);

[H,P,CI,STATS] = ttest2(frNSCells, frBSCells)

fprintf('%d/%d = %0.1f%%\n', numel(frNSCells), (numel(frNSCells) + numel(frBSCells)), ...
        numel(frNSCells) / (numel(frNSCells) + numel(frBSCells)) * 100);


%% generate spike time series
% why do spike times only have 4 significant figures?? -- suggests 10KHz
% sample rate. Doesn't really matter. Bin spikes anyway.
origSpikeFs = D.allSpikeStructs{1}.Fs;
spikeFs = 1000; % 1000 Hz, period 1 ms
spikeTimeSeriesAll = false(numel(D.allSpikeStructs), round(D.totalTicks/(origSpikeFs/spikeFs)));
for i = 1:numel(D.allSpikeStructs)
    spikeTimeSeriesAll(i,(ceil(spikeFs*D.allSpikeStructs{i}.ts))) = true;
end
% TODO the above is transposed from how spikeTimeSeries are used below.
% check this before converting the below

%% correlate each spike channel
corrMat = corr(spikeTimeSeriesAll');
corrMat(eye(size(corrMat))~=0) = NaN; % blank out diagonal

figure_tr_inch(10, 10);
imagesc(corrMat);
xlabel('Unit Number');
ylabel('Unit Number');
title(sprintf('Spike Correlation Matrix (All Time Periods, Bin %d ms)', round(1000/spikeFs)));
colorbar;
origCAxis = caxis;
caxis([origCAxis(1) min(origCAxis(2), 0.15)]);
set(gca, 'FontSize', 16);

plotFileName = sprintf('%s-spike_corr-bin%dms.png', sessionName, round(1000/spikeFs));
export_fig(plotFileName, '-nocrop');
% close;

%% correlate each spike channel (cells with inflections tp* only)
isUnitCell = cellfun(@(x) startsWith(x.inflectionPattern, 'tp'), D.allSpikeStructs);
spikeTimeSeriesCells = spikeTimeSeriesAll;
spikeTimeSeriesCells(~isUnitCell,:) = NaN;
% spikeTimeSeriesCells = spikeTimeSeriesAll(isUnitCell,:);
corrMat = corr(spikeTimeSeriesCells');
corrMat(eye(size(corrMat))~=0) = NaN; % blank out diagonal

figure_tr_inch(10, 10);
imagesc(corrMat);
xlabel('Unit Number');
ylabel('Unit Number');
title(sprintf('Spike Correlation Matrix (All Time Periods, Bin %d ms)', round(1000/spikeFs)));
colorbar;
origCAxis = caxis;
caxis([origCAxis(1) min(origCAxis(2), 0.15)]);
set(gca, 'FontSize', 16);

plotFileName = sprintf('%s-spike_corr_cellsOnly-bin%dms.png', sessionName, round(1000/spikeFs));
export_fig(plotFileName, '-nocrop');
% close;

%% autocorr unit 1
maxLag = 200; % 200 ms in either direction
x = (-maxLag:maxLag)/spikeFs;
nMovMeanPts = 10;
figure;
hold on;
for i = 1
    fprintf('xcorr(spikes(ch%d))...\n', i);
    xc = xcorr(spikeTimeSeriesAll(i,:), maxLag);
    xc(maxLag + 1) = NaN; % this is just the number of spikes
    plot(x, xc); 
    plot(x, movmean(xc, nMovMeanPts), 'LineWidth', 3);
    xlim([x(1) x(end)]);
    xlabel('Lag (s)');
    grid on;
end

%% xcorr unit x, y
maxLag = 200; % 200 ms in either direction
nXcPts = 2*maxLag + 1;
xcAll = nan(numel(D.allSpikeStructs), numel(D.allSpikeStructs), nXcPts);

for i = 1:numel(D.allSpikeStructs)
    for j = i:numel(D.allSpikeStructs)
        fprintf('xcorr(spikes(ch%d), spikes(ch%d)...\n', i, j);
        xcAll(i,j,:) = xcorr(spikeTimeSeriesAll(i,:), spikeTimeSeriesAll(j,:), maxLag);
        if i == j
            xcAll(i,j,maxLag+1) = NaN; % this is just the number of spikes
        end
    end
end

% fill in nans of xcAll - the other triangle
for i = 1:numel(D.allSpikeStructs)
    for j = 1:(i-1)
        fprintf('xcorr(spikes(ch%d), spikes(ch%d)...\n', i, j);
        xcAll(i,j,:) = xcAll(j,i,:);
    end
end

saveFileName = sprintf('%s-xcAll_maxLag%d_spikeFs%d.mat', sessionName, maxLag, spikeFs);
savefast(saveFileName, 'xcAll', 'maxLag', 'spikeFs');

%% plot xcorr unit a, b
nMovMeanPts = 10;
maxLag = 200; % 200 ms in either direction
% saveFileName = sprintf('%s-xcAll_maxLag%d_spikeFs%d.mat', sessionName, maxLag, spikeFs);
saveFileName = sprintf('%s-xcAll.mat', sessionName);
load(saveFileName, 'xcAll', 'maxLag', 'spikeFs');

unitInd1 = 112;
unitInd2 = 11;
figure;
hold on;
plot((-maxLag:maxLag)/spikeFs, movmean(squeeze(xcAll(unitInd1,unitInd2,:)), nMovMeanPts), 'LineWidth', 3);
xlim([-maxLag maxLag]/spikeFs);
xlabel('Lag (s)');
title(sprintf('Cross-correlation between Unit %d and Unit %d', unitInd1, unitInd2));

%% plot xcorr unit 1, x
nMovMeanPts = 10;

unitInd1 = 29;
figure;
hold on;
for j = 1:numel(D.allSpikeStructs)
    y = squeeze(xcAll(unitInd1,j,:));
    plot((-maxLag:maxLag)/spikeFs, movmean(y, nMovMeanPts), 'LineWidth', 3);
end
origYLim = ylim();
plot([0 0], origYLim, '-', 'Color', 0.5*ones(3, 1));
xlim([-maxLag maxLag]/spikeFs);
xlabel('Lag (s)');

%% plot xcorr unit 1, x
nMovMeanPts = 10;

unitInd1 = 51;
figure_tr_inch(5, 10);
hold on;
yCount = 1;
cols = lines(6);
for j = 1:numel(D.allSpikeStructs)
    y = squeeze(xcAll(unitInd1,j,:));
    yEnds = y([1:round(1/4*numel(y)) round(3/4*numel(y)):end]);
    yNorm = (y - nanmean(yEnds))/std(yEnds);
    if startsWith(D.allSpikeStructs{j}.inflectionPattern, 'tp')
        lineStyle = '-';
        if j == unitInd1
            col = 'k';
        else
            col = cols(mod(yCount, 6) + 1, :);
        end
        plot((-maxLag:maxLag)/spikeFs, movmean(yNorm, nMovMeanPts)*2 - yCount*4, ...
            'LineWidth', 2, 'LineStyle', lineStyle, 'Color', col);
        text(-0.6, -yCount*4, sprintf('%d', D.allSpikeStructs{j}.channelID));
        yCount = yCount + 1;
    else
        lineStyle = '-.';
    end
%     plot((-maxLag:maxLag)/spikeFs, movmean(yNorm, nMovMeanPts)*2 - yCount*4, ...
%             'LineWidth', 2, 'LineStyle', lineStyle);
%     yCount = yCount + 1;
end
origYLim = ylim();
plot([0 0], origYLim, '-', 'Color', 0.5*ones(3, 1));
xlim([-maxLag maxLag]/spikeFs);
xlabel('Lag (s)');
xlim([-0.1 0.1]);
set(gca, 'YTickLabel', []);
set(gca, 'YTick', []);
title(sprintf('Cross-correlation of Units with Unit %d', unitInd1));




%% image plot of all xcorr normalized
% could subtract the mean or subtract the tail
xcNormAll = xcAll - nanmean(xcAll, 3);

% could zero out the autocorr channel - a bursty cell will heavily skew the
% color axis

for i = 29;%1:numel(D.allSpikeStructs)
    fprintf('xcorr(spikes(ch%d), spikes(all)...\n', i);
    figure_tr_inch(10, 10);
    imagesc((-maxLag:maxLag)/spikeFs, 1:numel(D.allSpikeStructs), squeeze(xcNormAll(i,:,:)));
    xlim([-maxLag maxLag]/spikeFs);
    xlabel('Lag (s)');
    ylabel('Unit Number');
    title(sprintf('Normalized Spike Xcorr with Unit %d', i));
    caxis([-1000 1000])
    colorbar;
    set(gca, 'FontSize', 16);
    
%     plotFileName = sprintf('%s-spike_xcorr_norm_unit%d.png', ...
%             sessionName, i);
%     export_fig(plotFileName, '-nocrop');
%     close;
end

%% look for burstiness in cell x
minSilentPeriod = 0.02; % 0.1;
minISIWithinBurst = 0.004;
Fs = 1000;

for i = 1:numel(D.allSpikeStructs)
    spikeTimes = D.allSpikeStructs{i}.ts;
    [diffSpikeTimesAll, silentPeriodDurationsPrecedingSpikeAll, numSpikesTotalAll, spikesInBurstAll, averageISIWithinBurstAll] = ...
                burstAnalysis2(spikeTimes, min(spikeTimes), max(spikeTimes), minSilentPeriod, minISIWithinBurst);

    diffSpikeTimes = diffSpikeTimesAll{1}; % only one window
    fprintf('Note that there are %d spikes with ISI < 1ms.\n', sum(diffSpikeTimes < 0.001));
    
    otherUnitsThisChannel = findOtherUnitsSameCh(D.allSpikeStructs, i);
    
    titleName = sprintf('Burstiness - Unit %d', i);
    plotFileName = sprintf('%s-burstiness-unit%d.png', ...
                sessionName, i);

    plotBurstAnalysis(silentPeriodDurationsPrecedingSpikeAll, numel(spikeTimes), ...
            spikesInBurstAll, diffSpikeTimes, Fs, minSilentPeriod, ...
            minISIWithinBurst, i, D.allSpikeStructs, titleName, plotFileName);
    close;
end

%% create N "trials"
nTrials = 500;
minTrials = 5;
origSpikeFs = 40000;
spikeFs = 1000;
trialDuration = 1000; % samples in ms
% total duration examined = 500 seconds ~= 8 minutes
for i = 1:numel(D.allSpikeStructs)
    fprintf('creating trials for unit %d...\n', i);
    spikeTimeSeriesTrials = zeros(nTrials, trialDuration);
    spikeTimesInSpikeFs = round(spikeFs*D.allSpikeStructs{i}.ts);
    spikeTimesInWindow = cell(nTrials, 1);
    for j = 1:nTrials
        % "trials" are consecutive time segments starting with time point 1
        % which is not necessarily the best... but to start off
        binaryTimeSpikeStartTime = (j-1)*trialDuration + 1; 
        spikeTimesInWindow{j} = spikeTimesInSpikeFs(spikeTimesInSpikeFs >= binaryTimeSpikeStartTime & ...
                spikeTimesInSpikeFs < (binaryTimeSpikeStartTime + trialDuration)) ...
                - binaryTimeSpikeStartTime + 1;

        % create fake data
%         nBurstsPerTrial = 2;
%         nSpikesInBurst = 3;
%         randTimes = randi(size(spikeTimeSeriesTrials, 2) - (nSpikesInBurst - 1), nBurstsPerTrial);
%         if j < nTrials/2
%             if nSpikesInBurst == 2
%                 spikeTimesInWindow{j} = [randTimes; randTimes+1];
%             elseif nSpikesInBurst == 3
%                 spikeTimesInWindow{j} = [randTimes; randTimes+1; randTimes+2];
%             end
%         else
%             spikeTimesInWindow{j} = randTimes;
%         end

        spikeTimeSeriesTrials(j,spikeTimesInWindow{j}) = 1;
    end
    
    % remove all "trials" with no spikes
    spikeTimeSeriesTrials(all(spikeTimeSeriesTrials == 0, 2), :) = [];
    if size(spikeTimeSeriesTrials, 1) < minTrials
        fprintf('not enough spikes in trials to compute SPNA for unit %d. skipping...\n', i);
        continue;
    end
    
    fprintf('computing SPNA for unit %d...\n', i);
    % compute auto-correlation with up to 10ms lag in either direction
    % TODO: use coeff scaling so that autocorr at zero lag = 1? this will
    % typically reduce all the values
    maxLag = spikeFs / 100; % 10ms
    lagT = (1:maxLag) * 1000 / spikeFs;
    
    % SPNA is a poor, underpowered method for detecting bursts here:
    % 1) does not pick up on sparse bursting on say half of trials using 
    % fake data
    % 2) if 1 trial has a prolonged burst and no other trials do, the cell
    % is identified as bursting
    % 3) cells with obvious bursting patterns as revealed by the ISI
    % distribution are not identified as bursting
    % 4) it is computationally slow and requires splitting the data into
    % trials
    
    [SPNA, meanAutoCorr, meanCrossCorr, stdCrossCorr] = computeSPNA(spikeTimeSeriesTrials, maxLag);
    
    % BRI: burstiness-refractoriness index:
    % average of SPNA over lag of 1-4 ms
    BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
    BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
    BRI = mean(SPNA(BRILagStartInd:BRILagEndInd));
    
    titleName = sprintf('SPNA Burstiness - Unit %d', i);
    plotFileName = sprintf('%s-spna-burstiness-unit%d.png', ...
                sessionName, i);

    plotSPNABurstAnalysis(spikeTimesInWindow, spikeFs, lagT, maxLag, ...
            meanAutoCorr, meanCrossCorr, stdCrossCorr, SPNA, BRI, ...
            titleName, plotFileName);
    close;
end


%% TODO shuffle test

%% plot unit 1 and unit 138 - eh very hard to see by eye
unitInd1 = 27;
unitInd2 = 43;

maxLag = 200; % 200 ms in either direction
% figure;
% plot((-maxLag:maxLag)/spikeFs, movmean(xcAll(138,:), nMovMeanPts), 'LineWidth', 3);

spikeTimeSeries1 = zeros(round(D.totalTicks/(origSpikeFs/spikeFs)), 1);
spikeTimeSeries1(round(spikeFs*D.allSpikeStructs{unitInd1}.ts)) = 1; % look at unit 1
spikeTimeSeries2 = zeros(round(D.totalTicks/(origSpikeFs/spikeFs)), 1);
spikeTimeSeries2(round(spikeFs*D.allSpikeStructs{unitInd2}.ts)) = 1; % look at unit 138
tInd = (5000:6000)*spikeFs; % 10 seconds
figure;
hold on;
plot(tInd, spikeTimeSeries1(tInd));
plot(tInd, spikeTimeSeries2(tInd)+1);

D.allSpikeStructs{unitInd1}
D.allSpikeStructs{unitInd2}

%% plot isi distribution
bins = 0:0.002:0.2;
for i = 1:numel(D.allSpikeStructs)
    figure;
    hold on;
    grid on;
    histogram(diff(D.allSpikeStructs{i}.ts), bins);
    xlabel('Interspike Interval (s) (2ms bins)');
    ylabel('Frequency');
    title(sprintf('ISI Distribution for Unit %d', i));
    
    plotFileName = sprintf('%s-spike_isi_dist_unit%d.png', ...
            sessionName, i);
    export_fig(plotFileName, '-nocrop');
    close;
end


%%
lfpFs = 1000;
isLoadLfp = 1;
D = loadPL2(spikesFileName, sessionName, animalName, areaName, isLoadLfp); 

%% process FP channels

% this is a big variable - probably better to load from pl2 file one at a time as needed!!
% 1.8 GB for an hour of recording and 64 channels
allLfp = D.lfps;
biLfp = diff(allLfp); % local channel differences


%% xcorr between ch1 and all other lfps
maxLag = lfpFs; % one second in either direction
figure;
hold on;
for i = 1:size(allLfp, 1)
    fprintf('xcorr(lfp(ch1), lfp(ch%d))...\n', i);
    xc = xcorr(allLfp(1,:), allLfp(i,:), maxLag); 
    if i == 1
        lineWidth = 3;
    else
        lineWidth = 1;
    end
    plot((-maxLag:maxLag)/lfpFs, xc, 'LineWidth', lineWidth); 
    xlim([-maxLag maxLag]/lfpFs);
end

%% xcorr between bipolar lfp ch1 and all other bipolar lfp
maxLag = lfpFs; % one second in either direction
figure;
hold on;
for i = 1:size(biLfp, 1)
    fprintf('xcorr(biLfp(ch1), biLfp(ch%d))...\n', i);
    xc = xcorr(biLfp(1,:), biLfp(i,:), maxLag); 
    if i == 1
        lineWidth = 3;
    else
        lineWidth = 1;
    end
    plot((-maxLag:maxLag)/lfpFs, xc, 'LineWidth', lineWidth); 
    xlim([-maxLag maxLag]/lfpFs);
    xlabel('Lag (s)');
end

%% autocorr each LFP channel
maxLag = lfpFs; % one second in either direction
figure;
hold on;
for i = 1:size(allLfp, 1)
    fprintf('xcorr(lfp(ch%d))...\n', i);
    xc = xcorr(allLfp(i,:), maxLag);
    plot((-maxLag:maxLag)/lfpFs, xc); 
    xlim([-maxLag maxLag]/lfpFs);
    xlabel('Lag (s)');
end

%% correlate each LFP channel
corrMat = corr(allLfp');
corrMat(eye(size(corrMat))~=0) = NaN; % blank out diagonal
figure_tr_inch(10, 10);
imagesc(corrMat);
xlabel('Channel Number');
ylabel('Channel Number');
title('LFP Correlation Matrix (All Time Periods)');
colorbar;

plotFileName = sprintf('%s-lfp_corr.png', sessionName);
export_fig(plotFileName);
close;





%% plot spkcs -- for grant figure 5 (use 20170329)
chsToPlot = 1:6:32;
yScale = 8;
yOffset = 1;
tInd = (1:20000) + 40000*6;

figure_tr_inch(2, 3);
subaxis(1, 1, 1, 'MT', 0.02, 'MB', 0.02, 'ML', 0.26, 'MR', 0.03);
hold on;
box off;
for i = 1:numel(chsToPlot)
    plot(D.spkcs(chsToPlot(i),tInd)*yScale + i*yOffset, 'LineWidth', 2);
end
set(gca, 'YDir', 'reverse');
ylim([1-yOffset+0.5 numel(chsToPlot) + yOffset-0.5]);
ylabel('Channel Number');
set(gca, 'XTick', []);
set(gca, 'XColor', 'none');
set(gca, 'YTick', 1:numel(chsToPlot));
set(gca, 'YTickLabel', chsToPlot);
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');

plotFileName = sprintf('%s-spkc_series.png', sessionName);
export_fig(plotFileName, '-nocrop');


%% run cell assembly detection script
addpath(genpath('C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\Cell_assembly_detection'));

% find max number of spikes
maxNSpikes = 0;
for i = 1:numel(D.allSpikeStructs)
    if numel(D.allSpikeStructs{i}.ts) > maxNSpikes
        maxNSpikes = numel(D.allSpikeStructs{i}.ts);
    end
end

% create data matrix
nCellsToProcess = numel(D.allSpikeStructs); % 20 cells, all ts = 14 minutes
allSpikesMat = nan(nCellsToProcess, maxNSpikes);
for i = 1:nCellsToProcess
    allSpikesMat(i,1:numel(D.allSpikeStructs{i}.ts)) = D.allSpikeStructs{i}.ts;
end

% assembly search settings
binSizes = [0.015 0.025 0.04 0.06 0.085 0.15 0.25];
maxLags = 10*ones(size(binSizes));

tic;
[assembly] = Main_assemblies_detection(allSpikesMat, maxLags, binSizes);
toc;
% TODO save assembly to file

%%
[As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly, binSizes);

saveFileName = sprintf('%s-assembly.mat', sessionName);
savefast(saveFileName, 'allSpikesMat', 'maxLags', 'binSizes', 'assembly', ...
        'As_across_bins', 'As_across_bins_index');

figure_tr_inch(10, 10);
display = 'ordunit';
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins, nCellsToProcess, binSizes, display);
% note colorbar was removed from the created figure temporarily because the
% bars below did not align with them. TODO fix this

plotFileName = sprintf('%s-assembly_assignment.png', sessionName);
export_fig(plotFileName);
close;

%%
criteria = 'biggest';
[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins, As_across_bins_index, nCellsToProcess,criteria);

figure_tr_inch(10, 10);
display = 'ordunit';
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins_pr, nCellsToProcess, binSizes, display);

plotFileName = sprintf('%s-assembly_assignment_pruned.png', sessionName);
export_fig(plotFileName);
close;


%% 
% plot assembly time courses