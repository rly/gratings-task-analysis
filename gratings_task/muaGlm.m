% train GLM on N neurons' activity, test on N+1th neuron
% is performance better than GLM trained on shuffled n-1 data

% This GLM performs well, but so does a GLM trained on spike-jittered data
% with a jitter of +/- 400 ms. It seems that having 31 units in the GLM,
% whether activity is structured or not, may be sufficient to capture most
% of the variance in spiking activity of many individual units, though it
% fails to capture some sharp visual-evoked transients (11 units had p<0.05
% compared to null distribution.

clear;
sessionName = 'M20170311';
sessionInd = 7;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
channelsToLoad = 1:32;
muaChannelsToLoad = 1:32;
lfpChannelsToLoad = 1:32;
lfpChannels = lfpChannelsToLoad;
numRandomizations = 2;
isLoadMua = 1;
isLoadSortedSua = 0;

v = 12;

outputDir = sprintf('%s/%s/', processedDataRootDir, 'MUA_GLM');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%%
taskName = 'GRATINGS';
scriptName = 'MUA_GRATINGS';
isZeroDistractors = 0;
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, taskName, scriptName, isLoadSortedSua, isLoadMua, 0, 0, 0);
sessionName = R.sessionName;
areaName = R.areaName;

%%   
totalTimeOverall = sum(D.blockStopTimes(R.blockIndices) - D.blockStartTimes(R.blockIndices));
minFiringRateOverall = 0.2;

nUnits = numel(D.allUnitStructs);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing Session %s %s, Channels %d-%d\n', sessionName, areaName, channelsToLoad([1 end]));
fprintf('Processing %d units...\n', nUnits);

% ES.cueOnset is struct for data -0.8 to 0.8 ms from cue onset
% ES.cueOnset.spdf uses spikes -0.7 to 0.7 ms from cue onset
binWidth = 0.01;
binBounds = -0.8:binWidth:0.8;
spikeMat = zeros(nUnits, numel(binBounds));
spdfMat = zeros(nUnits, 1401);

unitNames = cell(nUnits, 1);
for j = 1:nUnits
    unitStruct = D.allUnitStructs{j};
    unitName = unitStruct.name;
    spikeTimes = unitStruct.ts;
    unitNames{j} = unitName;

    saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            processedDataDir, unitName, blockName, v);
    fprintf('Processing %s...\n', saveFileName);

    ES = load(saveFileName);

    for k = 1:numel(ES.cueOnset.spikeTimes)
        spikeTs = ES.cueOnset.spikeTimes(k).times;
        for l = 1:numel(spikeTs)
            spikeMat(j,ceil(spikeTs(l)/binWidth)) = spikeMat(j,ceil(spikeTs(l)/binWidth)) + 1;
        end
    end
    spdfMat(j,:) = ES.cueOnset.spdf;
end
spdfT = ES.cueOnset.t + ES.cueOnset.spdfWindowOffset(1);

%%
% tInds = binBounds > -0.2 & binBounds <= 0.6;
% spdfTInds = spdfT > -0.2 & spdfT <= 0.6;
% for j = 1:nUnits
%     % leave out j-th unit
%     inds = 1:nUnits;
%     inds(inds == j) = [];
%     
%     % predict binned spiking activity of left out multi-unit
%     b = glmfit(spikeMat(inds,:)', spikeMat(j,:)', 'normal');
%     yfit = glmval(b, spikeMat(inds,:)', 'identity');
%     
%     figure_tr_inch(10, 6);
%     subaxis(1, 2, 1);
%     hold on;
%     plot(binBounds(tInds), spikeMat(j,tInds), '-');
%     plot(binBounds(tInds), yfit(tInds), '-');
%     
%     % predict spike density function of left out multi-unit
%     [b,dev] = glmfit(spdfMat(inds,:)', spdfMat(j,:)', 'normal');
%     yfit = glmval(b, spdfMat(inds,:)', 'identity');
%     
%     subaxis(1, 2, 2);
%     hold on;
%     plot(spdfT(spdfTInds), spdfMat(j,spdfTInds), '-');
%     plot(spdfT(spdfTInds), yfit(spdfTInds), '-');
% end

%% shuffle test
nShuffle = 200;
jitterTime = 0.4; % plus or minus 400 ms
spdfMatShuffle = zeros(nShuffle, nUnits, 1401);
for j = 1:nUnits
    unitStruct = D.allUnitStructs{j};
    unitName = unitStruct.name;
    spikeTimes = unitStruct.ts;

    saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            processedDataDir, unitName, blockName, v);
    fprintf('Processing %s...\n', saveFileName);

    ES = load(saveFileName);

    for m = 1:nShuffle
        % jitter each spike by random uniform [-400,400] ms within bounds
        spikeTimesJitter = struct();
        for k = 1:numel(ES.cueOnset.spikeTimes) % for each trial
            spikeTs = ES.cueOnset.spikeTimes(k).times;
            spikeTimesJitter(k).times = addBoundedRandom(spikeTs, jitterTime, [0 sum(ES.cueOnset.window)]);
        end
        spdfMatShuffle(m,j,:) = fixedPsth(spikeTimesJitter, ES.kernelSigma, 0, ES.cueOnset.t);
    end
end
spdfT = ES.cueOnset.t + ES.cueOnset.spdfWindowOffset(1);

%% shuffle test comparison
spdfTInds = spdfT > -0.1 & spdfT <= 0.5;
saveFileNames = cell(nUnits, 1);
pVals = nan(nUnits, 1);
for j = 1:nUnits
    % leave out j-th unit
    inds = 1:nUnits;
    inds(inds == j) = [];
    
    % predict spike density function of left out multi-unit
    [b,dev] = glmfit(spdfMat(inds,spdfTInds)', spdfMat(j,spdfTInds)', 'normal');
    yfit = glmval(b, spdfMat(inds,spdfTInds)', 'identity');
    
    devShuffle = nan(nShuffle, 1);
%     figure;
%     hold on;
%     plot(spdfT(spdfTInds), spdfMat(j,spdfTInds), '-', 'Color', cols(3,:), 'LineWidth', 2);
    for m = 1:nShuffle
%         inds = inds(1:10);
        [b,devShuffle(m)] = glmfit(squeeze(spdfMatShuffle(m,inds,spdfTInds))', spdfMat(j,spdfTInds)', 'normal');
        yfitShuffle = glmval(b, squeeze(spdfMatShuffle(m,inds,spdfTInds))', 'identity');
%         plot(spdfT(spdfTInds), yfitShuffle, '-', 'Color', cols(4,:), 'LineWidth', 2);
    end
    
    pVals(j) = min([sum(dev < devShuffle) / numel(devShuffle) sum(dev > devShuffle) / numel(devShuffle)]) * 2; % two-tailed
    
    cols = lines(4);
    
    figure_tr_inch(12, 6);
    subaxis(1, 2, 1, 'SH', 0.08, 'MB', 0.13);
    hold on;
    histogram(devShuffle);
    plot([dev dev], ylim(), 'r-', 'LineWidth', 2);
    xlabel('Deviance');
    ylabel('Number of shuffles');
    legend({'Null Distribution', sprintf('Actual Data p=%0.3f', pVals(j))});
    set(gca, 'FontSize', 14);
    
    subaxis(1, 2, 2);
    hold on;
    plot(spdfT(spdfTInds), spdfMat(j,spdfTInds), '-', 'Color', cols(3,:), 'LineWidth', 2);
    plot(spdfT(spdfTInds), yfit, '-', 'Color', cols(4,:), 'LineWidth', 2);
    xlabel('Time from Cue Onset (s)');
    ylabel('Estimated Spike Rate (Hz)');
    legend({'Actual', 'Predicted'});
    set(gca, 'FontSize', 14);
    
    suptitle(sprintf('GLM Prediction for Unit %s (N=%d)', unitNames{j}, nUnits-1), 'FontSize', 18);
    
    % save as pdf
    saveFileNames{j} = sprintf('%s/%s-withinSessionGlmPredict-v%d.pdf', outputDir, unitNames{j}, v);
    export_fig(saveFileNames{j}, '-nocrop');
end

%%
figure;
histogram(pVals, 0:0.025:1);
xlabel('P-Value of Deviance from Null Distribution');
ylabel('Number of Units');

%% collate all pdfs into one
appendSaveFileName = sprintf('%s/%s-all-withinSessionGlmPredict-v%d.pdf', outputDir, sessionName, v);
deleteAndAppendPdfs(appendSaveFileName, saveFileNames);

%% pca
spdfTInds = spdfT > -0.2 & spdfT <= 0.6;
data = spdfMat(:,spdfTInds);
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca(data);

fprintf('PCA: %d variables, %d observations\n', size(data, 2), size(data, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC2 explains %0.1f%% of the variance.\n', pcaPctExplained(2));
fprintf('\tPC3 explains %0.1f%% of the variance.\n', pcaPctExplained(3));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', sum(pcaPctExplained(1:2)));
fprintf('\tPC1 + PC2 + PC3 explain %0.1f%% of the variance.\n', sum(pcaPctExplained(1:3)));
fprintf('\tFirst %d PCs explain >99%% of the variance.\n', find(cumsum(pcaPctExplained) > 99, 1, 'first'));
fprintf('\tFirst %d PCs explain >99.9%% of the variance.\n', find(cumsum(pcaPctExplained) > 99.9, 1, 'first'));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
sh = scatter(pcaScore(:,1), pcaScore(:,2), 100, 'k');
sh.MarkerFaceAlpha = 0.8;
sh.LineWidth = 2;
xlabel('PC1');
ylabel('PC2');

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
plot(spdfT(spdfTInds), pcaCoeff(:,1), 'LineWidth', 9);
plot(spdfT(spdfTInds), pcaCoeff(:,2), 'LineWidth', 5);
plot(spdfT(spdfTInds), pcaCoeff(:,3), 'LineWidth', 1);

