function lfpAnalysisSummary(processedDataRootDir, recordingInfoFileName, sessionInds)

% clear;
% readDataLocally;
% sessionInds = 1:23;

v = 12;

%%
recordingInfo = readRecordingInfo(recordingInfoFileName);

outputDir = sprintf('%s/%s/', processedDataRootDir, 'LFP_GRATINGS_SUMMARY');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

nSessions = numel(sessionInds);
nLfpsApprox = nSessions * 16; % should be equal or an underestimate
lfpCount = 0;

isInDPulvinar = false(nLfpsApprox, 1);
isInVPulvinar = false(nLfpsApprox, 1);

fAxis = NaN;
baselinePower = nan(nLfpsApprox, 199);
cueResponsePower = nan(nLfpsApprox, 199);
cueTargetDelayPowerAllLocs = nan(nLfpsApprox, 199);
cueTargetDelayPowerP3 = nan(nLfpsApprox, 199);
cueTargetDelayPowerP1 = nan(nLfpsApprox, 199);

baselineWindowOffset = [-0.25 0];
cueResponseOffset = [0 0.25];
cueTargetDelayOffset = [-0.25 0];

% chronux parameters
params.tapers = [2 3];
params.fpass = [5 200];
params.pad = 2;
params.Fs = 1000;
params.trialave = 1;

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    lfpChannelsToLoad = R.lfpChannelsToLoad;
    blockName = strjoin(R.blockNames(R.gratingsTask3DIndices), '-'); % temp, only 3D
    processedDataDir = sprintf('%s/%s/LFP_GRATINGS/', processedDataRootDir, sessionName);
    fileNamePrefix = sprintf('%s-ch%d-ch%d-%s', sessionName, lfpChannelsToLoad([1 end]), blockName);
    saveFileName = sprintf('%s/%s-evokedLfps-v%d.mat', processedDataDir, fileNamePrefix, v);
    fprintf('Loading file %s...\n', saveFileName);
    EL = load(saveFileName);
        
    for j = 1:numel(EL.channelInds)
        fprintf('Processing channel %d...\n', EL.channelInds(j));
        lfpCount = lfpCount + 1;
        isInVPulvinar(lfpCount) = 0;
        isInDPulvinar(lfpCount) = 0;
        if ismember(EL.channelInds(j), R.vPulChannels)
            isInVPulvinar(lfpCount) = 1;
        end
        if ismember(EL.channelInds(j), R.dPulChannels)
            isInDPulvinar(lfpCount) = 1;
        end
        
        baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);
        cueResponseInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, cueResponseOffset);
        cueTargetDelayInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, cueTargetDelayOffset);

        
        % compute power in baseline - should probably be bipolar reference
        % all locations together
        cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(j,:,:))';
        arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(j,:,:))';
        arrayOnsetLfpP3Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 3,:))';
        arrayOnsetLfpP1Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 1,:))';
        
        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        [baselinePower(lfpCount,:),fAxis] = mtspectrumc(preCueBaselineLfps, params);
        
        cueResponseLfps = cueOnsetLfpCurrent(cueResponseInd,:);
        [cueResponsePower(lfpCount,:),fAxis] = mtspectrumc(cueResponseLfps, params);
        
        cueTargetDelayLfps = arrayOnsetLfpCurrent(cueTargetDelayInd,:);
        [cueTargetDelayPowerAllLocs(lfpCount,:),fAxis] = mtspectrumc(cueTargetDelayLfps, params);
        
        cueTargetDelayLfpsP3 = arrayOnsetLfpP3Current(cueTargetDelayInd,:);
        [cueTargetDelayPowerP3(lfpCount,:),fAxis] = mtspectrumc(cueTargetDelayLfpsP3, params);
        
        cueTargetDelayLfpsP1 = arrayOnsetLfpP1Current(cueTargetDelayInd,:);
        [cueTargetDelayPowerP1(lfpCount,:),fAxis] = mtspectrumc(cueTargetDelayLfpsP1, params);
    end
    
    clear EL;
end

%%
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-v%d.mat', outputDir, v);
save(saveFileName);

%%
outputDirOrig = outputDir;
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-v%d.mat', outputDir, v);
load(saveFileName);
outputDir = outputDirOrig;

%% plot baseline power
baselinePowerDPul = (baselinePower(isInDPulvinar,:));
baselinePowerVPul = (baselinePower(isInVPulvinar,:));

meanBaselinePowerDPul = median(baselinePowerDPul);
meanBaselinePowerVPul = median(baselinePowerVPul);
seBaselinePowerDPul = std(baselinePowerDPul) / sqrt(sum(isInDPulvinar));
seBaselinePowerVPul = std(baselinePowerVPul) / sqrt(sum(isInVPulvinar));

cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

figure; 
hold on;
fillH = jbfill(fAxis, meanBaselinePowerVPul - seBaselinePowerVPul, meanBaselinePowerVPul + seBaselinePowerVPul, ...
        vPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanBaselinePowerDPul - seBaselinePowerDPul, meanBaselinePowerDPul + seBaselinePowerDPul, ...
        dPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
plot(fAxis, meanBaselinePowerDPul, 'LineWidth', 2, 'Color', dPulCol);
plot(fAxis, meanBaselinePowerVPul, 'LineWidth', 2, 'Color', vPulCol);
xlim([5 80]);
% ylim([-34 -14]);

plotFileName = sprintf('%s/allSessions-baselinePower-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
% regardless of cue location
cueTargetDelayRelativePowerDPul = ((cueTargetDelayPowerAllLocs(isInDPulvinar,:))) ./ ((baselinePower(isInDPulvinar,:)));
cueTargetDelayRelativePowerVPul = ((cueTargetDelayPowerAllLocs(isInVPulvinar,:))) ./ ((baselinePower(isInVPulvinar,:)));

% HACK outlier removal
fInd = fAxis >= 10 & fAxis <= 80;
cueTargetDelayRelativePowerDPul(any(cueTargetDelayRelativePowerDPul(:,fInd) > 1.4 | cueTargetDelayRelativePowerDPul(:,fInd) < 0.6, 2),:) = [];
cueTargetDelayRelativePowerVPul(any(cueTargetDelayRelativePowerVPul(:,fInd) > 1.4 | cueTargetDelayRelativePowerVPul(:,fInd) < 0.6, 2),:) = [];

meanCueTargetDelayRelativePowerDPul = mean(cueTargetDelayRelativePowerDPul);
meanCueTargetDelayRelativePowerVPul = mean(cueTargetDelayRelativePowerVPul);
seCueTargetDelayRelativePowerDPul = std(cueTargetDelayRelativePowerDPul) / sqrt(sum(isInDPulvinar));
seCueTargetDelayRelativePowerVPul = std(cueTargetDelayRelativePowerVPul) / sqrt(sum(isInVPulvinar));

figure; 
hold on;
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerVPul - seCueTargetDelayRelativePowerVPul, meanCueTargetDelayRelativePowerVPul + seCueTargetDelayRelativePowerVPul, ...
        vPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerDPul - seCueTargetDelayRelativePowerDPul, meanCueTargetDelayRelativePowerDPul + seCueTargetDelayRelativePowerDPul, ...
        dPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
plot(fAxis, meanCueTargetDelayRelativePowerDPul, 'LineWidth', 2, 'Color', dPulCol);
plot(fAxis, meanCueTargetDelayRelativePowerVPul, 'LineWidth', 2, 'Color', vPulCol);
xlim([5 80]);
% ylim([1 1.03]);

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-allLocs-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
% regardless of cue location
% TODO refactor
cueTargetDelayRelativePowerDPul = (((cueTargetDelayPowerP3(isInDPulvinar,:))) ./ ((baselinePower(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerVPul = (((cueTargetDelayPowerP3(isInVPulvinar,:))) ./ ((baselinePower(isInVPulvinar,:))) - 1) * 100;

% HACK outlier removal
fInd = fAxis >= 10 & fAxis <= 80;
cueTargetDelayRelativePowerDPul(any(cueTargetDelayRelativePowerDPul(:,fInd) > 40 | cueTargetDelayRelativePowerDPul(:,fInd) < -40, 2),:) = [];
cueTargetDelayRelativePowerVPul(any(cueTargetDelayRelativePowerVPul(:,fInd) > 40 | cueTargetDelayRelativePowerVPul(:,fInd) < -40, 2),:) = [];

meanCueTargetDelayRelativePowerDPul = mean(cueTargetDelayRelativePowerDPul);
meanCueTargetDelayRelativePowerVPul = mean(cueTargetDelayRelativePowerVPul);
seCueTargetDelayRelativePowerDPul = std(cueTargetDelayRelativePowerDPul) / sqrt(sum(isInDPulvinar));
seCueTargetDelayRelativePowerVPul = std(cueTargetDelayRelativePowerVPul) / sqrt(sum(isInVPulvinar));

figure_tr_inch(7, 5); 
subaxis(1, 1, 1, 'MB', 0.14);
hold on;
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerVPul - seCueTargetDelayRelativePowerVPul, meanCueTargetDelayRelativePowerVPul + seCueTargetDelayRelativePowerVPul, ...
        vPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerDPul - seCueTargetDelayRelativePowerDPul, meanCueTargetDelayRelativePowerDPul + seCueTargetDelayRelativePowerDPul, ...
        dPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
h1 = plot(fAxis, meanCueTargetDelayRelativePowerDPul, 'LineWidth', 2, 'Color', dPulCol);
h2 = plot(fAxis, meanCueTargetDelayRelativePowerVPul, 'LineWidth', 2, 'Color', vPulCol);
legend([h1 h2], {' dPul', ' vPul'}, 'box', 'off', 'Location', 'SouthEast');
xlim([5 80]);
ylim([-20 0]);
xlabel('Frequency (Hz)');
ylabel('Percent Change in Power Rel. to Baseline');
set(gca, 'box', 'off');
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-P3-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
% regardless of cue location
cueTargetDelayRelativePowerDPul = ((cueTargetDelayPowerP1(isInDPulvinar,:))) ./ ((baselinePower(isInDPulvinar,:)));
cueTargetDelayRelativePowerVPul = ((cueTargetDelayPowerP1(isInVPulvinar,:))) ./ ((baselinePower(isInVPulvinar,:)));

% HACK outlier removal
fInd = fAxis >= 10 & fAxis <= 80;
cueTargetDelayRelativePowerDPul(any(cueTargetDelayRelativePowerDPul(:,fInd) > 1.4 | cueTargetDelayRelativePowerDPul(:,fInd) < 0.6, 2),:) = [];
cueTargetDelayRelativePowerVPul(any(cueTargetDelayRelativePowerVPul(:,fInd) > 1.4 | cueTargetDelayRelativePowerVPul(:,fInd) < 0.6, 2),:) = [];

meanCueTargetDelayRelativePowerDPul = mean(cueTargetDelayRelativePowerDPul);
meanCueTargetDelayRelativePowerVPul = mean(cueTargetDelayRelativePowerVPul);
seCueTargetDelayRelativePowerDPul = std(cueTargetDelayRelativePowerDPul) / sqrt(sum(isInDPulvinar));
seCueTargetDelayRelativePowerVPul = std(cueTargetDelayRelativePowerVPul) / sqrt(sum(isInVPulvinar));

figure; 
hold on;
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerVPul - seCueTargetDelayRelativePowerVPul, meanCueTargetDelayRelativePowerVPul + seCueTargetDelayRelativePowerVPul, ...
        vPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerDPul - seCueTargetDelayRelativePowerDPul, meanCueTargetDelayRelativePowerDPul + seCueTargetDelayRelativePowerDPul, ...
        dPulCol, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
plot(fAxis, meanCueTargetDelayRelativePowerDPul, 'LineWidth', 2, 'Color', dPulCol);
plot(fAxis, meanCueTargetDelayRelativePowerVPul, 'LineWidth', 2, 'Color', vPulCol);
xlim([5 80]);
ylim([0.8 1]);

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-P1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


%% plot power in cue-target delay relative to baseline
% regardless of cue location
cueTargetDelayRelativePowerP3 = ((cueTargetDelayPowerP3(isInDPulvinar,:))) ./ ((baselinePower(isInDPulvinar,:)));
cueTargetDelayRelativePowerP1 = ((cueTargetDelayPowerP1(isInDPulvinar,:))) ./ ((baselinePower(isInDPulvinar,:)));

% HACK outlier removal
fInd = fAxis >= 10 & fAxis <= 80;
goodSessions = any(cueTargetDelayRelativePowerP3(:,fInd) > 1.4 | cueTargetDelayRelativePowerP3(:,fInd) < 0.6, 2) & ...
        any(cueTargetDelayRelativePowerP1(:,fInd) > 1.4 | cueTargetDelayRelativePowerP1(:,fInd) < 0.6, 2);
cueTargetDelayRelativePowerP3(goodSessions,:) = [];
cueTargetDelayRelativePowerP1(goodSessions,:) = [];

meanCueTargetDelayRelativePowerP3 = median(cueTargetDelayRelativePowerP3);
meanCueTargetDelayRelativePowerP1 = median(cueTargetDelayRelativePowerP1);
seCueTargetDelayRelativePowerP3 = std(cueTargetDelayRelativePowerP3) / sqrt(sum(isInDPulvinar));
seCueTargetDelayRelativePowerP1 = std(cueTargetDelayRelativePowerP1) / sqrt(sum(isInDPulvinar));

p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

figure; 
hold on;
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerP1 - seCueTargetDelayRelativePowerP1, meanCueTargetDelayRelativePowerP1 + seCueTargetDelayRelativePowerP1, ...
        p1Col, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerP3 - seCueTargetDelayRelativePowerP3, meanCueTargetDelayRelativePowerP3 + seCueTargetDelayRelativePowerP3, ...
        p3Col, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
plot(fAxis, meanCueTargetDelayRelativePowerP3, 'LineWidth', 2, 'Color', p3Col);
plot(fAxis, meanCueTargetDelayRelativePowerP1, 'LineWidth', 2, 'Color', p1Col);
xlim([5 80]);
ylim([0.8 1]);

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-dPul-P3vsP1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


%% plot power in cue-target delay relative to baseline
% regardless of cue location
cueTargetDelayRelativePowerP3 = ((cueTargetDelayPowerP3(isInVPulvinar,:))) ./ ((baselinePower(isInVPulvinar,:)));
cueTargetDelayRelativePowerP1 = ((cueTargetDelayPowerP1(isInVPulvinar,:))) ./ ((baselinePower(isInVPulvinar,:)));

% HACK outlier removal
fInd = fAxis >= 10 & fAxis <= 80;
goodSessions = any(cueTargetDelayRelativePowerP3(:,fInd) > 1.4 | cueTargetDelayRelativePowerP3(:,fInd) < 0.6, 2) & ...
        any(cueTargetDelayRelativePowerP1(:,fInd) > 1.4 | cueTargetDelayRelativePowerP1(:,fInd) < 0.6, 2);
cueTargetDelayRelativePowerP3(goodSessions,:) = [];
cueTargetDelayRelativePowerP1(goodSessions,:) = [];

meanCueTargetDelayRelativePowerP3 = mean(cueTargetDelayRelativePowerP3);
meanCueTargetDelayRelativePowerP1 = mean(cueTargetDelayRelativePowerP1);
seCueTargetDelayRelativePowerP3 = std(cueTargetDelayRelativePowerP3) / sqrt(sum(isInVPulvinar));
seCueTargetDelayRelativePowerP1 = std(cueTargetDelayRelativePowerP1) / sqrt(sum(isInVPulvinar));

p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

figure; 
hold on;
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerP1 - seCueTargetDelayRelativePowerP1, meanCueTargetDelayRelativePowerP1 + seCueTargetDelayRelativePowerP1, ...
        p1Col, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanCueTargetDelayRelativePowerP3 - seCueTargetDelayRelativePowerP3, meanCueTargetDelayRelativePowerP3 + seCueTargetDelayRelativePowerP3, ...
        p3Col, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
plot(fAxis, meanCueTargetDelayRelativePowerP3, 'LineWidth', 2, 'Color', p3Col);
plot(fAxis, meanCueTargetDelayRelativePowerP1, 'LineWidth', 2, 'Color', p1Col);
xlim([5 80]);
ylim([0.8 1]);

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-vPul-P3vsP1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
