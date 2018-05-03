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
cueTargetDelayPower = nan(nLfpsApprox, 199);

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
        
        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        [baselinePower(lfpCount,:),fAxis] = mtspectrumc(preCueBaselineLfps, params);
        
        cueResponseLfps = cueOnsetLfpCurrent(cueResponseInd,:);
        [cueResponsePower(lfpCount,:),fAxis] = mtspectrumc(cueResponseLfps, params);
        
        cueTargetDelayLfps = arrayOnsetLfpCurrent(cueTargetDelayInd,:);
        [cueTargetDelayPower(lfpCount,:),fAxis] = mtspectrumc(cueTargetDelayLfps, params);
    end
    
    clear EL;
end

%%
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-v%d.mat', outputDir, v);
save(saveFileName);

%%

baselinePowerDPul = 10*log10(baselinePower(isInDPulvinar,:));
baselinePowerVPul = 10*log10(baselinePower(isInVPulvinar,:));

meanBaselinePowerDPul = mean(baselinePowerDPul);
meanBaselinePowerVPul = mean(baselinePowerVPul);
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

plotFileName = sprintf('%s/allSessions-baselinePower-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%
cueTargetDelayRelativePowerDPul = (10*log10(cueTargetDelayPower(isInDPulvinar,:))) ./ (10*log10(baselinePower(isInDPulvinar,:)));
cueTargetDelayRelativePowerVPul = (10*log10(cueTargetDelayPower(isInVPulvinar,:))) ./ (10*log10(baselinePower(isInVPulvinar,:)));

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

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
