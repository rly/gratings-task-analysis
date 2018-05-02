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
end

%%
figure; 
hold on;
plot(fAxis, baselinePower(isInDPulvinar,:), 'LineWidth', 2);
plot(fAxis, baselinePower(isInVPulvinar,:), 'LineWidth', 2);
plot(fAxis, baselinePower, '--', 'LineWidth', 2);

plotFileName = sprintf('%s/%s-baselinePower-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

figure; 
hold on;
plot(fAxis, nanmean(cueTargetDelayPower(isInDPulvinar,:) ./ baselinePower(isInDPulvinar,:)), 'LineWidth', 2);
plot(fAxis, nanmean(cueTargetDelayPower(isInVPulvinar,:) ./ baselinePower(isInVPulvinar,:)), 'LineWidth', 2);
plot(fAxis, nanmean(cueTargetDelayPower ./ baselinePower), '--', 'LineWidth', 2);

plotFileName = sprintf('%s/%s-cueTargetDelayPower-v%d.png', outputDir, sessionName, sessionInd, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
