function saveFileName = computeEvokedSpiking(saveFileName, spikeStruct, nLoc, UE, numRandomizations)
% spikeStruct = struct of spiking data
% UE = useful events struct

% compute SPDF around each event
% compute average firing rate, Fano Factor, peak rate in an analysis window
% using both SPDF values and raw spike counts
% determine InRF and ExRF locations
% compute max firing rate for normalization
% compute latency of responses to certain events
% compute task modulation stats
% - activity by location > baseline
% --- rank sum test
% --- resampled baseline test
% - visual activity > previous delay period
% --- sign rank test
% --- resampled baseline
% compute attention effect in CT delay, array hold response, TD delay
% - InRF > ExRF
% --- rank sum test
% --- permutation test
% --- attention index
% compute spatial selectivity (use all locs)
% - info rate > 0 
% --- permutation test


%%

spikeTs = spikeStruct.ts;
kernelSigma = 0.01;

clear spikeStruct;

nTrials = numel(UE.cueOnset);
nTrialsHoldBal = numel(UE.arrayOnsetHoldBal);
nTrialsByLoc = cellfun(@numel, UE.cueOnsetByLoc);
isLocUsed = nTrialsByLoc > 0;
fprintf('Computing evoked spiking with SPDF sigma %0.3f seconds and %d randomizations over %d trials...\n', ...
        kernelSigma, numRandomizations, nTrials);
fprintf('\t');

%% align spikes to cue onset, compute spdf
cueOnset.window = [0.8 0.8]; % seconds before, after
cueOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
cueOnsetRel = cueOnset; % copy
cueOnsetHold = cueOnset; % copy

cueOnset = createTimeLockedSpdf(spikeTs, UE.cueOnset, UE.cueOnsetByLoc, cueOnset, kernelSigma);
cueOnsetRel = createTimeLockedSpdf(spikeTs, UE.cueOnsetRel, UE.cueOnsetRelByLoc, cueOnsetRel, kernelSigma);
cueOnsetHold = createTimeLockedSpdf(spikeTs, UE.cueOnsetHold, UE.cueOnsetHoldByLoc, cueOnsetHold, kernelSigma);

%% align spikes to array onset, compute spdf
arrayOnset.window = [0.8 0.8]; % seconds before, after
arrayOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
arrayOnsetRel = arrayOnset; % copy
arrayOnsetHold = arrayOnset; % copy
arrayOnsetRelBal = arrayOnset; % copy
arrayOnsetHoldBal = arrayOnset; % copy

arrayOnset = createTimeLockedSpdf(spikeTs, UE.arrayOnset, UE.arrayOnsetByLoc, arrayOnset, kernelSigma);
arrayOnsetRel = createTimeLockedSpdf(spikeTs, UE.arrayOnsetRel, UE.arrayOnsetRelByLoc, arrayOnsetRel, kernelSigma);
arrayOnsetHold = createTimeLockedSpdf(spikeTs, UE.arrayOnsetHold, UE.arrayOnsetHoldByLoc, arrayOnsetHold, kernelSigma);
arrayOnsetRelBal = createTimeLockedSpdf(spikeTs, UE.arrayOnsetRelBal, UE.arrayOnsetRelBalByLoc, arrayOnsetRelBal, kernelSigma);
arrayOnsetHoldBal = createTimeLockedSpdf(spikeTs, UE.arrayOnsetHoldBal, UE.arrayOnsetHoldBalByLoc, arrayOnsetHoldBal, kernelSigma);

%% align spikes to target dimming, compute spdf
targetDim.window = [0.8 0.8]; % seconds before, after
targetDim.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
targetDimBalShortHold = targetDim;
targetDimBalLongHold = targetDim;
targetDimBal = targetDim;

targetDim = createTimeLockedSpdf(spikeTs, UE.targetDim, UE.targetDimByLoc, targetDim, kernelSigma);
targetDimBalShortHold = createTimeLockedSpdf(spikeTs, UE.targetDimBalShortHold, UE.targetDimBalShortHoldByLoc, targetDimBalShortHold, kernelSigma);
targetDimBalLongHold = createTimeLockedSpdf(spikeTs, UE.targetDimBalLongHold, UE.targetDimBalLongHoldByLoc, targetDimBalLongHold, kernelSigma);
targetDimBal = createTimeLockedSpdf(spikeTs, UE.targetDimBal, UE.targetDimBalByLoc, targetDimBal, kernelSigma);

%% look at lever-locked responses and saccade-locked responses
% note that this includes both hold and release trials
assert(~isempty(UE.fixationAndLeverTimes))

%% align spikes to enter/exit fixation, compute spdf
enterFixation.window = [0.8 0.8]; % seconds before, after
enterFixation.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
enterFixationRel = enterFixation; % copy
enterFixationHold = enterFixation; % copy

enterFixation = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCueByLoc, enterFixation, kernelSigma);
enterFixationRel = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationRelTimesPreCue, UE.fixationAndLeverTimes.firstEnterFixationRelTimesPreCueByLoc, enterFixationRel, kernelSigma);
enterFixationHold = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationHoldTimesPreCue, UE.fixationAndLeverTimes.firstEnterFixationHoldTimesPreCueByLoc, enterFixationHold, kernelSigma);

exitFixation.window = [0.8 0.8]; % seconds before, after
exitFixation.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
exitFixationLeft = exitFixation; % copy
exitFixationRight = exitFixation; % copy
exitFixationRel = exitFixation; % copy
exitFixationHold = exitFixation; % copy

exitFixation = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuiceByLoc, exitFixation, kernelSigma);
exitFixationLeft = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationLeftTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationLeftTimesAroundJuiceByLoc, exitFixationLeft, kernelSigma);
exitFixationRight = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationRightTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationRightTimesAroundJuiceByLoc, exitFixationRight, kernelSigma);
exitFixationRel = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationRelTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationRelTimesAroundJuiceByLoc, exitFixationRel, kernelSigma);
exitFixationHold = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationHoldTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationHoldTimesAroundJuiceByLoc, exitFixationHold, kernelSigma);

%% align spikes to lever press and release, compute spdf
leverPress.window = [0.8 0.8]; % seconds before, after
leverPress.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects

leverPress = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstLeverPressTimesPreCue, UE.fixationAndLeverTimes.firstLeverPressTimesPreCueByLoc, leverPress, kernelSigma);

leverRelease.window = [0.8 0.8]; % seconds before, after
leverRelease.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects

leverRelease = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuiceByLoc, leverRelease, kernelSigma);

%% calc time between motor events
assert(numel(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue) == numel(UE.fixationAndLeverTimes.firstLeverPressTimesPreCue));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice));

arrayOnsetRelToJuiceEventTime = UE.firstJuiceEvent(~UE.isHoldTrial) - UE.arrayOnsetRel; % almost same as UE.rt(~isHoldTrial)?
targetDimToJuiceEventTime = UE.firstJuiceEvent(UE.isHoldTrial) - UE.targetDim; % almost same as UE.rt(isHoldTrial)?
enterFixationToLeverPressTime = UE.fixationAndLeverTimes.firstLeverPressTimesPreCue - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
enterFixationToCueOnsetTime = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
exitFixationToJuiceEventTime = UE.firstJuiceEvent - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
exitFixationToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
arrayOnsetRelToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
targetDimToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;
arrayOnsetRelToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
targetDimToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;

%% analysis time windows
preLeverReleaseWindowOffset = [-0.175 0];
preCueBaselineWindowOffset = [-0.175 0];
cueResponseWindowOffset = [0.025 0.2];
cueTargetDelayWindowOffset = [-0.175 0];
cueTargetDelayLongWindowOffset = [-0.4 0];
arrayResponseWindowOffset = [0.025 0.2];
targetDimDelayWindowOffset = [-0.175 0];
targetDimDelayLongWindowOffset = [-0.4 0];
targetDimResponseWindowOffset = [0.025 0.2]; 
postTargetDimMotorResponseWindowOffset = [0.425 0.6];
preEnterFixationWindowOffset = [-0.175 0];
postEnterFixationWindowOffset = [0.025 0.2];
postEnterFixationLateWindowOffset = [0.15 0.325];
preExitFixationWindowOffset = [-0.05 0]; % note smaller window - TODO does this matter much?
preExitFixationEarlyWindowOffset = [-0.2 -0.05]; % note smaller window - TODO does this matter much?
postExitFixationWindowOffset = [0.0 0.1]; % note smaller window - TODO does this matter much?
postLeverReleaseWindowOffset = [0.025 0.2];

%% average number of spikes in each window
averageFiringRatesByCount = struct();

averageFiringRatesByCount.preEnterFixation = computeAverageFiringRateByCount(...
        preEnterFixationWindowOffset, enterFixation);
averageFiringRatesByCount.postEnterFixation = computeAverageFiringRateByCount(...
        postEnterFixationWindowOffset, enterFixation);
averageFiringRatesByCount.postEnterFixationLate = computeAverageFiringRateByCount(...
        postEnterFixationLateWindowOffset, enterFixation);
averageFiringRatesByCount.preCueBaseline = computeAverageFiringRateByCount(...
        preCueBaselineWindowOffset, cueOnset);
averageFiringRatesByCount.cueResponse = computeAverageFiringRateByCount(...
        cueResponseWindowOffset, cueOnset);
averageFiringRatesByCount.cueTargetDelay = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnset);
averageFiringRatesByCount.cueTargetDelayRelBal = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnsetRelBal);
averageFiringRatesByCount.cueTargetDelayHoldBal = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnsetHoldBal);
averageFiringRatesByCount.cueTargetDelayLong = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnset);
averageFiringRatesByCount.cueTargetDelayLongRelBal = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnsetRelBal);
averageFiringRatesByCount.cueTargetDelayLongHoldBal = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnsetHoldBal);
averageFiringRatesByCount.arrayResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnset);
averageFiringRatesByCount.arrayRelBalResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnsetRelBal);
averageFiringRatesByCount.arrayHoldBalResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnsetHoldBal);
averageFiringRatesByCount.targetDimBalDelay = computeAverageFiringRateByCount(...
        targetDimDelayWindowOffset, targetDimBal);
averageFiringRatesByCount.targetDimBalDelayLong = computeAverageFiringRateByCount(...
        targetDimDelayLongWindowOffset, targetDimBal);
averageFiringRatesByCount.targetDimBalResponse = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimBal);
averageFiringRatesByCount.targetDimBalShortHoldResponse = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimBalShortHold);
averageFiringRatesByCount.targetDimBalLongHoldResponse = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimBalLongHold);
averageFiringRatesByCount.postTargetDimBalMotorResponse = computeAverageFiringRateByCount(...
        postTargetDimMotorResponseWindowOffset, targetDimBal);
averageFiringRatesByCount.preExitFixation = computeAverageFiringRateByCount(...
        preExitFixationWindowOffset, exitFixation);
averageFiringRatesByCount.preExitFixationEarly = computeAverageFiringRateByCount(...
        preExitFixationEarlyWindowOffset, exitFixation);
averageFiringRatesByCount.postExitFixation = computeAverageFiringRateByCount(...
        postExitFixationWindowOffset, exitFixation);

%% average spdf value in each window
averageFiringRatesBySpdf = struct();
averageFiringRatesBySpdf.kernelSigma = kernelSigma;
averageFiringRatesBySpdf.preEnterFixation = computeAverageFiringRateBySpdf(...
        preEnterFixationWindowOffset, enterFixation);
averageFiringRatesBySpdf.postEnterFixation = computeAverageFiringRateBySpdf(...
        postEnterFixationWindowOffset, enterFixation);
averageFiringRatesBySpdf.postEnterFixationLate = computeAverageFiringRateBySpdf(...
        postEnterFixationLateWindowOffset, enterFixation);
averageFiringRatesBySpdf.preCueBaseline = computeAverageFiringRateBySpdf(...
        preCueBaselineWindowOffset, cueOnset);
averageFiringRatesBySpdf.cueResponse = computeAverageFiringRateBySpdf(...
        cueResponseWindowOffset, cueOnset);
averageFiringRatesBySpdf.cueTargetDelay = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnset);
averageFiringRatesBySpdf.cueTargetDelayRelBal = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnsetRelBal);
averageFiringRatesBySpdf.cueTargetDelayHoldBal = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnsetHoldBal);
averageFiringRatesBySpdf.cueTargetDelayLong = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnset);
averageFiringRatesBySpdf.cueTargetDelayLongRelBal = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnsetRelBal);
averageFiringRatesBySpdf.cueTargetDelayLongHoldBal = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnsetHoldBal);
averageFiringRatesBySpdf.arrayResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnset);
averageFiringRatesBySpdf.arrayRelBalResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnsetRelBal);
averageFiringRatesBySpdf.arrayHoldBalResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnsetHoldBal);
averageFiringRatesBySpdf.targetDimBalDelay = computeAverageFiringRateBySpdf(...
        targetDimDelayWindowOffset, targetDimBal);
averageFiringRatesBySpdf.targetDimBalDelayLong = computeAverageFiringRateBySpdf(...
        targetDimDelayLongWindowOffset, targetDimBal);
averageFiringRatesBySpdf.targetDimBalResponse = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimBal);
averageFiringRatesBySpdf.targetDimBalShortHoldResponse = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimBalShortHold);
averageFiringRatesBySpdf.targetDimBalLongHoldResponse = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimBalLongHold);
averageFiringRatesBySpdf.postTargetDimBalMotorResponse = computeAverageFiringRateBySpdf(...
        postTargetDimMotorResponseWindowOffset, targetDimBal);
averageFiringRatesBySpdf.preExitFixation = computeAverageFiringRateBySpdf(...
        preExitFixationWindowOffset, exitFixation);
averageFiringRatesBySpdf.preExitFixationEarly = computeAverageFiringRateBySpdf(...
        preExitFixationEarlyWindowOffset, exitFixation);
averageFiringRatesBySpdf.postExitFixation = computeAverageFiringRateBySpdf(...
        postExitFixationWindowOffset, exitFixation);

%% compute max firing rate across RF locations and conditions
fn = fieldnames(averageFiringRatesByCount);
maxFiringRateByCount = -Inf;
minFiringRateByCount = Inf;
for i = 1:numel(fn)
    maxFiringRate = max(averageFiringRatesByCount.(fn{i}).byLoc);
    minFiringRate = min(averageFiringRatesByCount.(fn{i}).byLoc);
    if maxFiringRate > maxFiringRateByCount
        maxFiringRateByCount = maxFiringRate;
    end
    if minFiringRate < minFiringRateByCount
        minFiringRateByCount = minFiringRate;
    end
end
clear fn maxFiringRate minFiringRate;

%% compute max firing rate across RF locations and conditions
% use spdf time courses
cueOnsetNormalizationWindowOffset = [-0.3 0.6];
arrayOnsetHoldNormalizationWindowOffset = [-0.6 0.65];
arrayOnsetRelNormalizationWindowOffset = [-0.6 0.25];
targetDimNormalizationWindowOffset = [-0.65 0.25];

cueOnsetNormalizationLogical = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + cueOnsetNormalizationWindowOffset);
arrayOnsetHoldNormalizationLogical = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + arrayOnsetHoldNormalizationWindowOffset);
arrayOnsetRelNormalizationLogical = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + arrayOnsetRelNormalizationWindowOffset);
targetDimNormalizationLogical = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimNormalizationWindowOffset);
% exclude motor response here

% chain time courses into long vectors for each location
allSpdfsByLoc = [cueOnset.spdfByLoc(:,cueOnsetNormalizationLogical) ...
        arrayOnset.spdfByLoc(:,arrayOnsetRelNormalizationLogical) ... % be conservative and use rel here
        arrayOnsetRelBal.spdfByLoc(:,arrayOnsetRelNormalizationLogical) ...
        arrayOnsetHoldBal.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
        targetDimBal.spdfByLoc(:,targetDimNormalizationLogical) ...
        targetDimBalShortHold.spdfByLoc(:,targetDimNormalizationLogical) ...
        targetDimBalLongHold.spdfByLoc(:,targetDimNormalizationLogical)];
maxFiringRateBySpdf = max(max(allSpdfsByLoc));
minFiringRateBySpdf = min(min(allSpdfsByLoc));

enterFixationInclMotorNormalizationWindowOffset = [-0.4 0.4];
cueOnsetInclMotorNormalizationWindowOffset = [-0.6 0.6];
arrayOnsetRelInclMotorNormalizationWindowOffset = [-0.6 0.6];
targetDimInclMotorNormalizationWindowOffset = [-0.65 0.6];
exitFixationInclMotorNormalizationWindowOffset = [-0.4 0.4];

enterFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(enterFixation.t, enterFixation.window(1) + enterFixationInclMotorNormalizationWindowOffset);
cueOnsetInclMotorNormalizationLogical = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + cueOnsetInclMotorNormalizationWindowOffset);
arrayOnsetRelInclMotorNormalizationLogical = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + arrayOnsetRelInclMotorNormalizationWindowOffset);
targetDimInclMotorNormalizationLogical = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimInclMotorNormalizationWindowOffset);
exitFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(exitFixation.t, exitFixation.window(1) + exitFixationInclMotorNormalizationWindowOffset);

allSpdfsByLoc = [enterFixation.spdfByLoc(:,enterFixationInclMotorNormalizationLogical) ...
        cueOnset.spdfByLoc(:,cueOnsetInclMotorNormalizationLogical) ...
        arrayOnset.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ... % be conservative and use hold here
        arrayOnsetRelBal.spdfByLoc(:,arrayOnsetRelInclMotorNormalizationLogical) ...
        arrayOnsetHoldBal.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
        targetDimBal.spdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
        targetDimBalShortHold.spdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
        targetDimBalLongHold.spdfByLoc(:,targetDimInclMotorNormalizationLogical)...
        exitFixation.spdfByLoc(:,exitFixationInclMotorNormalizationLogical)];
maxFiringRateBySpdfInclMotor = max(max(allSpdfsByLoc));
minFiringRateBySpdfInclMotor = min(min(allSpdfsByLoc));

clear allSpdfsByLoc;
clear cueOnsetNormalizationLogical arrayOnsetHoldNormalizationLogical ...
        arrayOnsetRelNormalizationLogical targetDimNormalizationLogical ...
        enterFixationInclMotorNormalizationLogical cueOnsetInclMotorNormalizationLogical ...
        arrayOnsetRelInclMotorNormalizationLogical targetDimInclMotorNormalizationLogical ...
        exitFixationInclMotorNormalizationLogical;

%% compute RF as largest mean cue response over baseline
meanCueResponseBaselineCorrByLoc = nan(nLoc, 1);
for i = 1:nLoc
    meanCueResponseBaselineCorrByLoc(i) = averageFiringRatesBySpdf.cueResponse.byLoc(i) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i);
end
[~,inRFLoc] = max(meanCueResponseBaselineCorrByLoc);
assert(nLoc == 4); % next line based on nLoc == 4
exRFLoc = mod(inRFLoc + 1, 4) + 1; % opposite location

%% compute RF by max per-condition z-scored cue response
cueResponseWindowIndices = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + cueResponseWindowOffset);
maxCueResponseBaselineCorrByLoc = nan(nLoc, 1);
extremeCueResponseBaselineCorrByLoc = nan(nLoc, 1);
for i = 1:nLoc
    zScoredCueResponseBaselineCorrByLoc = (cueOnset.spdfByLoc(i,cueResponseWindowIndices) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i)) / ...
            averageFiringRatesBySpdf.preCueBaseline.byLocSDOverTime(i);
    maxCueResponseBaselineCorrByLoc(i) = max(zScoredCueResponseBaselineCorrByLoc);
    extremeCueResponseBaselineCorrByLoc(i) = max(abs(zScoredCueResponseBaselineCorrByLoc));
end
[~,inRFLocByMax] = max(maxCueResponseBaselineCorrByLoc);
exRFLocByMax = mod(inRFLocByMax + 1, 4) + 1; % opposite location
[~,inRFLocByExtreme] = max(extremeCueResponseBaselineCorrByLoc);
exRFLocByExtreme = mod(inRFLocByExtreme + 1, 4) + 1; % opposite location

clear cueResponseWindowIndices zScoredCueResponseBaselineCorrByLoc ...
        maxCueResponseBaselineCorrByLoc extremeCueResponseBaselineCorrByLoc;

%% compute cue response at InRF > baseline
% bootstrap on mean of baseline SPDF with 500 shuffles
% compare actual mean SPDF cue response to distribution of bootstrapped
% pre-cue baselines
% TODO using count should be faster than psth
bootstrappedMeanPreCueBaselines = zeros(numRandomizations, 1);
preCueBaselineWindowIndices = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + preCueBaselineWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrials, 1, nTrials);
    bootstrapPsth = fixedPsth(cueOnset.spikeTimes(bootRandInd), kernelSigma, 0, cueOnset.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanPreCueBaselines(m) = mean(bootstrapPsth(preCueBaselineWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanPreCueBaselines = mean(bootstrappedMeanPreCueBaselines);
clear m bootRandInd bootstrapPsthResponse preCueBaselineWindowIndices;

% inRFLoc is defined as location with largest mean > baseline
% directional test
cueResponsePValueByBootstrapBaselineSpdfInRF = sum(averageFiringRatesBySpdf.cueResponse.byLoc(inRFLoc) < bootstrappedMeanPreCueBaselines) / numRandomizations;

%% compute cue response at any location ~= baseline
cueResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.cueResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        cueResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        cueResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute cue target delay at any location ~= baseline
cueTargetDelayPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute array response (HOLD) at any location ~= baseline
arrayHoldBalResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        arrayHoldBalResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        arrayHoldBalResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute target dim delay at any location ~= baseline
targetDimBalDelayPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimBalDelay.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        targetDimBalDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimBalDelay.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        targetDimBalDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimBalDelay.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute target dim response at any location ~= baseline
targetDimBalResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimBalResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        targetDimBalResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimBalResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        targetDimBalResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimBalResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute pre exit fixation response at any location ~= baseline
preExitFixationPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.preExitFixation.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        preExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixation.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        preExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixation.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute post exit fixation response at any location ~= baseline
postExitFixationPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.postExitFixation.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        postExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.postExitFixation.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        postExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.postExitFixation.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute array response at InRF > cue-target delay
% bootstrap on mean of cue-target delay SPDF with 500 shuffles
% compare actual mean SPDF array response to distribution of bootstrapped
% cue-target delay
% TODO using count should be faster than psth
bootstrappedMeanCueTargetDelays = zeros(numRandomizations, 1);
cueTargetDelayWindowIndices = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + cueTargetDelayWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrialsHoldBal, 1, nTrialsHoldBal);
    bootstrapPsth = fixedPsth(arrayOnset.spikeTimes(bootRandInd), kernelSigma, 0, arrayOnset.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanCueTargetDelays(m) = mean(bootstrapPsth(cueTargetDelayWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanCueTargetDelays = mean(bootstrappedMeanCueTargetDelays);
clear m bootRandInd bootstrapPsthResponse cueTargetDelayWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
arrayHoldBalResponsePValueByBootstrapCueTargetDelaySpdfInRF = sum(averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(inRFLoc) < bootstrappedMeanCueTargetDelays) / numRandomizations;

%% compute array response (HOLD) at any location ~= cue-target delay
arrayHoldBalResponsePValueByBootstrapCueTargetDelaySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(i) > meanBootstrappedMeanCueTargetDelays
        arrayHoldBalResponsePValueByBootstrapCueTargetDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(i) < bootstrappedMeanCueTargetDelays) / numRandomizations * 2;
    else
        arrayHoldBalResponsePValueByBootstrapCueTargetDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(i) > bootstrappedMeanCueTargetDelays) / numRandomizations * 2;
    end
end

%% compute target-dim response at InRF > target-dim delay
% bootstrap on mean of target-dim delay SPDF with 500 shuffles
% compare actual mean SPDF target-dim response to distribution of bootstrapped
% target-dim delay
% TODO using count should be faster than psth
bootstrappedMeanTargetDimBalDelays = zeros(numRandomizations, 1);
targetDimBalDelayWindowIndices = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimDelayWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrialsHoldBal, 1, nTrialsHoldBal);
    bootstrapPsth = fixedPsth(targetDimBal.spikeTimes(bootRandInd), kernelSigma, 0, targetDim.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanTargetDimBalDelays(m) = mean(bootstrapPsth(targetDimBalDelayWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanTargetDimBalDelays = mean(bootstrappedMeanTargetDimBalDelays);
clear m bootRandInd bootstrapPsthResponse targetDimBalDelayWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
targetDimBalResponsePValueByBootstrapTargetDimBalDelaySpdfInRF = sum(averageFiringRatesBySpdf.targetDimBalResponse.byLoc(inRFLoc) < bootstrappedMeanTargetDimBalDelays) / numRandomizations;

%% compute target dim response at any location ~= target-dim delay
targetDimBalResponsePValueByBootstrapTargetDimBalDelaySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimBalResponse.byLoc(i) > meanBootstrappedMeanTargetDimBalDelays
        targetDimBalResponsePValueByBootstrapTargetDimBalDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimBalResponse.byLoc(i) < bootstrappedMeanTargetDimBalDelays) / numRandomizations * 2;
    else
        targetDimBalResponsePValueByBootstrapTargetDimBalDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimBalResponse.byLoc(i) > bootstrappedMeanTargetDimBalDelays) / numRandomizations * 2;
    end
end

%% compute pre-exit fixation response at InRF > earlier pre-exit fixation window
% bootstrap on mean of pre-exit fixation early SPDF with 500 shuffles
% compare actual mean SPDF pre-exit fixation response to distribution of bootstrapped
% earlier pre-exit fixation responses
% TODO using count should be faster than psth
bootstrappedMeanPreExitFixationEarlys = zeros(numRandomizations, 1);
preExitFixationEarlyWindowIndices = getTimeLogicalWithTolerance(exitFixation.t, exitFixation.window(1) + preExitFixationEarlyWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrials, 1, nTrials);
    bootstrapPsth = fixedPsth(exitFixation.spikeTimes(bootRandInd), kernelSigma, 0, exitFixation.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanPreExitFixationEarlys(m) = mean(bootstrapPsth(preExitFixationEarlyWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanPreExitFixationEarlys = mean(bootstrappedMeanPreExitFixationEarlys);
clear m bootRandInd bootstrapPsthResponse preExitFixationEarlyWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
preExitFixationPValueByBootstrapPreExitFixationEarlySpdfInRF = sum(averageFiringRatesBySpdf.preExitFixationEarly.byLoc(inRFLoc) < bootstrappedMeanPreExitFixationEarlys) / numRandomizations;

%% compute pre-exit fixation response at any location ~= earlier pre-exit fixation window
preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.preExitFixationEarly.byLoc(i) > meanBootstrappedMeanPreExitFixationEarlys
        preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixationEarly.byLoc(i) < bootstrappedMeanPreExitFixationEarlys) / numRandomizations * 2;
    else
        preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixationEarly.byLoc(i) > bootstrappedMeanPreExitFixationEarlys) / numRandomizations * 2;
    end
end

%% rank-sum test on x vs baseline using individual trial spike counts
% spike rates per trial are not normally distributed
% they may also not have equal variance
% use nonparametric or resampling methods
% in fact, spike rates per trial are quantized, and because the length of
% the time periods differ between the baseline and the other periods,
% distributions may have different means (or maybe even medians) when that
% is due to the low resolution of the data

% could change the analysis windows to be equal in duration -- which should
% make the periods interchangeable, and make the rank sum test or
% permutation tests work

% also, consider a paired test may be more powerful TODO test this

% does it matter much that the cue response by location is compared to all
% the baseline trials together? shouldn't matter much. more data here, but 
% lose statistical power in that it's not a paired test
cueResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.cueResponse.trialRateByLoc);
cueTargetDelayVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc);
arrayHoldBalResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.arrayHoldBalResponse.trialRateByLoc);
targetDimBalDelayVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.targetDimBalDelay.trialRateByLoc);
targetDimBalResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.targetDimBalResponse.trialRateByLoc);
preExitFixationVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.preExitFixation.trialRateByLoc);
postExitFixationVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.postExitFixation.trialRateByLoc);

% by location -- use sign-rank test here
cueResponseVsBaselineSignRankTestStatsByLoc = computeSignRankTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRateByLoc, ...
        averageFiringRatesByCount.cueResponse.trialRateByLoc);
arrayHoldBalResponseVsCTDelayHoldBalSignRankTestStatsByLoc = computeSignRankTestByLoc(averageFiringRatesByCount.cueTargetDelayHoldBal.trialRateByLoc, ...
        averageFiringRatesByCount.arrayHoldBalResponse.trialRateByLoc);
targetDimBalResponseVsTargetDimBalDelaySignRankTestStatsByLoc = computeSignRankTestByLoc(averageFiringRatesByCount.targetDimBalDelay.trialRateByLoc, ...
        averageFiringRatesByCount.targetDimBalResponse.trialRateByLoc);
preExitFixationVsPreExitFixationEarlySignRankTestStatsByLoc = computeSignRankTestByLoc(averageFiringRatesByCount.preExitFixationEarly.trialRateByLoc, ...
        averageFiringRatesByCount.preExitFixation.trialRateByLoc);

%% shuffle test on cue onset time as another test for cue-evoked change in activity
% concatenate the two time periods in a trial and circularly permute the
% event time across time bins 
% I think the event response time period has to be reasonably short (e.g.
% doesn't have contamination from other evoked responses), but the
% preceding time can be however long. 
% This is similar to the permutation test of having group 1 being mean 
% baseline activity over N1 trials and group 2 being mean cue response 
% activity over N2 trials, pooling them and assigning N1 trials to group 1
% and N2 trials ot group 2, and re-computing the mean activity across
% trials P times. 

%% permutation test on cue-target delay period
shuffleCueTargetDelayDiff = zeros(numRandomizations, 1);
shuffleCueTargetDelayAI = zeros(numRandomizations, 1);
cueTargetDelayWindowIndices = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + cueTargetDelayWindowOffset);
nTrialsInRF = numel(arrayOnset.spikeTimesByLoc{inRFLoc});
nTrialsExRF = numel(arrayOnset.spikeTimesByLoc{exRFLoc});
nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
arrayOnsetHoldBalSpikeTimesInRFExRF = [arrayOnset.spikeTimesByLoc{inRFLoc} arrayOnset.spikeTimesByLoc{exRFLoc}];
for m = 1:numRandomizations
    shuffleIndices = randperm(nTrialsInRFExRF);
    randSpikeTimesInRF = arrayOnsetHoldBalSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
    randSpikeTimesExRF = arrayOnsetHoldBalSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
    
    shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, arrayOnset.t);
    shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, arrayOnset.t);
    meanResponseShufflePsthInRF = mean(shufflePsthInRF(cueTargetDelayWindowIndices));
    meanResponseShufflePsthExRF = mean(shufflePsthExRF(cueTargetDelayWindowIndices));
    shuffleCueTargetDelayDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
    shuffleCueTargetDelayAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
            (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
end
fprintf('.');
cueTargetDelayDiff = averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc);
cueTargetDelayAI = (averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc)) / ...
        (averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) + averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc));
% inRFLoc is defined as location with largest mean > baseline
if cueTargetDelayDiff > mean(shuffleCueTargetDelayDiff)
    cueTargetDelayDiffPValueByShuffleSpdf = sum(cueTargetDelayDiff < shuffleCueTargetDelayDiff) / numRandomizations * 2;
else
    cueTargetDelayDiffPValueByShuffleSpdf = sum(cueTargetDelayDiff > shuffleCueTargetDelayDiff) / numRandomizations * 2;
end

if cueTargetDelayAI > mean(shuffleCueTargetDelayAI)
    cueTargetDelayAIPValueByShuffleSpdf = sum(cueTargetDelayAI < shuffleCueTargetDelayAI) / numRandomizations * 2;
else
    cueTargetDelayAIPValueByShuffleSpdf = sum(cueTargetDelayAI > shuffleCueTargetDelayAI) / numRandomizations * 2;
end
clear cueTargetDelayWindowIndices arrayOnsetSpikeTimesInRFExRF;

%% rank sum test on cue-target delay period using individual trial spike counts
cueTargetDelayDiffPValueByRankSumTest = ranksum(...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{exRFLoc});

%% permutation test on array onset response period - hold balanced trials only
shuffleArrayHoldBalResponseDiff = zeros(numRandomizations, 1);
shuffleArrayHoldBalResponseAI = zeros(numRandomizations, 1);
arrayHoldBalResponseWindowIndices = getTimeLogicalWithTolerance(arrayOnsetHoldBal.t, arrayOnsetHoldBal.window(1) + arrayResponseWindowOffset);
nTrialsInRF = numel(arrayOnsetHoldBal.spikeTimesByLoc{inRFLoc});
nTrialsExRF = numel(arrayOnsetHoldBal.spikeTimesByLoc{exRFLoc});
nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
arrayOnsetHoldBalSpikeTimesInRFExRF = [arrayOnsetHoldBal.spikeTimesByLoc{inRFLoc} arrayOnsetHoldBal.spikeTimesByLoc{exRFLoc}];
for m = 1:numRandomizations
    shuffleIndices = randperm(nTrialsInRFExRF);
    randSpikeTimesInRF = arrayOnsetHoldBalSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
    randSpikeTimesExRF = arrayOnsetHoldBalSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
    
    shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, arrayOnsetHoldBal.t);
    shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, arrayOnsetHoldBal.t);
    meanResponseShufflePsthInRF = mean(shufflePsthInRF(arrayHoldBalResponseWindowIndices));
    meanResponseShufflePsthExRF = mean(shufflePsthExRF(arrayHoldBalResponseWindowIndices));
    shuffleArrayHoldBalResponseDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
    shuffleArrayHoldBalResponseAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
            (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
end
fprintf('.');
arrayHoldBalResponseDiff = averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(inRFLoc) - averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(exRFLoc);
arrayHoldBalResponseAI = (averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(inRFLoc) - averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(exRFLoc)) / ...
        (averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(inRFLoc) + averageFiringRatesBySpdf.arrayHoldBalResponse.byLoc(exRFLoc));
% inRFLoc is defined as location with largest mean > baseline
if arrayHoldBalResponseDiff > mean(shuffleArrayHoldBalResponseDiff)
    arrayHoldBalResponseDiffPValueByShuffleSpdf = sum(arrayHoldBalResponseDiff < shuffleArrayHoldBalResponseDiff) / numRandomizations * 2;
else
    arrayHoldBalResponseDiffPValueByShuffleSpdf = sum(arrayHoldBalResponseDiff > shuffleArrayHoldBalResponseDiff) / numRandomizations * 2;
end

if arrayHoldBalResponseAI > mean(shuffleArrayHoldBalResponseAI)
    arrayHoldBalResponseAIPValueByShuffleSpdf = sum(arrayHoldBalResponseAI < shuffleArrayHoldBalResponseAI) / numRandomizations * 2;
else
    arrayHoldBalResponseAIPValueByShuffleSpdf = sum(arrayHoldBalResponseAI > shuffleArrayHoldBalResponseAI) / numRandomizations * 2;
end
clear arrayHoldBalResponseWindowIndices arrayOnsetHoldBalSpikeTimesInRFExRF;

%% rank sum test on array hold response using individual trial spike counts
arrayHoldBalResponseDiffPValueByRankSumTest = ranksum(...
        averageFiringRatesByCount.arrayHoldBalResponse.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.arrayHoldBalResponse.trialRateByLoc{exRFLoc});
    
%% permutation test on target-dim delay period
shuffleTargetDimBalDelayDiff = zeros(numRandomizations, 1);
shuffleTargetDimBalDelayAI = zeros(numRandomizations, 1);
targetDimBalDelayWindowIndices = getTimeLogicalWithTolerance(targetDimBal.t, targetDimBal.window(1) + targetDimDelayWindowOffset);
nTrialsInRF = numel(targetDimBal.spikeTimesByLoc{inRFLoc});
nTrialsExRF = numel(targetDimBal.spikeTimesByLoc{exRFLoc});
nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
targetDimBalSpikeTimesInRFExRF = [targetDimBal.spikeTimesByLoc{inRFLoc} targetDimBal.spikeTimesByLoc{exRFLoc}];
for m = 1:numRandomizations    
    shuffleIndices = randperm(nTrialsInRFExRF);
    randSpikeTimesInRF = targetDimBalSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
    randSpikeTimesExRF = targetDimBalSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
    
    shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, targetDimBal.t);
    shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, targetDimBal.t);
    meanResponseShufflePsthInRF = mean(shufflePsthInRF(targetDimBalDelayWindowIndices));
    meanResponseShufflePsthExRF = mean(shufflePsthExRF(targetDimBalDelayWindowIndices));
    shuffleTargetDimBalDelayDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
    shuffleTargetDimBalDelayAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
            (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
end
fprintf('.');
targetDimBalDelayDiff = averageFiringRatesBySpdf.targetDimBalDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.targetDimBalDelay.byLoc(exRFLoc);
targetDimBalDelayAI = (averageFiringRatesBySpdf.targetDimBalDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.targetDimBalDelay.byLoc(exRFLoc)) / ...
        (averageFiringRatesBySpdf.targetDimBalDelay.byLoc(inRFLoc) + averageFiringRatesBySpdf.targetDimBalDelay.byLoc(exRFLoc));
% inRFLoc is defined as location with largest mean > baseline
if targetDimBalDelayDiff > mean(shuffleTargetDimBalDelayDiff)
    targetDimBalDelayDiffPValueByShuffleSpdf = sum(targetDimBalDelayDiff < shuffleTargetDimBalDelayDiff) / numRandomizations * 2;
else
    targetDimBalDelayDiffPValueByShuffleSpdf = sum(targetDimBalDelayDiff > shuffleTargetDimBalDelayDiff) / numRandomizations * 2;
end

if targetDimBalDelayAI > mean(shuffleTargetDimBalDelayAI)
    targetDimBalDelayAIPValueByShuffleSpdf = sum(targetDimBalDelayAI < shuffleTargetDimBalDelayAI) / numRandomizations * 2;
else
    targetDimBalDelayAIPValueByShuffleSpdf = sum(targetDimBalDelayAI > shuffleTargetDimBalDelayAI) / numRandomizations * 2;
end
clear targetDimBalDelayWindowIndices targetDimBalSpikeTimesInRFExRF;

%% rank-sum test on target-dim delay period using individual trial spike counts
targetDimBalDelayDiffPValueByRankSumTest = ranksum(...
        averageFiringRatesByCount.targetDimBalDelay.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.targetDimBalDelay.trialRateByLoc{exRFLoc});

%% spatial selectivity during pre-cue baseline using info rate (control)
preCueBaselineInfoRateStruct = computeInfoRatePValueByShuffle(...
        cueOnset, averageFiringRatesBySpdf.preCueBaseline, preCueBaselineWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during cue response using info rate
cueResponseInfoRateStruct = computeInfoRatePValueByShuffle(...
        cueOnset, averageFiringRatesBySpdf.cueResponse, cueResponseWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during cue target delay period using info rate
cueTargetDelayInfoRateStruct = computeInfoRatePValueByShuffle(...
        arrayOnset, averageFiringRatesBySpdf.cueTargetDelay, cueTargetDelayWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during array response using info rate
arrayHoldBalResponseInfoRateStruct = computeInfoRatePValueByShuffle(...
        arrayOnsetHoldBal, averageFiringRatesBySpdf.arrayHoldBalResponse, arrayResponseWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during target dim delay period using info rate
targetDimBalDelayInfoRateStruct = computeInfoRatePValueByShuffle(...
        targetDimBal, averageFiringRatesBySpdf.targetDimBalDelay, targetDimDelayWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during target dimming using info rate
targetDimBalResponseInfoRateStruct = computeInfoRatePValueByShuffle(...
        targetDimBal, averageFiringRatesBySpdf.targetDimBalResponse, targetDimResponseWindowOffset, numRandomizations);
fprintf('.');

%% latency of response by location
cueOnset = computeResponseLatencyByLoc(cueOnset, isLocUsed);
arrayOnsetRelBal = computeResponseLatencyByLoc(arrayOnsetRelBal, isLocUsed);
arrayOnsetHoldBal = computeResponseLatencyByLoc(arrayOnsetHoldBal, isLocUsed);
targetDimBal = computeResponseLatencyByLoc(targetDimBal, isLocUsed);

%% clean up unused locations
% temp until the change is made in computeEvokedSpiking
cueResponsePValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
arrayHoldBalResponsePValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
targetDimBalDelayPValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
targetDimBalResponsePValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
preExitFixationPValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
arrayHoldBalResponsePValueByBootstrapCueTargetDelaySpdfByLoc(~isLocUsed) = NaN;
targetDimBalResponsePValueByBootstrapTargetDimBalDelaySpdfByLoc(~isLocUsed) = NaN;
preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc(~isLocUsed) = NaN;

unusedLocs = find(~isLocUsed);
for i = 1:numel(unusedLocs)
    k = unusedLocs(i);
    if k <= nLoc
        cueResponseVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        cueTargetDelayVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        arrayHoldBalResponseVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        targetDimBalDelayVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        targetDimBalResponseVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        preExitFixationVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        cueResponseVsBaselineSignRankTestStatsByLoc(k).p = NaN;
        arrayHoldBalResponseVsCTDelayHoldBalSignRankTestStatsByLoc(k).p = NaN;
        targetDimBalResponseVsTargetDimBalDelaySignRankTestStatsByLoc(k).p = NaN;
        preExitFixationVsPreExitFixationEarlySignRankTestStatsByLoc(k).p = NaN;
    end
end

%% save
clear i m k;
clear nTrialsInRF nTrialsExRF nTrialsInRFExRF shuffleIndices;
clear randSpikeTimesInRF randSpikeTimesExRF shufflePsthInRF shufflePsthExRF;
clear meanResponseShufflePsthInRF meanResponseShufflePsthExRF;

% save all workspace variables not cleared above into mat file
fprintf('\n');
fprintf('Saving to %s.\n', saveFileName);
save(saveFileName);

end
