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

% Match spikeTimes with first fixation event
if isnan(spikeStruct.unitStartTime)
    startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>spikeTs(1), 1, 'first'));
else
    startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>spikeStruct.unitStartTime, 1, 'first'));
end
% Match spikeTimes with last lever release
if isnan(spikeStruct.unitEndTime)
    endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<spikeTs(end),1,'last'));
else
    endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<spikeStruct.unitEndTime,1,'last'));
end
assert(startTime < endTime)
    
kernelSigma = 0.01;

clear spikeStruct;

nTrials = numel(UE.cueOnset);
nTrialsHoldBal = numel(UE.arrayOnsetHoldBal);
nTrialsByLoc = cellfun(@numel, UE.cueOnsetByLoc);
isLocUsed = nTrialsByLoc > 0;
unusedLocs = find(~isLocUsed);

fprintf('Computing evoked spiking with SPDF sigma %0.3f seconds and %d randomizations over %d trials...\n', ...
        kernelSigma, numRandomizations, nTrials);
fprintf('\t');

%%
isSpikeTsInTrial = zeros(size(spikeTs));
for i = 1:nTrials
    inTrial = spikeTs >= UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(i) & ...
              spikeTs < UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(i);
    isSpikeTsInTrial = isSpikeTsInTrial | inTrial;
end
clear inTrial;

%% align spikes to cue onset, compute spdf
cueOnset.window = [0.8 0.8]; % seconds before, after
cueOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
% cueOnsetRel = cueOnset; % copy
% cueOnsetHold = cueOnset; % copy
cueOnsetError = cueOnset; % copy

cueOnset = createTimeLockedSpdf(spikeTs, UE.cueOnset, UE.cueOnsetByLoc, cueOnset, kernelSigma, startTime, endTime);
% cueOnsetRel = createTimeLockedSpdf(spikeTs, UE.cueOnsetRel, UE.cueOnsetRelByLoc, cueOnsetRel, kernelSigma);
% cueOnsetHold = createTimeLockedSpdf(spikeTs, UE.cueOnsetHold, UE.cueOnsetHoldByLoc, cueOnsetHold, kernelSigma);
cueOnsetError = createTimeLockedSpdf(spikeTs, UE.cueOnsetError, UE.cueOnsetErrorByLoc, cueOnsetError, kernelSigma, startTime, endTime);

fprintf('Using %d valid trials out of %d total trials.\n', numel(cueOnset.validEvents), nTrials)

%% align spikes to array onset, compute spdf
arrayOnset.window = [0.8 0.8]; % seconds before, after
arrayOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
% arrayOnsetRel = arrayOnset; % copy
% arrayOnsetHold = arrayOnset; % copy
arrayOnsetRelBal = arrayOnset; % copy
arrayOnsetHoldBal = arrayOnset; % copy
arrayOnsetError = arrayOnset;
arrayOnsetRelBalError = arrayOnset; % copy
arrayOnsetHoldBalError = arrayOnset; % copy

arrayOnset = createTimeLockedSpdf(spikeTs, UE.arrayOnset, UE.arrayOnsetByLoc, arrayOnset, kernelSigma, startTime, endTime);
% arrayOnsetRel = createTimeLockedSpdf(spikeTs, UE.arrayOnsetRel, UE.arrayOnsetRelByLoc, arrayOnsetRel, kernelSigma);
% arrayOnsetHold = createTimeLockedSpdf(spikeTs, UE.arrayOnsetHold, UE.arrayOnsetHoldByLoc, arrayOnsetHold, kernelSigma);
arrayOnsetRelBal = createTimeLockedSpdf(spikeTs, UE.arrayOnsetRelBal, UE.arrayOnsetRelBalByLoc, arrayOnsetRelBal, kernelSigma, startTime, endTime);
arrayOnsetHoldBal = createTimeLockedSpdf(spikeTs, UE.arrayOnsetHoldBal, UE.arrayOnsetHoldBalByLoc, arrayOnsetHoldBal, kernelSigma, startTime, endTime);
arrayOnsetError = createTimeLockedSpdf(spikeTs, UE.arrayOnsetError, UE.arrayOnsetErrorByLoc, arrayOnsetError, kernelSigma, startTime, endTime);
arrayOnsetRelBalError = createTimeLockedSpdf(spikeTs, UE.arrayOnsetRelBalError, UE.arrayOnsetRelBalErrorByLoc, arrayOnsetRelBalError, kernelSigma, startTime, endTime);
arrayOnsetHoldBalError = createTimeLockedSpdf(spikeTs, UE.arrayOnsetHoldBalError, UE.arrayOnsetHoldBalErrorByLoc, arrayOnsetHoldBalError, kernelSigma, startTime, endTime);

%% align spikes to target dimming, compute spdf
targetDimBal.window = [0.8 0.8]; % seconds before, after
targetDimBal.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
targetDimShortHoldBal = targetDimBal; % copy
targetDimLongHoldBal = targetDimBal; % copy

% targetDim = createTimeLockedSpdf(spikeTs, UE.targetDim, UE.targetDimByLoc, targetDim, kernelSigma);
targetDimShortHoldBal = createTimeLockedSpdf(spikeTs, UE.targetDimShortHoldBal, UE.targetDimShortHoldBalByLoc, targetDimShortHoldBal, kernelSigma, startTime, endTime);
targetDimLongHoldBal = createTimeLockedSpdf(spikeTs, UE.targetDimLongHoldBal, UE.targetDimLongHoldBalByLoc, targetDimLongHoldBal, kernelSigma, startTime, endTime);
targetDimBal = createTimeLockedSpdf(spikeTs, UE.targetDimBal, UE.targetDimBalByLoc, targetDimBal, kernelSigma, startTime, endTime);

%% look at lever-locked responses and saccade-locked responses
% note that this includes both hold and release trials
assert(~isempty(UE.fixationAndLeverTimes))

%% align spikes to enter/exit fixation, compute spdf
enterFixation.window = [0.8 0.8]; % seconds before, after
enterFixation.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
enterFixationRel = enterFixation; % copy
enterFixationHold = enterFixation; % copy

enterFixation = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue, ...
        UE.fixationAndLeverTimes.firstEnterFixationTimesPreCueByLoc, enterFixation, kernelSigma, startTime, endTime);
enterFixationRel = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationRelTimesPreCue, ...
        UE.fixationAndLeverTimes.firstEnterFixationRelTimesPreCueByLoc, enterFixationRel, kernelSigma, startTime, endTime);
enterFixationHold = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationHoldTimesPreCue, ...
        UE.fixationAndLeverTimes.firstEnterFixationHoldTimesPreCueByLoc, enterFixationHold, kernelSigma, startTime, endTime);

exitFixation.window = [0.8 0.8]; % seconds before, after
exitFixation.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
exitFixationLeft = exitFixation; % copy
exitFixationRight = exitFixation; % copy
exitFixationRel = exitFixation; % copy
exitFixationHold = exitFixation; % copy

exitFixation = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice, ...
        UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuiceByLoc, exitFixation, kernelSigma, startTime, endTime);
exitFixationLeft = exitFixation; %createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationLeftTimesAroundJuice, ...
%         UE.fixationAndLeverTimes.firstExitFixationLeftTimesAroundJuiceByLoc, exitFixationLeft, kernelSigma, startTime, endTime);
exitFixationRight = exitFixation; %createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationRightTimesAroundJuice, ...
%         UE.fixationAndLeverTimes.firstExitFixationRightTimesAroundJuiceByLoc, exitFixationRight, kernelSigma, startTime, endTime);
exitFixationRel = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationRelTimesAroundJuice, ...
        UE.fixationAndLeverTimes.firstExitFixationRelTimesAroundJuiceByLoc, exitFixationRel, kernelSigma, startTime, endTime);
exitFixationHold = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationHoldTimesAroundJuice, ...
        UE.fixationAndLeverTimes.firstExitFixationHoldTimesAroundJuiceByLoc, exitFixationHold, kernelSigma, startTime, endTime);

%% align spikes to lever press and release, compute spdf
leverPress.window = [0.8 0.8]; % seconds before, after
leverPress.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects

leverPress = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstLeverPressTimesPreCue, UE.fixationAndLeverTimes.firstLeverPressTimesPreCueByLoc, leverPress, kernelSigma, startTime, endTime);

leverRelease.window = [0.8 0.8]; % seconds before, after
leverRelease.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects

leverRelease = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuiceByLoc, leverRelease, kernelSigma, startTime, endTime);

%% calc time between motor events
assert(numel(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue) == numel(UE.fixationAndLeverTimes.firstLeverPressTimesPreCue));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice));

arrayOnsetRelToJuiceEventTime = UE.firstJuiceEvent(UE.isRelBal) - UE.arrayOnsetRelBal; % almost same as UE.rt(~isHoldTrial)?
targetDimToJuiceEventTime = UE.firstJuiceEvent(UE.isHoldBal) - UE.targetDimBal; % almost same as UE.rt(isHoldTrial)?
enterFixationToLeverPressTime = UE.fixationAndLeverTimes.firstLeverPressTimesPreCue - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
enterFixationToCueOnsetTime = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
exitFixationToJuiceEventTime = UE.firstJuiceEvent - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
exitFixationToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
arrayOnsetRelToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(UE.isRelBal) - UE.arrayOnsetRelBal;
targetDimToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(UE.isHoldBal) - UE.targetDimBal;
arrayOnsetRelToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(UE.isRelBal) - UE.arrayOnsetRelBal;
targetDimToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(UE.isHoldBal) - UE.targetDimBal;

%% analysis time windows
preLeverReleaseWindowOffset = [-0.175 0];
preCueBaselineWindowOffset = [-0.175 0];
cueResponseWindowOffset = [0.025 0.2];
cueTargetDelayWindowOffset = [-0.175 0];
cueTargetDelayLongWindowOffset = [-0.4 0];
arrayResponseWindowOffset = [0.025 0.2];
arrayResponseMidWindowOffset = [0.1 0.2];
arrayResponseLateWindowOffset = [0.2 0.375];
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
        preCueBaselineWindowOffset, cueOnset, cueOnsetError);
averageFiringRatesByCount.cueResponse = computeAverageFiringRateByCount(...
        cueResponseWindowOffset, cueOnset, cueOnsetError);
averageFiringRatesByCount.cueTargetDelay = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnset, arrayOnsetError);
averageFiringRatesByCount.cueTargetDelayRelBal = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnsetRelBal, arrayOnsetRelBalError);
averageFiringRatesByCount.cueTargetDelayHoldBal = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesByCount.cueTargetDelayLong = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnset, arrayOnsetError);
averageFiringRatesByCount.cueTargetDelayRelBalLong = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnsetRelBal, arrayOnsetRelBalError);
averageFiringRatesByCount.cueTargetDelayHoldBalLong = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesByCount.arrayResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnset, arrayOnsetError);
averageFiringRatesByCount.arrayResponseRelBal = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnsetRelBal, arrayOnsetRelBalError);
averageFiringRatesByCount.arrayResponseHoldBal = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesByCount.arrayResponseHoldMidBal = computeAverageFiringRateByCount(...
        arrayResponseMidWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesByCount.arrayResponseHoldLateBal = computeAverageFiringRateByCount(...
        arrayResponseLateWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesByCount.targetDimDelayBal = computeAverageFiringRateByCount(...
        targetDimDelayWindowOffset, targetDimBal);
averageFiringRatesByCount.targetDimDelayBalLong = computeAverageFiringRateByCount(...
        targetDimDelayLongWindowOffset, targetDimBal);
averageFiringRatesByCount.targetDimResponseBal = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimBal);
averageFiringRatesByCount.targetDimResponseShortHoldBal = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimShortHoldBal);
averageFiringRatesByCount.targetDimResponseLongHoldBal = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimLongHoldBal);
averageFiringRatesByCount.postTargetDimMotorResponseBal = computeAverageFiringRateByCount(...
        postTargetDimMotorResponseWindowOffset, targetDimBal);
averageFiringRatesByCount.preExitFixation = computeAverageFiringRateByCount(...
        preExitFixationWindowOffset, exitFixation);
averageFiringRatesByCount.preExitFixationEarly = computeAverageFiringRateByCount(...
        preExitFixationEarlyWindowOffset, exitFixation);

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
        preCueBaselineWindowOffset, cueOnset, cueOnsetError);
averageFiringRatesBySpdf.cueResponse = computeAverageFiringRateBySpdf(...
        cueResponseWindowOffset, cueOnset, cueOnsetError);
averageFiringRatesBySpdf.cueTargetDelay = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnset, arrayOnsetError);
averageFiringRatesBySpdf.cueTargetDelayRelBal = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnsetRelBal, arrayOnsetRelBalError);
averageFiringRatesBySpdf.cueTargetDelayHoldBal = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesBySpdf.cueTargetDelayLong = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnset, arrayOnsetError);
averageFiringRatesBySpdf.cueTargetDelayRelBalLong = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnsetRelBal, arrayOnsetRelBalError);
averageFiringRatesBySpdf.cueTargetDelayHoldBalLong = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesBySpdf.arrayResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnset, arrayOnsetError);
averageFiringRatesBySpdf.arrayResponseRelBal = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnsetRelBal, arrayOnsetRelBalError);
averageFiringRatesBySpdf.arrayResponseHoldBal = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesBySpdf.arrayResponseHoldMidBal = computeAverageFiringRateBySpdf(...
        arrayResponseMidWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesBySpdf.arrayResponseHoldLateBal = computeAverageFiringRateBySpdf(...
        arrayResponseLateWindowOffset, arrayOnsetHoldBal, arrayOnsetHoldBalError);
averageFiringRatesBySpdf.targetDimDelayBal = computeAverageFiringRateBySpdf(...
        targetDimDelayWindowOffset, targetDimBal);
averageFiringRatesBySpdf.targetDimDelayBalLong = computeAverageFiringRateBySpdf(...
        targetDimDelayLongWindowOffset, targetDimBal);
averageFiringRatesBySpdf.targetDimResponseBal = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimBal);
averageFiringRatesBySpdf.targetDimResponseShortHoldBal = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimShortHoldBal);
averageFiringRatesBySpdf.targetDimResponseLongHoldHal = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimLongHoldBal);
averageFiringRatesBySpdf.postTargetDimMotorResponseBal = computeAverageFiringRateBySpdf(...
        postTargetDimMotorResponseWindowOffset, targetDimBal);
averageFiringRatesBySpdf.preExitFixation = computeAverageFiringRateBySpdf(...
        preExitFixationWindowOffset, exitFixation);
averageFiringRatesBySpdf.preExitFixationEarly = computeAverageFiringRateBySpdf(...
        preExitFixationEarlyWindowOffset, exitFixation);

%% latency of response by location
cueOnset = computeResponseLatencyByLoc(cueOnset, isLocUsed);
arrayOnsetRelBal = computeResponseLatencyByLoc(arrayOnsetRelBal, isLocUsed);
arrayOnsetHoldBal = computeResponseLatencyByLoc(arrayOnsetHoldBal, isLocUsed);
targetDimBal = computeResponseLatencyByLoc(targetDimBal, isLocUsed);

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
targetDimNormalizationLogical = getTimeLogicalWithTolerance(targetDimBal.t, targetDimBal.window(1) + targetDimNormalizationWindowOffset);
% exclude motor response here

% chain time courses into long vectors for each location
allSpdfsByLoc = [cueOnset.spdfByLoc(:,cueOnsetNormalizationLogical) ...
        arrayOnset.spdfByLoc(:,arrayOnsetRelNormalizationLogical) ... % be conservative and use rel here
        arrayOnsetRelBal.spdfByLoc(:,arrayOnsetRelNormalizationLogical) ...
        arrayOnsetHoldBal.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
        targetDimBal.spdfByLoc(:,targetDimNormalizationLogical) ...
        targetDimShortHoldBal.spdfByLoc(:,targetDimNormalizationLogical) ...
        targetDimLongHoldBal.spdfByLoc(:,targetDimNormalizationLogical)];
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
targetDimInclMotorNormalizationLogical = getTimeLogicalWithTolerance(targetDimBal.t, targetDimBal.window(1) + targetDimInclMotorNormalizationWindowOffset);
exitFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(exitFixation.t, exitFixation.window(1) + exitFixationInclMotorNormalizationWindowOffset);

allSpdfsByLoc = [enterFixation.spdfByLoc(:,enterFixationInclMotorNormalizationLogical) ...
        cueOnset.spdfByLoc(:,cueOnsetInclMotorNormalizationLogical) ...
        arrayOnset.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ... % be conservative and use hold here
        arrayOnsetRelBal.spdfByLoc(:,arrayOnsetRelInclMotorNormalizationLogical) ...
        arrayOnsetHoldBal.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
        targetDimBal.spdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
        targetDimShortHoldBal.spdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
        targetDimLongHoldBal.spdfByLoc(:,targetDimInclMotorNormalizationLogical)...
        exitFixation.spdfByLoc(:,exitFixationInclMotorNormalizationLogical)];
maxFiringRateBySpdfInclMotor = max(max(allSpdfsByLoc));
minFiringRateBySpdfInclMotor = min(min(allSpdfsByLoc));

clear allSpdfsByLoc;
clear cueOnsetNormalizationLogical arrayOnsetHoldNormalizationLogical ...
        arrayOnsetRelNormalizationLogical targetDimNormalizationLogical ...
        enterFixationInclMotorNormalizationLogical cueOnsetInclMotorNormalizationLogical ...
        arrayOnsetRelInclMotorNormalizationLogical targetDimInclMotorNormalizationLogical ...
        exitFixationInclMotorNormalizationLogical;

%% compute InRF as largest mean cue response over per-condition baseline
meanCueResponseBaselineCorrByLoc = nan(nLoc, 1);
for i = 1:nLoc
    meanCueResponseBaselineCorrByLoc(i) = averageFiringRatesBySpdf.cueResponse.byLoc(i) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i);
end
[~,inRFLocByMean] = max(meanCueResponseBaselineCorrByLoc);
assert(nLoc == 4); % next line based on nLoc == 4
exRFLocByMean = mod(inRFLocByMean + 1, 4) + 1; % opposite location

% InRF and ExRF are always defined, even if there is no significant
% response

% FOR NOW, INRF is P3 or P4
[~,inRFLoc] = max(meanCueResponseBaselineCorrByLoc(3:4));
inRFLoc = inRFLoc + 2;
exRFLoc = mod(inRFLoc + 1, 4) + 1; % opposite location

%% compute InRF by max per-condition z-scored cue response
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

%% rank-sum test on x vs y using average SPIKE COUNTS
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

% versus the pre-cue baseline
cueResponseVsBaselineStatsByLoc = computeRankSumTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.cueResponse.trialRateByLoc);
    
cueTargetDelayVsBaselineStatsByLoc = computeRankSumTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc);
    
% balanced trials only
arrayResponseHoldVsBaselineStatsByLoc = computeRankSumTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.arrayResponseHoldBal.trialRateByLoc);
    
targetDimDelayVsBaselineStatsByLoc = computeRankSumTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.targetDimDelayBal.trialRateByLoc);
    
targetDimResponseVsBaselineStatsByLoc = computeRankSumTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.targetDimResponseBal.trialRateByLoc);

% all trials: hold, release, balanced, unbalanced
preExitFixationVsBaselineStatsByLoc = computeRankSumTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.preExitFixation.trialRateByLoc);

% by location -- use sign-rank test here
cueResponseVsBaselineStatsByLoc = computeSignRankTestByLoc(...
        averageFiringRatesByCount.preCueBaseline.trialRateByLoc, ...
        averageFiringRatesByCount.cueResponse.trialRateByLoc, ...
        cueResponseVsBaselineStatsByLoc); % add to existing struct
    
arrayResponseHoldVsPrevStatsByLoc = computeSignRankTestByLoc(...
        averageFiringRatesByCount.cueTargetDelayHoldBal.trialRateByLoc, ...
        averageFiringRatesByCount.arrayResponseHoldBal.trialRateByLoc);
    
targetDimResponseVsPrevStatsByLoc = computeSignRankTestByLoc(...
        averageFiringRatesByCount.targetDimDelayBal.trialRateByLoc, ...
        averageFiringRatesByCount.targetDimResponseBal.trialRateByLoc);
    
preExitFixationVsPrevStatsByLoc = computeSignRankTestByLoc(...
        averageFiringRatesByCount.preExitFixationEarly.trialRateByLoc, ...
        averageFiringRatesByCount.preExitFixation.trialRateByLoc);

%% compute cue response at InRF > baseline (SLOW)
% bootstrap on mean of baseline SPDF with 500 shuffles
% compare actual mean SPDF cue response to distribution of bootstrapped
% pre-cue baselines
% TODO using count should be faster than psth
% bootstrappedMeanPreCueBaselines = generateNullDistMeanPsth(cueOnset, preCueBaselineWindowOffset, kernelSigma, numRandomizations);
% fprintf('.');
% 
% % compute array response > cue-target delay
% bootstrappedMeanCueTargetDelayHoldBals = generateNullDistMeanPsth(arrayOnsetHoldBal, cueTargetDelayWindowOffset, kernelSigma, numRandomizations);
% fprintf('.');
% 
% % compute target-dim response > target-dim delay
% bootstrappedMeanTargetDimDelayBals = generateNullDistMeanPsth(targetDimBal, targetDimDelayWindowOffset, kernelSigma, numRandomizations);
% fprintf('.');
% 
% % compute pre-exit fixation response > earlier pre-exit fixation window
% bootstrappedMeanPreExitFixationEarlys = generateNullDistMeanPsth(exitFixation, preExitFixationEarlyWindowOffset, kernelSigma, numRandomizations);
% fprintf('.');

%% bootstrap test on x vs y using average SPDFs
% versus the pre-cue baseline
% cueResponseVsBaselineStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.cueResponse.byLoc, bootstrappedMeanPreCueBaselines, cueResponseVsBaselineStatsByLoc);
% cueTargetDelayVsBaselineStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.cueTargetDelay.byLoc, bootstrappedMeanPreCueBaselines, cueTargetDelayVsBaselineStatsByLoc);
% 
% % balanced trials only
% arrayResponseHoldVsBaselineStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.arrayResponseHoldBal.byLoc, bootstrappedMeanPreCueBaselines, arrayResponseHoldVsBaselineStatsByLoc);
% targetDimDelayVsBaselineStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.targetDimDelayBal.byLoc, bootstrappedMeanPreCueBaselines, targetDimDelayVsBaselineStatsByLoc);
% targetDimResponseVsBaselineStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.targetDimResponseBal.byLoc, bootstrappedMeanPreCueBaselines, targetDimResponseVsBaselineStatsByLoc);
% 
% % all trials: hold, release, balanced, unbalanced
% preExitFixationVsBaselineStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.preExitFixation.byLoc, bootstrappedMeanPreCueBaselines, preExitFixationVsBaselineStatsByLoc);
% 
% % versus the previous window
% arrayResponseHoldVsPrevStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.arrayResponseHoldBal.byLoc, bootstrappedMeanCueTargetDelayHoldBals, arrayResponseHoldVsPrevStatsByLoc);
% targetDimResponseVsPrevStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.targetDimResponseBal.byLoc, bootstrappedMeanTargetDimDelayBals, targetDimResponseVsPrevStatsByLoc);
% preExitFixationVsPrevStatsByLoc = getNullDistPValByLoc(averageFiringRatesBySpdf.preExitFixationEarly.byLoc, bootstrappedMeanPreExitFixationEarlys, preExitFixationVsPrevStatsByLoc);

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
cueTargetDelayAttnStats = permutationTestAttn(arrayOnset, cueTargetDelayWindowOffset, ...
        inRFLoc, exRFLoc, averageFiringRatesBySpdf.cueTargetDelay.byLoc, kernelSigma, numRandomizations);
fprintf('.');

%% rank sum test on cue-target delay period using individual trial spike counts
cueTargetDelayAttnStats.diffCount.rankSum.p = ranksum(...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{exRFLoc});
    
%% permutation test on cue-target delay period LONG
cueTargetDelayLongAttnStats = permutationTestAttn(arrayOnset, cueTargetDelayLongWindowOffset, ...
        inRFLoc, exRFLoc, averageFiringRatesBySpdf.cueTargetDelayLong.byLoc, kernelSigma, numRandomizations);
fprintf('.');

%% rank sum test on cue-target delay period LONG using individual trial spike counts
cueTargetDelayLongAttnStats.diffCount.rankSum.p = ranksum(...
        averageFiringRatesByCount.cueTargetDelayLong.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.cueTargetDelayLong.trialRateByLoc{exRFLoc});

%% permutation test on array onset response period - hold balanced trials only
arrayResponseHoldAttnStats = permutationTestAttn(arrayOnsetHoldBal, arrayResponseWindowOffset, ...
        inRFLoc, exRFLoc, averageFiringRatesBySpdf.arrayResponseHoldBal.byLoc, kernelSigma, numRandomizations);
fprintf('.');

%% rank sum test on array hold response using individual trial spike counts
arrayResponseHoldAttnStats.diffCount.rankSum.p = ranksum(...
        averageFiringRatesByCount.arrayResponseHoldBal.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.arrayResponseHoldBal.trialRateByLoc{exRFLoc});
    
%% permutation test on target-dim delay period - balanced trials only
targetDimDelayAttnStats = permutationTestAttn(targetDimBal, targetDimDelayWindowOffset, ...
        inRFLoc, exRFLoc, averageFiringRatesBySpdf.targetDimDelayBal.byLoc, kernelSigma, numRandomizations);
fprintf('.');

%% rank-sum test on target-dim delay period using individual trial spike counts
targetDimDelayAttnStats.diffCount.rankSum.p = ranksum(...
        averageFiringRatesByCount.targetDimDelayBal.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.targetDimDelayBal.trialRateByLoc{exRFLoc});

%% spatial selectivity during different periods 
% TODO address array balance issues 

% pre-cue baseline using info rate (control)
preCueBaselineInfoRate = computeInfoRatePValueByShuffle(...
        cueOnset, averageFiringRatesBySpdf.preCueBaseline, preCueBaselineWindowOffset, numRandomizations);
fprintf('.');

% cue response using info rate
cueResponseInfoRate = computeInfoRatePValueByShuffle(...
        cueOnset, averageFiringRatesBySpdf.cueResponse, cueResponseWindowOffset, numRandomizations);
fprintf('.');

% cue target delay period using info rate
cueTargetDelayInfoRate = computeInfoRatePValueByShuffle(...
        arrayOnset, averageFiringRatesBySpdf.cueTargetDelay, cueTargetDelayWindowOffset, numRandomizations);
fprintf('.');

% cue target delay period long using info rate
cueTargetDelayLongInfoRate = computeInfoRatePValueByShuffle(...
        arrayOnset, averageFiringRatesBySpdf.cueTargetDelayLong, cueTargetDelayLongWindowOffset, numRandomizations);
fprintf('.');

% array response using info rate
arrayResponseHoldInfoRate = computeInfoRatePValueByShuffle(...
        arrayOnsetHoldBal, averageFiringRatesBySpdf.arrayResponseHoldBal, arrayResponseWindowOffset, numRandomizations);
fprintf('.');

% target dim delay period using info rate
targetDimDelayInfoRate = computeInfoRatePValueByShuffle(...
        targetDimBal, averageFiringRatesBySpdf.targetDimDelayBal, targetDimDelayWindowOffset, numRandomizations);
fprintf('.');

% target dimming using info rate
targetDimResponseInfoRate = computeInfoRatePValueByShuffle(...
        targetDimBal, averageFiringRatesBySpdf.targetDimResponseBal, targetDimResponseWindowOffset, numRandomizations);
fprintf('.');

%% save
clear i;

% save ALL workspace variables not cleared above into mat file
fprintf('\n');
fprintf('Saving to %s.\n', saveFileName);
save(saveFileName);

end
