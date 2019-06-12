
% M20170615_40b needs to be removed. 
clear;
readDataLocally;
sessionInds = 39:53;

v = 13;

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);

summaryDataDir = sprintf('%s/%s/', processedDataRootDir, 'SUA_MUA_GRATINGS_SUMMARY'); % THIS IS THE KEY DIFFERENCE
if ~exist(summaryDataDir, 'dir')
    error('No directory %s\n', summaryDataDir);
end

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

nSessions = numel(sessionInds);
nUnitsApprox = nSessions * 2; % should be equal or an underestimate

esFileNames = cell(nUnitsApprox, 1);
unitNames = cell(nUnitsApprox, 1);
unitStructs = cell(nUnitsApprox, 1);
meanWfs = cell(nUnitsApprox, 1);
physClass = cell(nUnitsApprox, 1);
isSignificantResponseVsBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
isSignificantResponseVsPreviousPeriod = false(nUnitsApprox, 4);
% isSignificantResponseVsBootstrapBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
% isSignificantResponseVsBootstrapPreviousPeriod = false(nUnitsApprox, 4);
isSignificantSelectivity = false(nUnitsApprox, 6); % 5 periods info rate
cueResponseVsBaselineDirection = zeros(nUnitsApprox, 1);
infoRates = nan(nUnitsApprox, 6); % 5 periods
diffRates = nan(nUnitsApprox, 4); % 2 delay periods + array response
attnIndices = nan(nUnitsApprox, 4); % 2 delay periods + array response
localization = cell(nUnitsApprox, 1);
isInVPulvinar = false(nUnitsApprox, 1);
isInDPulvinar = false(nUnitsApprox, 1);
earlyPreExitFixationSlope = nan(nUnitsApprox, 1);
latePreExitFixationSlope = nan(nUnitsApprox, 1);
spdfInfo = struct();

nUnitsBySessionWithDelaySelectivity = zeros(nSessions, 1);
nUnitsBySessionWithCueTargetDelaySelectivity = zeros(nSessions, 1);
nUnitsBySessionWithTargetDimDelaySelectivity = zeros(nSessions, 1);

rtFiringRateStruct = struct();

arrayHoldBalLatencyInRF = nan(nUnitsApprox, 1);
arrayHoldBalLatencyExRF = nan(nUnitsApprox, 1);
arrayRelBalLatencyInRF = nan(nUnitsApprox, 1);
arrayRelBalLatencyExRF = nan(nUnitsApprox, 1);
targetDimBalLatencyInRF = nan(nUnitsApprox, 1);
targetDimBalLatencyExRF = nan(nUnitsApprox, 1);

averageFiringRatesBySpdf = struct();
averageFiringRatesByCount = struct();
inRFLocs = nan(nUnitsApprox, 1);
exRFLocs = nan(nUnitsApprox, 1);
inRFCountNormFactor = nan(nUnitsApprox, 1);
exRFCountNormFactor = nan(nUnitsApprox, 1);

% consider sparse matrix instead
nLoc = 4;
cueTargetDelayNoiseCorrAll = cell(nSessions, 1);
arrayResponseHoldMidNoiseCorrAll = cell(nSessions, 1);
arrayResponseHoldLateNoiseCorrAll = cell(nSessions, 1);
targetDimDelayNoiseCorrAll = cell(nSessions, 1);

currentUnitIndsAll = cell(nSessions, 1);
unitCount = 0;
% should also be running a lot of shuffle tests given the number of trials

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    sessionName = recordingInfo(sessionInd).sessionName;
    saveFileName = sprintf('%s/%s-sessionInd%d-muaAnalysisSummaryData-v%d.mat', summaryDataDir, sessionName, sessionInd, v);
    fprintf('Loading file %s ...\n', saveFileName);
    S = load(saveFileName);
    
    fprintf('Found %d units...\n', numel(S.unitNames));
    % no pre-allocation here
    currentUnitInds = (unitCount + 1):(unitCount + numel(S.unitNames));
    currentUnitIndsAll{i} = currentUnitInds;
    unitCount = unitCount + numel(S.unitNames);
    
    esFileNames(currentUnitInds) = S.esFileNames;
    unitNames(currentUnitInds) = S.unitNames; % TODO include sessionInd
    unitStructs(currentUnitInds) = S.unitStructs;
    meanWfs(currentUnitInds,:) = S.meanWfs;
    physClass(currentUnitInds,:) = S.physClass;
    isSignificantResponseVsBaseline(currentUnitInds,:) = S.isSignificantResponseVsBaseline;
    isSignificantResponseVsPreviousPeriod(currentUnitInds,:) = S.isSignificantResponseVsPreviousPeriod;
%     isSignificantResponseVsBootstrapBaseline(currentUnitInds,:) = S.isSignificantResponseVsBootstrapBaseline;
%     isSignificantResponseVsBootstrapPreviousPeriod(currentUnitInds,:) = S.isSignificantResponseVsBootstrapPreviousPeriod;
    isSignificantSelectivity(currentUnitInds,:) = S.isSignificantSelectivity;
    cueResponseVsBaselineDirection(currentUnitInds,:) = S.cueResponseVsBaselineDirection;
    infoRates(currentUnitInds,:) = S.infoRates;
    diffRates(currentUnitInds,:) = S.diffRates;
    attnIndices(currentUnitInds,:) = S.attnIndices;
    localization(currentUnitInds) = S.localization;
    isInVPulvinar(currentUnitInds) = S.isInVPulvinar;
    isInDPulvinar(currentUnitInds) = S.isInDPulvinar;
    earlyPreExitFixationSlope(currentUnitInds) = S.earlyPreExitFixationSlope;
    latePreExitFixationSlope(currentUnitInds) = S.latePreExitFixationSlope;
    
    fn = fieldnames(S.spdfInfo);
    for j = 1:numel(fn)
        if isfield(spdfInfo, fn{j}) && size(S.spdfInfo.(fn{j}), 1) == numel(currentUnitInds)
            spdfInfo.(fn{j})(currentUnitInds,:) = [S.spdfInfo.(fn{j})];
        else
            spdfInfo.(fn{j}) = [S.spdfInfo.(fn{j})]; % no pre-allocation
        end
    end
    
    % test regardless of which delay period
    nUnitsBySessionWithDelaySelectivity(i) = sum(any(S.isSignificantSelectivity(:,[2 4]), 2));
    nUnitsBySessionWithCueTargetDelaySelectivity(i) = sum(S.isSignificantSelectivity(:,2));
    nUnitsBySessionWithTargetDimDelaySelectivity(i) = sum(S.isSignificantSelectivity(:,4));
    
    % overwrite each session but that's ok. they should all be the same
    enterFixationT = S.enterFixationT;
    cueOnsetT = S.cueOnsetT;
    arrayOnsetT = S.arrayOnsetT;
    targetDimT = S.targetDimT;
    exitFixationT = S.exitFixationT;
    
    % note: rtFiringRateStruct is structured differently than spdfInfo
    fn = fieldnames(S.rtFiringRateStruct);
    for j = 1:numel(fn)
        if isfield(rtFiringRateStruct, fn{j})
            rtFiringRateStruct.(fn{j})(currentUnitInds,:) = [S.rtFiringRateStruct.(fn{j})]';
        else
            rtFiringRateStruct.(fn{j}) = [S.rtFiringRateStruct.(fn{j})]'; % no pre-allocation
        end
    end
    
    arrayHoldBalLatencyInRF(currentUnitInds) = S.arrayHoldBalLatencyInRF;
    arrayHoldBalLatencyExRF(currentUnitInds) = S.arrayHoldBalLatencyExRF;
    arrayRelBalLatencyInRF(currentUnitInds) = S.arrayRelBalLatencyInRF;
    arrayRelBalLatencyExRF(currentUnitInds) = S.arrayRelBalLatencyExRF;
    targetDimBalLatencyInRF(currentUnitInds) = S.targetDimBalLatencyInRF;
    targetDimBalLatencyExRF(currentUnitInds) = S.targetDimBalLatencyExRF;
    
    fn = fieldnames(S.averageFiringRatesBySpdf);
    for j = 1:numel(fn)
        if isfield(averageFiringRatesBySpdf, fn{j})
            averageFiringRatesBySpdf.(fn{j})(currentUnitInds,:) = [S.averageFiringRatesBySpdf.(fn{j})];
        else
            averageFiringRatesBySpdf.(fn{j}) = [S.averageFiringRatesBySpdf.(fn{j})]; % no pre-allocation
        end
    end
    
    fn = fieldnames(S.averageFiringRatesByCount);
    for j = 1:numel(fn)
        if isfield(averageFiringRatesByCount, fn{j})
            averageFiringRatesByCount.(fn{j})(currentUnitInds,:) = [S.averageFiringRatesByCount.(fn{j})];
        else
            averageFiringRatesByCount.(fn{j}) = [S.averageFiringRatesByCount.(fn{j})]; % no pre-allocation
        end
    end
    inRFLocs(currentUnitInds) = S.inRFLocs;
    exRFLocs(currentUnitInds) = S.exRFLocs;
    
    inRFCountNormFactor(currentUnitInds) = S.inRFCountNormFactor;
    exRFCountNormFactor(currentUnitInds) = S.exRFCountNormFactor;    
    
    % consider sparse matrix
    cueTargetDelayNoiseCorrAll{i} = S.cueTargetDelayNoiseCorr;
    arrayResponseHoldMidNoiseCorrAll{i} = S.arrayResponseHoldMidNoiseCorr;
    arrayResponseHoldLateNoiseCorrAll{i} = S.arrayResponseHoldLateNoiseCorr;
    targetDimDelayNoiseCorrAll{i} = S.targetDimDelayNoiseCorr;
end
clear S;

% convert cell array to 2d array
meanWfs = cell2mat(meanWfs);

% initialize using exact unit count. otherwise lots of 0s
cueTargetDelayNoiseCorr = nan(unitCount, unitCount, nLoc);
arrayResponseHoldMidNoiseCorr = nan(unitCount, unitCount, nLoc);
arrayResponseHoldLateNoiseCorr = nan(unitCount, unitCount, nLoc);
targetDimDelayNoiseCorr = nan(unitCount, unitCount, nLoc);

for i = 1:nSessions
    currentUnitInds = currentUnitIndsAll{i};
    % consider sparse matrix
    cueTargetDelayNoiseCorr(currentUnitInds,currentUnitInds,:) = cueTargetDelayNoiseCorrAll{i};
    arrayResponseHoldMidNoiseCorr(currentUnitInds,currentUnitInds,:) = arrayResponseHoldMidNoiseCorrAll{i};
    arrayResponseHoldLateNoiseCorr(currentUnitInds,currentUnitInds,:) = arrayResponseHoldLateNoiseCorrAll{i};
    targetDimDelayNoiseCorr(currentUnitInds,currentUnitInds,:) = targetDimDelayNoiseCorrAll{i};
end

%% print session-wise presence of delay selectivity in MUA firing
fprintf('--------------\n');
fprintf('\n');
for i = 1:nSessions
    sessionInd = sessionInds(i);
    sessionName = recordingInfo(sessionInd).sessionName;
    if nUnitsBySessionWithCueTargetDelaySelectivity(i) > 1
        fprintf('Session %s (index %d) has %d units with cue-target delay selectivity\n', sessionName, sessionInd, nUnitsBySessionWithCueTargetDelaySelectivity(i));
    else
        warning('Session %s (index %d) has %d units with cue-target delay selectivity\n', sessionName, sessionInd, nUnitsBySessionWithCueTargetDelaySelectivity(i));
    end
end
fprintf('--------------\n');
fprintf('\n');
for i = 1:nSessions
    sessionInd = sessionInds(i);
    sessionName = recordingInfo(sessionInd).sessionName;
    if nUnitsBySessionWithTargetDimDelaySelectivity(i) > 1
        fprintf('Session %s (index %d) has %d units with target-dim delay selectivity\n', sessionName, sessionInd, nUnitsBySessionWithTargetDimDelaySelectivity(i));
    else
        warning('Session %s (index %d) has %d units with target-dim delay selectivity\n', sessionName, sessionInd, nUnitsBySessionWithTargetDimDelaySelectivity(i));
    end
end
fprintf('--------------\n');
fprintf('\n');

%% summarize
nUnitsAll = unitCount;
isCell = true(unitCount, 1); % for MUA, cannot distinguish between cell and not cell

% TODO compute these per unit above so that the right alpha is used
isSignificantAnyTaskMod = isCell & any(isSignificantResponseVsBaseline, 2);
isSignificantCueResponse = isCell & isSignificantResponseVsPreviousPeriod(:,1);
isSignificantArrayResponse = isCell & isSignificantResponseVsPreviousPeriod(:,2);
isSignificantTargetDimResponse = isCell & isSignificantResponseVsPreviousPeriod(:,3);
isSignificantVisualResponse = isCell & (isSignificantResponseVsPreviousPeriod(:,1) | ...
        isSignificantResponseVsPreviousPeriod(:,2) | isSignificantResponseVsPreviousPeriod(:,3));
isSignificantPreExitFixation = isCell & isSignificantResponseVsPreviousPeriod(:,4);

isSignificantCueResponseInc = isCell & isSignificantCueResponse & cueResponseVsBaselineDirection == 1;
isSignificantCueResponseDec = isCell & isSignificantCueResponse & cueResponseVsBaselineDirection == -1;

isCTDelayBelowBaseline = ([averageFiringRatesBySpdf.cueTargetDelay.all] < [averageFiringRatesBySpdf.preCueBaseline.all])';
isSignificantCTDelayBelowBaseline = isCell & isSignificantResponseVsBaseline(:,2) & isCTDelayBelowBaseline;

isSignificantAnySpatialSelectivity = isCell & any(isSignificantSelectivity, 2);
isSignificantEvokedSelectivity = isCell & any(isSignificantSelectivity(:,[1 3 5]), 2);
isSignificantDelaySelectivity = isCell & any(isSignificantSelectivity(:,[2 4]), 2);
isSignificantSelectivityCueTargetDelay = isCell & isSignificantSelectivity(:,2);
isSignificantSelectivityArrayResponseHold = isCell & isSignificantSelectivity(:,3);
isSignificantSelectivityTargetDimDelay = isCell & isSignificantSelectivity(:,4);
isSignificantSelectivityTargetDimResponse = isCell & isSignificantSelectivity(:,5);
isSignificantSelectivityCueTargetDelayLong = isCell & isSignificantSelectivity(:,6);

isSignificantSelectivityCueTargetDelayInc = isCell & isSignificantSelectivity(:,2) & diffRates(:,1) > 0;
isSignificantSelectivityCueTargetDelayDec = isCell & isSignificantSelectivity(:,2) & diffRates(:,1) < 0;
isSignificantSelectivityArrayHoldResponseInc = isCell & isSignificantSelectivity(:,3) & diffRates(:,2) > 0;
isSignificantSelectivityArrayHoldResponseDec = isCell & isSignificantSelectivity(:,3) & diffRates(:,2) < 0;
isSignificantSelectivityTargetDimDelayInc = isCell & isSignificantSelectivity(:,4) & diffRates(:,3) > 0;
isSignificantSelectivityTargetDimDelayDec = isCell & isSignificantSelectivity(:,4) & diffRates(:,3) < 0;
isSignificantSelectivityCueTargetDelayLongInc = isCell & isSignificantSelectivity(:,6) & diffRates(:,4) > 0;
isSignificantSelectivityCueTargetDelayLongDec = isCell & isSignificantSelectivity(:,6) & diffRates(:,4) < 0;

isInRFP1 = isCell & inRFLocs == 1;
isInRFP3 = isCell & inRFLocs == 3;

isInPulvinar = strcmp(localization, 'vPul') | strcmp(localization, 'dPul');
% strcmp(localization, 'PLd') | strcmp(localization, 'PLv') | strcmp(localization, 'PM') | strcmp(localization, 'PI');

fprintf('-------------------------------\n');
fprintf('%d/%d = %d%% units were localized to the pulvinar.\n', ...
        sum(isCell & isInPulvinar), nUnitsAll, ...
        round(sum(isCell & isInPulvinar)/nUnitsAll * 100));
fprintf('\n');

fprintf('%d/%d = %d%% units show significant task modulation (6 periods) compared to baseline (not corrected for multiple comparisons).\n', ...
        sum(isSignificantAnyTaskMod), nUnitsAll, ...
        round(sum(isSignificantAnyTaskMod)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant cue response compared to baseline.\n', ...
        sum(isSignificantCueResponse), nUnitsAll, ...
        round(sum(isSignificantCueResponse)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant array response compared to cue-target delay.\n', ...
        sum(isSignificantArrayResponse), nUnitsAll, ...
        round(sum(isSignificantArrayResponse)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant target dim response compared to target-dim delay.\n', ...
        sum(isSignificantTargetDimResponse), nUnitsAll, ...
        round(sum(isSignificantTargetDimResponse)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant pre-saccadic activity compared to previous period.\n', ...
        sum(isSignificantPreExitFixation), nUnitsAll, ...
        round(sum(isSignificantPreExitFixation)/nUnitsAll * 100));
fprintf('\n');

fprintf('%d/%d = %d%% units show significant spatial selectivity during some task period.\n', ...
        sum(isSignificantAnySpatialSelectivity), nUnitsAll, ...
        round(sum(isSignificantAnySpatialSelectivity)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant spatial selectivity in a task period involving a visual change.\n', ...
        sum(isSignificantEvokedSelectivity), nUnitsAll, ...
        round(sum(isSignificantEvokedSelectivity)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant spatial selectivity in a task period involving sustained attention.\n', ...
        sum(isSignificantDelaySelectivity), nUnitsAll, ...
        round(sum(isSignificantDelaySelectivity)/nUnitsAll * 100));
fprintf('\n');

fprintf('%d/%d = %d%% pulvinar units show significant task modulation compared to baseline.\n', ...
        sum(isSignificantAnyTaskMod & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantAnyTaskMod & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show significant cue response compared to baseline.\n', ...
        sum(isSignificantCueResponse & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantCueResponse & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show significant array response compared to cue-target delay.\n', ...
        sum(isSignificantArrayResponse & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantArrayResponse & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show significant target dim response compared to target-dim delay.\n', ...
        sum(isSignificantTargetDimResponse & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantTargetDimResponse & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show significant visual response in any period compared to the previous baseline period.\n', ...
        sum(isSignificantVisualResponse & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantVisualResponse & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show significant pre-saccadic activity compared to previous period.\n', ...
        sum(isSignificantPreExitFixation & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantPreExitFixation & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show significant cue response vs baseline, \n\t%d (%d%%) are increases, %d (%d%%) are decreases\n', ...
        sum(isSignificantCueResponse & isInPulvinar), ...
        sum(isSignificantCueResponseInc & isInPulvinar), ...
        round(sum(isSignificantCueResponseInc & isInPulvinar)/sum(isSignificantCueResponse & isInPulvinar) * 100), ...
        sum(isSignificantCueResponseDec & isInPulvinar), ...
        round(sum(isSignificantCueResponseDec & isInPulvinar)/sum(isSignificantCueResponse & isInPulvinar) * 100));
fprintf('Note: the increase / decrease labeling is based on the location with the largest significant response. There may be both significant increases and decreases.\n');
fprintf('\n');

fprintf('%d/%d = %d%% pulvinar units show significant spatial selectivity during some task period.\n', ...
        sum(isSignificantAnySpatialSelectivity & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantAnySpatialSelectivity & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show spatial selectivity in a task period involving a visual change.\n', ...
        sum(isSignificantEvokedSelectivity & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantEvokedSelectivity & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show spatial selectivity in a task period involving sustained attention.\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantDelaySelectivity & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show spatial selectivity during the cue-target delay, \n\t%d (%d%%) are InRF > ExRF, %d (%d%%) are InRF < ExRF\n', ...
        sum(isSignificantSelectivityCueTargetDelay & isInPulvinar), ...
        sum(isSignificantSelectivityCueTargetDelayInc & isInPulvinar), ...
        round(sum(isSignificantSelectivityCueTargetDelayInc & isInPulvinar)/sum(isSignificantSelectivityCueTargetDelay & isInPulvinar) * 100), ...
        sum(isSignificantSelectivityCueTargetDelayDec & isInPulvinar), ...
        round(sum(isSignificantSelectivityCueTargetDelayDec & isInPulvinar)/sum(isSignificantSelectivityCueTargetDelay & isInPulvinar) * 100));
fprintf('Of the %d units in the pulvinar that show spatial selectivity during the array hold response, \n\t%d (%d%%) are InRF > ExRF, %d (%d%%) are InRF < ExRF\n', ...
        sum(isSignificantSelectivityArrayResponseHold & isInPulvinar), ...
        sum(isSignificantSelectivityArrayHoldResponseInc & isInPulvinar), ...
        round(sum(isSignificantSelectivityArrayHoldResponseInc & isInPulvinar)/sum(isSignificantSelectivityArrayResponseHold & isInPulvinar) * 100), ...
        sum(isSignificantSelectivityArrayHoldResponseDec & isInPulvinar), ...
        round(sum(isSignificantSelectivityArrayHoldResponseDec & isInPulvinar)/sum(isSignificantSelectivityArrayResponseHold & isInPulvinar) * 100));
fprintf('Of the %d units in the pulvinar that show spatial selectivity during the target-dim delay, \n\t%d (%d%%) are InRF > ExRF, %d (%d%%) are InRF < ExRF\n', ...
        sum(isSignificantSelectivityTargetDimDelay & isInPulvinar), ...
        sum(isSignificantSelectivityTargetDimDelayInc & isInPulvinar), ...
        round(sum(isSignificantSelectivityTargetDimDelayInc & isInPulvinar)/sum(isSignificantSelectivityTargetDimDelay & isInPulvinar) * 100), ...
        sum(isSignificantSelectivityTargetDimDelayDec & isInPulvinar), ...
        round(sum(isSignificantSelectivityTargetDimDelayDec & isInPulvinar)/sum(isSignificantSelectivityTargetDimDelay & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show significant cue response vs baseline, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantCueResponse & isInPulvinar), ...
        sum(isSignificantCueResponse & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantCueResponse & strcmp(localization, 'vPul'))/sum(isSignificantCueResponse & isInPulvinar) * 100), ...
        sum(isSignificantCueResponse & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantCueResponse & strcmp(localization, 'dPul'))/sum(isSignificantCueResponse & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show significant pre-saccadic activity vs baseline, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantPreExitFixation & isInPulvinar), ...
        sum(isSignificantPreExitFixation & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantPreExitFixation & strcmp(localization, 'vPul'))/sum(isSignificantPreExitFixation & isInPulvinar) * 100), ...
        sum(isSignificantPreExitFixation & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantPreExitFixation & strcmp(localization, 'dPul'))/sum(isSignificantPreExitFixation & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show spatial selectivity during a visual change, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantEvokedSelectivity & isInPulvinar), ...
        sum(isSignificantEvokedSelectivity & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantEvokedSelectivity & strcmp(localization, 'vPul'))/sum(isSignificantEvokedSelectivity & isInPulvinar) * 100), ...
        sum(isSignificantEvokedSelectivity & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantEvokedSelectivity & strcmp(localization, 'dPul'))/sum(isSignificantEvokedSelectivity & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show spatial selectivity during attention, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100));
fprintf('\n');
fprintf('\n');

fprintf('Significant cue response compared to baseline:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantCueResponse & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantCueResponse & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantCueResponse & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantCueResponse & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant pre-saccadic activity compared to baseline:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantPreExitFixation & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantPreExitFixation & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantPreExitFixation & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantPreExitFixation & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in cue-target delay:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantSelectivityCueTargetDelay & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantSelectivityCueTargetDelay & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantSelectivityCueTargetDelay & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantSelectivityCueTargetDelay & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in array hold response:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantSelectivityArrayResponseHold & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantSelectivityArrayResponseHold & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantSelectivityArrayResponseHold & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantSelectivityArrayResponseHold & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in target-dim delay:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in target-dim response:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantSelectivityTargetDimResponse & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantSelectivityTargetDimResponse & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantSelectivityTargetDimResponse & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantSelectivityTargetDimResponse & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in both delay periods:\n');
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantSelectivityCueTargetDelay & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantSelectivityCueTargetDelay & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantSelectivityCueTargetDelay & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantSelectivityCueTargetDelay & isSignificantSelectivityTargetDimDelay & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');
fprintf('\n');

fprintf('-----------------------------\n');
precondition = isInPulvinar & isSignificantCueResponseInc;
fprintf('Of the %d units in the pulvinar that show significantly increased cue response compared to baseline:\n', sum(precondition));
fprintf('\t%d (%d%%) show significant pre-saccadic activity compared to baseline\n', sum(precondition & isSignificantPreExitFixation), ...
        round(sum(precondition & isSignificantPreExitFixation)/sum(precondition) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the cue-target delay\n', sum(precondition & isSignificantSelectivityCueTargetDelay), ...
        round(sum(precondition & isSignificantSelectivityCueTargetDelay)/sum(precondition) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the array hold response\n', sum(precondition & isSignificantSelectivityArrayResponseHold), ...
        round(sum(precondition & isSignificantSelectivityArrayResponseHold)/sum(precondition) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the target-dim delay\n', sum(precondition & isSignificantSelectivityTargetDimDelay), ...
        round(sum(precondition & isSignificantSelectivityTargetDimDelay)/sum(precondition) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the target-dim response\n', sum(precondition & isSignificantSelectivityTargetDimResponse), ...
        round(sum(precondition & isSignificantSelectivityTargetDimResponse)/sum(precondition) * 100));
fprintf('\n');

fprintf('-----------------------------\n');
precondition = isInPulvinar & isSignificantCueResponseInc;
fprintf('Of the %d units in the pulvinar that show significantly increased cue response compared to baseline:\n', sum(precondition));
fprintf('\t%d (%d%%) has InRF P1\n', sum(precondition & isInRFP1), ...
        round(sum(precondition & isInRFP1)/sum(precondition) * 100));
fprintf('\t%d (%d%%) has InRF P3\n', sum(precondition & isInRFP3), ...
        round(sum(precondition & isInRFP3)/sum(precondition) * 100));
precondition = isInPulvinar & isSignificantCueResponseDec;
fprintf('Of the %d units in the pulvinar that show significantly decreased cue response compared to baseline:\n', sum(precondition));
fprintf('\t%d (%d%%) has InRF P1 (significantly suppressed response to P3 for the 2 loc sessions)\n', sum(precondition & isInRFP1), ...
        round(sum(precondition & isInRFP1)/sum(precondition) * 100));
fprintf('\t%d (%d%%) has InRF P3 (significantly suppressed response to P1 for the 2 loc sessions)\n', sum(precondition & isInRFP3), ...
        round(sum(precondition & isInRFP3)/sum(precondition) * 100));

cols = lines(6);

unitNamesDPul = unitNames(isInDPulvinar);
unitNamesVPul = unitNames(isInVPulvinar);

save(sprintf('%s/unitNamesPul-v%d.mat', summaryDataDir, v), 'unitNamesDPul', 'unitNamesVPul');

%%
cols = lines(6);
allCol = [153 51 255]/255;
bsCol = [0 176 80]/255; %cols(3,:);
ambCol = [0.3 0.3 0.3];
nsCol = [1 0 0]; %cols(6,:);
minBaselineFR = 1.5;
troughToPeakTimes = cellfun(@(x) x.troughToPeakTimeFine * 1000, unitStructs);
indBS = isSignificantCueResponseInc & strcmp(physClass, 'Broad-Spiking') & troughToPeakTimes > 0.4001 & [averageFiringRatesBySpdf.preCueBaseline.all]' >= minBaselineFR;
indAmb = isSignificantCueResponseInc & strcmp(physClass, 'Broad-Spiking') & troughToPeakTimes > 0.3751 & troughToPeakTimes <= 0.4001 & [averageFiringRatesBySpdf.preCueBaseline.all]' >= minBaselineFR;
indNS = isSignificantCueResponseInc & (strcmp(physClass, 'Narrow-Spiking') | strcmp(physClass, 'Broad-Spiking')) & troughToPeakTimes <= 0.3751 & [averageFiringRatesBySpdf.preCueBaseline.all]' >= minBaselineFR;
indEither = (indBS | indNS | indAmb);
dPulCol = cols(6,:);
vPulCol = cols(5,:);
indDPul = isInDPulvinar & isSignificantCueResponseInc & indEither;
indVPul = isInVPulvinar & isSignificantCueResponseInc & indEither;

%% mean and image plots per-condition baseline-corrected normalized
fprintf('\n');
fprintf('Plotting normalized mean SPDFs...\n');
subdivisionNames = {'PulCueInc'};
for j = 1:numel(subdivisionNames)
    subdivisionName = subdivisionNames{j};
    yBounds = [-0.15 0.32];
    isShowLabels = 1;
    if strcmp(subdivisionName, 'PulCueInc')
        isInSubdivision = indEither;
        yBounds = [-0.14 0.3];
    elseif strcmp(subdivisionName, 'BSCueInc')
        isInSubdivision = indBS;
    elseif strcmp(subdivisionName, 'NSCueInc')
        isInSubdivision = indNS;
    elseif strcmp(subdivisionName, 'DPulCueInc')
        isInSubdivision = indDPul;
    elseif strcmp(subdivisionName, 'VPulCueInc')
        isInSubdivision = indVPul;
    else
        isInSubdivision = strcmp(localization, subdivisionName);
    end
    
    plotFileName = sprintf('%s/allSessions-%s-meanSpdfs4-diss-v%d.png', summaryDataDir, subdivisionName, v);
    quickSpdfAllEventsRowPopMeanRunner(subdivisionName, isInSubdivision, spdfInfo, ...
            cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
            yBounds, isShowLabels, plotFileName);
    
    diffPlotFileName = sprintf('%s/allSessions-%s-imagePopDiffSpdfs4-diss-v%d.png', summaryDataDir, subdivisionName, v);
    inRFPlotFileName = sprintf('%s/allSessions-%s-imagePopInRFSpdfs4-diss-v%d.png', summaryDataDir, subdivisionName, v);
    quickImagePlotAllEventsRowRunner(subdivisionName, isInSubdivision, spdfInfo, ...
            cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
            diffPlotFileName, inRFPlotFileName);
end


%% cue target delay, AI_fr
fprintf('\nCue target delay -- Firing Rate Attention Modulation Index\n');

goodUnits = indBS;
sub = 'BSCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), bsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indAmb;
sub = 'AmbCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), ambCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'NSCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), nsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

% TODO refine dPul and vPul based on distance from boundary

goodUnits = indDPul;
sub = 'DPulCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), dPulCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indVPul;
sub = 'VPulCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), vPulCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indEither;
sub = 'PULCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), allCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%
goodUnits = indNS;
sub = 'NSCueInc';

unitNamesSub = unitNames(goodUnits);
cueOnsetSpdfInRFNormSub = (spdfInfo.cueOnsetSpdfInRFNorm(goodUnits,:));
cueOnsetSpdfExRFNormSub = (spdfInfo.cueOnsetSpdfExRFNorm(goodUnits,:));
cueOnsetSpdfInRFNormErrSub = (spdfInfo.cueOnsetSpdfInRFNormErr(goodUnits,:));
cueOnsetSpdfExRFNormErrSub = (spdfInfo.cueOnsetSpdfExRFNormErr(goodUnits,:));
arrayOnsetHoldSpdfInRFNormSub = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(goodUnits,:));
arrayOnsetHoldSpdfExRFNormSub = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(goodUnits,:));
arrayOnsetHoldSpdfInRFNormErrSub = (spdfInfo.arrayOnsetHoldSpdfInRFNormErr(goodUnits,:));
arrayOnsetHoldSpdfExRFNormErrSub = (spdfInfo.arrayOnsetHoldSpdfExRFNormErr(goodUnits,:));

titleBase = sprintf('%s: Cue Onset', sub);
plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf2-cueOnset-v%d', summaryDataDir, sub, v);
makeTinyPlotsOfPopulation(cueOnsetSpdfInRFNormSub, cueOnsetSpdfInRFNormErrSub, ...
        cueOnsetSpdfExRFNormSub, cueOnsetSpdfExRFNormErrSub, arrayOnsetT, unitNamesSub, titleBase, plotFileBaseName);

titleBase = sprintf('%s: Array Onset', sub);
plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf3-arrayOnsetHold-v%d', summaryDataDir, sub, v);
makeTinyPlotsOfPopulation(arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfInRFNormErrSub, ...
        arrayOnsetHoldSpdfExRFNormSub, arrayOnsetHoldSpdfExRFNormErrSub, arrayOnsetT, unitNamesSub, titleBase, plotFileBaseName);

%% array response hold, AI_fr
fprintf('\nArray response hold -- Firing Rate Attention Modulation Index\n');

goodUnits = indBS;
sub = 'BSCueInc';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), bsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indAmb;
sub = 'AmbCueInc';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), ambCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'NSCueInc';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), nsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indDPul;
sub = 'DPulCueInc';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), dPulCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indVPul;
sub = 'VPulCueInc';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), vPulCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indEither;
sub = 'PULCueInc';
aiSub = attnIndices(goodUnits,2);
isSigUnit = isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), allCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% target-dim delay, AI_fr
fprintf('\nTarget dim delay -- Firing Rate Attention Modulation Index\n');

goodUnits = indBS;
sub = 'BSCueInc';
aiSub = attnIndices(goodUnits,3);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), bsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indAmb;
sub = 'AmbCueInc';
aiSub = attnIndices(goodUnits,3);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), ambCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'NSCueInc';
aiSub = attnIndices(goodUnits,3);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), nsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indDPul;
sub = 'DPulCueInc';
aiSub = attnIndices(goodUnits,3);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), dPulCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indVPul;
sub = 'VPulCueInc';
aiSub = attnIndices(goodUnits,3);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), vPulCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indEither;
sub = 'PULCueInc';
aiSub = attnIndices(goodUnits,3);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), allCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIfr-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%

fprintf('-----------------------------\n');
precondition = isInPulvinar & isSignificantCueResponseInc & indEither;
fprintf('Of the %d units in the pulvinar that show significantly increased cue response compared to baseline:\n', sum(precondition));
fprintf('\t%d (%d%%) show significant pre-saccadic activity compared to baseline\n', sum(precondition & isSignificantPreExitFixation), ...
        round(sum(precondition & isSignificantPreExitFixation)/sum(precondition) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the cue-target delay\n', ...
        sum(precondition & isSignificantSelectivityCueTargetDelay), ...
        round(sum(precondition & isSignificantSelectivityCueTargetDelay)/sum(precondition) * 100));
fprintf('\t\t%d attention-related increases\n', sum(precondition & isSignificantSelectivityCueTargetDelay & attnIndices(:,1) > 0));
fprintf('\t\t%d attention-related decreases\n', sum(precondition & isSignificantSelectivityCueTargetDelay & attnIndices(:,1) < 0));
fprintf('\t%d (%d%%) show significant selectivity during the array hold response\n', ...
        sum(precondition & isSignificantSelectivityArrayResponseHold), ...
        round(sum(precondition & isSignificantSelectivityArrayResponseHold)/sum(precondition) * 100));
fprintf('\t\t%d attention-related increases\n', sum(precondition & isSignificantSelectivityArrayResponseHold & attnIndices(:,2) > 0));
fprintf('\t\t%d attention-related decreases\n', sum(precondition & isSignificantSelectivityArrayResponseHold & attnIndices(:,2) < 0));
fprintf('\t%d (%d%%) show significant selectivity during the target-dim delay\n', ...
        sum(precondition & isSignificantSelectivityTargetDimDelay), ...
        round(sum(precondition & isSignificantSelectivityTargetDimDelay)/sum(precondition) * 100));
fprintf('\t\t%d attention-related increases\n', sum(precondition & isSignificantSelectivityTargetDimDelay & attnIndices(:,3) > 0));
fprintf('\t\t%d attention-related decreases\n', sum(precondition & isSignificantSelectivityTargetDimDelay & attnIndices(:,3) < 0));
fprintf('\t%d (%d%%) show significant selectivity during the target-dim response\n', ...
        sum(precondition & isSignificantSelectivityTargetDimResponse), ...
        round(sum(precondition & isSignificantSelectivityTargetDimResponse)/sum(precondition) * 100));
fprintf('\n');

goodUnits = indEither;
sum(isSignificantSelectivityCueTargetDelay(goodUnits) & ...
        isSignificantSelectivityArrayResponseHold(goodUnits))
sum(isSignificantSelectivityCueTargetDelay(goodUnits) & ...
        isSignificantSelectivityArrayResponseHold(goodUnits) & ...
        isSignificantSelectivityTargetDimDelay(goodUnits))

sum(attnIndices(indNS,1) < 0)
sum(indBS & isInDPulvinar)

fprintf('%d/%d (%d%%) units in the ventral pulvinar that show significantly increased cue response compared to baseline are narrow-spiking\n', ...
        sum(indNS & indVPul), sum(indVPul), round(sum(indNS & indVPul)/sum(indVPul)*100)); 
fprintf('%d/%d (%d%%) units that are narrow-spiking in the pulvinar and show significantly increased cue response compared to baseline are in the ventral pulvinar\n', ...
        sum(indNS & indVPul), sum(indNS), round(sum(indNS & indVPul)/sum(indNS)*100)); 
    
fprintf('%d/%d (%d%%) units that are narrow-spiking in the pulvinar and show significant cue response compared to baseline show an increase\n', ...
        sum(isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Narrow-Spiking')), ...
        sum(isInPulvinar & isSignificantCueResponse & strcmp(physClass, 'Narrow-Spiking')), ...
        round(sum(isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Narrow-Spiking')) / ...
        sum(isInPulvinar & isSignificantCueResponse & strcmp(physClass, 'Narrow-Spiking'))*100)); 

fprintf('%d/%d (%d%%) units that are broad-spiking in the pulvinar show significant cue response compared to baseline show an increase\n', ...
        sum(isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Broad-Spiking')), ...
        sum(isInPulvinar & isSignificantCueResponse & strcmp(physClass, 'Broad-Spiking')), ...
        round(sum(isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Broad-Spiking')) / ...
        sum(isInPulvinar & isSignificantCueResponse & strcmp(physClass, 'Broad-Spiking'))*100));

%% correlation AI_fr cue target delay vs AI_fr array response hold
goodUnits = indEither & isSignificantSelectivityCueTargetDelay & isSignificantSelectivityArrayResponseHold;
figure_tr_inch(5, 5); 
subaxis(1, 1, 1, 'MB', 0.15, 'ML', 0.17);
hold on;
plot([0 0], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
plot(attnIndices(goodUnits,1), attnIndices(goodUnits,2), '.', 'MarkerSize', 20, 'Color', 'k');
xlabel('AI Cue-Target Delay');
ylabel('AI Array Response Hold');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 16);
box off;
minVal = min([attnIndices(goodUnits,1); attnIndices(goodUnits,2)]);
maxVal = max([attnIndices(goodUnits,1); attnIndices(goodUnits,2)]);
xlim([minVal maxVal] + [-0.05 0.05]);
ylim([minVal maxVal] + [-0.05 0.05]);
[r,p] = corr(attnIndices(goodUnits,1), attnIndices(goodUnits,2))

%% correlation AI_fr cue target delay vs AI_fr array response hold
% & [averageFiringRatesBySpdf.preCueBaseline.all]' > 1; %
goodUnits = indEither & (isSignificantSelectivityCueTargetDelay | isSignificantSelectivityArrayResponseHold);
figure_tr_inch(5, 5); 
subaxis(1, 1, 1, 'MB', 0.15, 'ML', 0.17);
hold on;
plot([0 0], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
plot(attnIndices(goodUnits,1), attnIndices(goodUnits,2), '.', 'MarkerSize', 20, 'Color', 'k');
xlabel('AI Cue-Target Delay');
ylabel('AI Array Response Hold');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 16);
box off;
minVal = min([attnIndices(goodUnits,1); attnIndices(goodUnits,2)]);
maxVal = max([attnIndices(goodUnits,1); attnIndices(goodUnits,2)]);
xlim([minVal maxVal] + [-0.05 0.05]);
ylim([minVal maxVal] + [-0.05 0.05]);
[r,p] = corr(attnIndices(goodUnits,1), attnIndices(goodUnits,2), 'type', 'Pearson')
sum(goodUnits)

%% correlation AI_fr cue target delay vs AI_fr array response hold
% & [averageFiringRatesBySpdf.preCueBaseline.all]' > 1; %
goodUnits = indEither;
sub = 'PULCueInc';
figure_tr_inch(5, 5); 
subaxis(1, 1, 1, 'MB', 0.16, 'ML', 0.19, 'MT', 0.05);
hold on;
plot([0 0], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
h = scatter(attnIndices(goodUnits,1), attnIndices(goodUnits,2), 40, allCol, 'filled');
h.MarkerFaceAlpha = 0.8;
h.MarkerEdgeColor = zeros(3, 1);
isSig = isSignificantSelectivityArrayResponseHold & indEither;
h = scatter(attnIndices(isSig,1), attnIndices(isSig,2), 40, zeros(1, 3), 'filled');
h.MarkerFaceAlpha = 0.8;
h.MarkerEdgeColor = zeros(3, 1);
xlabel('AI_{FR} Cue-Target Delay');
ylabel('AI_{FR} Array Response Hold');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 16);
box off;
minVal = min([attnIndices(goodUnits,1); attnIndices(goodUnits,2)]);
maxVal = max([attnIndices(goodUnits,1); attnIndices(goodUnits,2)]);
xlim([minVal maxVal] + [-0.05 0.05]);
ylim([minVal maxVal] + [-0.05 0.05]);
[r,p] = corr(attnIndices(goodUnits,1), attnIndices(goodUnits,2), 'type', 'Spearman')
text(0.98, 0.02, ...
        {sprintf('Pearson r = %0.02f', r) ...
        sprintf('N = %d', sum(goodUnits))}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'Units', 'normalized', 'FontSize', 14); 
    
plotFileName = sprintf('%s/allSessions-%s-corr-cueTargetDelayAIFR-arrayResponseHoldAIFR-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% correlation AI_ff cue target delay vs AI_ff array response hold
goodUnits = indEither;
sub = 'PULCueInc';

firingRates = averageFiringRatesByCount.cueTargetDelay(goodUnits);
aiFFCueTargetDelay = computeFanoFactorAIDiss(firingRates, inRFLocs, exRFLocs);

firingRates = averageFiringRatesByCount.arrayResponseHoldBal(goodUnits);
aiFFArrayResponseHoldBal = computeFanoFactorAIDiss(firingRates, inRFLocs, exRFLocs);

toRemove = isnan(aiFFCueTargetDelay) | isnan(aiFFArrayResponseHoldBal);
aiFFCueTargetDelay(toRemove) = [];
aiFFArrayResponseHoldBal(toRemove) = [];

figure_tr_inch(5, 5); 
subaxis(1, 1, 1, 'MB', 0.16, 'ML', 0.19, 'MT', 0.05);
hold on;
plot([0 0], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-1 1], [-1 1], '-', 'Color', 0.3 * ones(3, 1));
h = scatter(aiFFCueTargetDelay, aiFFArrayResponseHoldBal, 40, allCol, 'filled');
h.MarkerFaceAlpha = 0.8;
h.MarkerEdgeColor = zeros(3, 1);
% isSig = isSignificantSelectivityArrayResponseHold & indEither;
% h = scatter(attnIndices(isSig,1), attnIndices(isSig,2), 40, zeros(1, 3), 'filled');
% h.MarkerFaceAlpha = 0.8;
% h.MarkerEdgeColor = zeros(3, 1);
xlabel('AI_{FF} Cue-Target Delay');
ylabel('AI_{FF} Array Response Hold');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 16);
box off;
minVal = min([aiFFCueTargetDelay; aiFFArrayResponseHoldBal]);
maxVal = max([aiFFCueTargetDelay; aiFFArrayResponseHoldBal]);
xlim([minVal maxVal] + [-0.05 0.05]);
ylim([minVal maxVal] + [-0.05 0.05]);
[r,p] = corr(aiFFCueTargetDelay, aiFFArrayResponseHoldBal, 'type', 'Spearman')
text(0.98, 0.02, ...
        {sprintf('Pearson r = %0.02f', r) ...
        sprintf('N = %d', sum(goodUnits))}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'Units', 'normalized', 'FontSize', 14); 

%% test
goodUnits = indBS & (isSignificantSelectivityCueTargetDelay | isSignificantSelectivityArrayResponseHold);
sub = 'BSCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), bsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');

goodUnits = indNS & [averageFiringRatesBySpdf.preCueBaseline.all]' > 5;
sub = 'BSCueInc';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelay(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), nsCol, isSigUnit, 0.1);
xlabel(ax, 'Firing Rate Attention Index');


%% for these units, does AI correlate with RT on release trials?

%% cue target delay, AI_ff
fprintf('\nCue target delay -- Fano Factor Attention Modulation Index\n');

goodUnits = indBS;
sub = 'BSCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), bsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'NSCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), nsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indDPul;
sub = 'DPulCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), dPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indVPul;
sub = 'VPulCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), vPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indEither;
sub = 'PULCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), allCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% array response hold, AI_ff
fprintf('\nArray response hold -- Fano Factor Attention Modulation Index\n');

goodUnits = indBS;
sub = 'BSCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), bsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'NSCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), nsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indDPul;
sub = 'DPulCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), dPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indVPul;
sub = 'VPulCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), vPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indEither;
sub = 'PULCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), allCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHold-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% target dim delay, AI_ff
fprintf('\nTarget dim delay -- Fano Factor Attention Modulation Index\n');

goodUnits = indBS;
sub = 'BSCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.targetDimDelayBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), bsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'NSCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.targetDimDelayBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), nsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indDPul;
sub = 'DPulCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.targetDimDelayBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), dPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indVPul;
sub = 'VPulCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.targetDimDelayBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), vPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indEither;
sub = 'PULCueInc';
ax = computeFanoFactorAI(averageFiringRatesByCount.targetDimDelayBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), allCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-targetDimDelay-AIff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');









%% compare pre-cue baseline firing rates
goodUnits = isInDPulvinar & isSignificantCueResponseInc;
dPulBaselineFR = [averageFiringRatesByCount.preCueBaseline(goodUnits).all];
goodUnits = isInVPulvinar & isSignificantCueResponseInc;
vPulBaselineFR = [averageFiringRatesByCount.preCueBaseline(goodUnits).all];

figure;
hold on;
histogram(dPulBaselineFR);
histogram(vPulBaselineFR);

p = ranksum(dPulBaselineFR, vPulBaselineFR)

%% compare pre-cue baseline firing rates
goodUnits = indBS & isSignificantCueResponseInc;
bsBaselineFR = [averageFiringRatesByCount.preCueBaseline(goodUnits).all];
goodUnits = indNS & isSignificantCueResponseInc;
nsBaselineFR = [averageFiringRatesByCount.preCueBaseline(goodUnits).all];

figure;
hold on;
histogram(bsBaselineFR);
histogram(nsBaselineFR);

p = ranksum(bsBaselineFR, nsBaselineFR)
median(bsBaselineFR) - median(nsBaselineFR)

%%
%% plot all waveforms
tWf = (-16:40) / 40000 * 1000; % ms

figure_tr_inch(6, 6);
plot(tWf, meanWfs');
title(sprintf('Mean Waveforms (ALL; N=%d)', size(meanWfs, 1)));
grid on;

% precondition = isInPulvinar & isSignificantCueResponseInc & (strcmp(physClass, 'Broad-Spiking') | strcmp(physClass, 'Narrow-Spiking'));
precondition = isSignificantCueResponseInc & (strcmp(physClass, 'Broad-Spiking') | strcmp(physClass, 'Narrow-Spiking'));
% figure_tr_inch(6, 6);
% plot(tWf, meanWfs(precondition,:)');
% title(sprintf('Mean Waveforms (N=%d)', size(meanWfs(precondition,:), 1)));
% grid on;
% set(gca, 'XMinorGrid', 'on');
troughToPeakTimeSub = cellfun(@(x) x.troughToPeakTimeFine * 1000, unitStructs(precondition));

meanWfSub = meanWfs(precondition,:);
physClassSub = physClass(precondition);
figure_tr_inch(8, 5);
hold on;
for i = 1:size(meanWfSub, 1)
    [mx,mi] = max(meanWfSub(i,:));
    if strcmp(physClassSub{i}, 'Narrow-Spiking')
        plot(tWf, meanWfSub(i,:)', 'r-');
        plot(troughToPeakTimeSub(i), mx, 'r+');
    elseif strcmp(physClassSub{i}, 'Broad-Spiking')
        plot(tWf, meanWfSub(i,:)', 'b-');
        plot(troughToPeakTimeSub(i), mx, 'b+');
    else
        warning('unknown phys class: %s\n', physClassSub{i});
%         plot(tWf, meanWfSub(i,:)', 'k:');
    end
end
plot([0.351 0.351], [-0.2 0.15], 'k--');
title(sprintf('Mean Waveforms (N=%d)', size(meanWfs(precondition,:), 1)));
grid on;
set(gca, 'XTick', -0.4:0.1:1);
xlabel('Time from Trough (ms)');

%%
figure_tr_inch(4, 4);
histogram(troughToPeakTimeSub, (0.2:0.01251:0.6));
xlabel('Trough to Peak Time (ms)');
ylabel('Number of Cells');

[dip,p_value] = HartigansDipSignifTest(sort(troughToPeakTimeSub), 1000)

%%
troughToPeakTimeAll = cellfun(@(x) x.troughToPeakTimeFine * 1000, unitStructs);

figure_tr_inch(5, 5);
subaxis(1, 1, 1, 'MB', 0.17);
hold on;
h1 = histogram(troughToPeakTimeAll(indBS), (0:0.0251:0.81));
h2 = histogram(troughToPeakTimeAll(indAmb), (0:0.0251:0.81));
h3 = histogram(troughToPeakTimeAll(indNS), (0:0.0251:0.81));
h1.FaceColor = [0 176 80]/255;
h2.FaceColor = [0.5 0.5 0.5];
h3.FaceColor = [1 0 0];

xlim([0 0.8]);
xlabel('Time from Minimum (ms)');
ylabel('Number of Units');
box off;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');


plotFileName = sprintf('%s/ns-amb-bs-width-histogram-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%
yBounds = [-0.1 0.07];
figure_tr_inch(4.5, 5.5);
subaxis(2, 1, 1, 'SV', 0.1, 'ML', 0.23, 'MB', 0.15, 'MT', 0.03);
hold on;
i = 42;
if strcmp(physClassSub{i}, 'Broad-Spiking')
    plot(tWf, meanWfSub(i,:)', '-', 'Color', bsCol, 'LineWidth', 4);
end

ylim(yBounds);
set(gca, 'XTick', -0.4:0.4:1);
set(gca, 'YTick', -0.1:0.05:0.1);
ylabel('Voltage (mV)');
set(gca, 'FontSize', 18);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 1);

subaxis(2, 1, 2);
hold on;
i = 2;
if strcmp(physClassSub{i}, 'Narrow-Spiking')
    plot(tWf, meanWfSub(i,:)', '-', 'Color', nsCol, 'LineWidth', 4);
end
ylim(yBounds);
set(gca, 'XTick', -0.4:0.4:1);
set(gca, 'YTick', -0.1:0.05:0.1);
xlabel('Time from Minimum (ms)');
ylabel('Voltage (mV)');
set(gca, 'FontSize', 18);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 1);

plotFileName = sprintf('%s/allSessions-example-bs-vs-ns-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


%% investigate cross-correlation between an NS unit and other units in session
subcondition = isInPulvinar & indEither;
matchInd = find(strcmp(unitNames(subcondition), 'M20170615_PUL_51a'));
assert(numel(matchInd) == 1);

sameSessionInds = find(cell2mat(strfind(unitNames(subcondition), 'M20170615_PUL_')));
sameSessionInds(sameSessionInds == matchInd) = [];

unitStructsSub = unitStructs(subcondition);
matchStruct = unitStructsSub{matchInd};
sameSessionStructs = unitStructsSub(sameSessionInds);

% inefficient way to do this
xcorrBounds = [-0.2 0.2];
xcorrFs = 200;
l = xcorrBounds(1):1/xcorrFs:xcorrBounds(2);
xc = zeros(numel(l), 1);
for j = 1:numel(sameSessionStructs)
    for i = 1:numel(matchStruct.ts)
        d = sameSessionStructs{j}.ts - matchStruct.ts(i);
        d(d < xcorrBounds(1) | d > xcorrBounds(2)) = [];
        if ~isempty(d)
            xc(round((d - xcorrBounds(1))*xcorrFs)+1) = xc(round((d - xcorrBounds(1))*xcorrFs)+1) + 1;
        end
    end
    figure; plot(l, xc);
end

% result: no clear relationships across simultaneously recorded NS cell and
% other cells. maybe something for M20170130_PUL_17a
