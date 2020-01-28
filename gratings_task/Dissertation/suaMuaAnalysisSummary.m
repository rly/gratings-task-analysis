
clear;
readDataLocally;
sessionInds = 1:37;

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
nsCol = [1 0 0]; %cols(6,:);
indBS = isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Broad-Spiking');
indNS = isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Narrow-Spiking');

%% cue target delay 400 ms window, AI InRF vs ExRF
goodUnits = isInPulvinar & indBS;
sub = 'BS';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelayLong(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), bsCol, isSigUnit, 0.1);
xlabel(ax, 'Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongAIDiff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = isInPulvinar & indNS;
sub = 'NS';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelayLong(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), nsCol, isSigUnit, 0.1);
xlabel(ax, 'Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongAIDiff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = isInPulvinar;
sub = 'PUL';
aiSub = attnIndices(goodUnits,1);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityCueTargetDelayLong(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), allCol, isSigUnit, 0.1);
xlabel(ax, 'Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongAIDiff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% array response hold, AI InRF vs ExRF
goodUnits = isInPulvinar & indBS;
sub = 'BS';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), bsCol, isSigUnit, 0.1);
xlabel(ax, 'Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldAIDiff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = isInPulvinar & indNS;
sub = 'NS';
aiSub = attnIndices(goodUnits,2);
isSigUnit = false(size(goodUnits));%isSignificantSelectivityArrayResponseHold(goodUnits);

ax = plotMetricDiffHistDiss(aiSub, zeros(size(aiSub)), nsCol, isSigUnit, 0.1);
xlabel(ax, 'Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldAIDiff-diss-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');




%%
cols = lines(6);
dPulCol = cols(4,:);
vPulCol = cols(5,:);

%% cue target delay 400 ms window, FF InRF vs ExRF
goodUnits = isInDPulvinar & isSignificantCueResponseInc;
sub = 'dPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), dPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = isInVPulvinar & isSignificantCueResponseInc;
sub = 'vPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), vPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% array response, AI InRF vs ExRF
goodUnits = isInDPulvinar & isSignificantCueResponseInc;
sub = 'dPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), dPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = isInVPulvinar & isSignificantCueResponseInc;
sub = 'vPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), vPulCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');



%%
cols = lines(6);
bsCol = cols(3,:);
nsCol = cols(6,:);
indBS = isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Broad-Spiking');
indNS = isInPulvinar & isSignificantCueResponseInc & strcmp(physClass, 'Narrow-Spiking');

%% cue target delay 400 ms window, FF InRF vs ExRF
goodUnits = indBS;
sub = 'dPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), bsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'vPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), nsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayLongFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% array response, AI InRF vs ExRF
goodUnits = indBS;
sub = 'dPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), bsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

goodUnits = indNS;
sub = 'vPul';

ax = computeFanoFactorAI(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
        inRFLocs(goodUnits), exRFLocs(goodUnits), nsCol, false(size(goodUnits)));
xlabel(ax, 'Fano Factor Attention Index');

plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldFFAIDiff-sfn-v%d.png', summaryDataDir, sub, v);
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

