% function muaAnalysisSummary(processedDataRootDir, recordingInfoFileName, sessionInds)

v = 11;

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);

summaryDataDir = sprintf('%s/%s/', processedDataRootDir, 'MUA_GRATINGS_SUMMARY');
if ~exist(summaryDataDir, 'dir')
    error('No directory %s\n', summaryDataDir);
end

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

nSessions = numel(sessionInds);
nUnitsApprox = nSessions * 16; % should be equal or an underestimate

unitNames = cell(nUnitsApprox, 1);
isSignificantResponseVsBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
isSignificantResponseVsPreviousPeriod = false(nUnitsApprox, 4);
isSignificantResponseVsBootstrapBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
isSignificantResponseVsBootstrapPreviousPeriod = false(nUnitsApprox, 4);
isSignificantSelectivity = false(nUnitsApprox, 5); % 5 periods info rate
cueResponseVsBootstrapBaselineDirection = zeros(nUnitsApprox, 1);
infoRates = nan(nUnitsApprox, 5); % 5 periods
diffRates = nan(nUnitsApprox, 3); % 2 delay periods + array response
attnIndices = nan(nUnitsApprox, 3); % 2 delay periods + array response
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
meanRTRelInRF = nan(nSessions, 1); % currently unused
meanRTRelExRF = nan(nSessions, 1);
meanRTHoldInRF = nan(nSessions, 1);
meanRTHoldExRF = nan(nSessions, 1);

arrayHoldResponseLatencyInRF = nan(nUnitsApprox, 1);
arrayHoldResponseLatencyExRF = nan(nUnitsApprox, 1);
arrayResponseLatencyInRF = nan(nUnitsApprox, 1);
arrayResponseLatencyExRF = nan(nUnitsApprox, 1);
targetDimResponseLatencyInRF = nan(nUnitsApprox, 1);
targetDimResponseLatencyExRF = nan(nUnitsApprox, 1);

averageFiringRatesBySpdf = struct();
averageFiringRatesByCount = struct();
inRFLocs = nan(nUnitsApprox, 1);
exRFLocs = nan(nUnitsApprox, 1);

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
    unitCount = unitCount + numel(S.unitNames);
    
    unitNames(currentUnitInds) = S.unitNames;
    isSignificantResponseVsBaseline(currentUnitInds,:) = S.isSignificantResponseVsBaseline;
    isSignificantResponseVsPreviousPeriod(currentUnitInds,:) = S.isSignificantResponseVsPreviousPeriod;
    isSignificantResponseVsBootstrapBaseline(currentUnitInds,:) = S.isSignificantResponseVsBootstrapBaseline;
    isSignificantResponseVsBootstrapPreviousPeriod(currentUnitInds,:) = S.isSignificantResponseVsBootstrapPreviousPeriod;
    isSignificantSelectivity(currentUnitInds,:) = S.isSignificantSelectivity;
    cueResponseVsBootstrapBaselineDirection(currentUnitInds,:) = S.cueResponseVsBootstrapBaselineDirection;
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
    
    arrayHoldResponseLatencyInRF(currentUnitInds) = S.arrayHoldResponseLatencyInRF;
    arrayHoldResponseLatencyExRF(currentUnitInds) = S.arrayHoldResponseLatencyExRF;
    arrayResponseLatencyInRF(currentUnitInds) = S.arrayResponseLatencyInRF;
    arrayResponseLatencyExRF(currentUnitInds) = S.arrayResponseLatencyExRF;
    targetDimResponseLatencyInRF(currentUnitInds) = S.targetDimResponseLatencyInRF;
    targetDimResponseLatencyExRF(currentUnitInds) = S.targetDimResponseLatencyExRF;
    
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
    
    
end
clear S;

%% print session-wise presence of delay selectivity
fprintf('--------------\n');
fprintf('\n');
for i = 1:nSessions
    sessionInd = sessionInds(i);
    sessionName = recordingInfo(sessionInd).sessionName;
    if nUnitsBySessionWithDelaySelectivity(i) > 1
        fprintf('Session %s (index %d) has %d units with delay selectivity\n', sessionName, sessionInd, nUnitsBySessionWithDelaySelectivity(i));
    else
        warning('Session %s (index %d) has %d units with delay selectivity\n', sessionName, sessionInd, nUnitsBySessionWithDelaySelectivity(i));
    end
end
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
isSignificantCueResponse = isCell & isSignificantResponseVsBootstrapBaseline(:,1);
isSignificantArrayResponse = isCell & isSignificantResponseVsPreviousPeriod(:,2);
isSignificantTargetDimResponse = isCell & isSignificantResponseVsPreviousPeriod(:,3);
isSignificantPreExitFixation = isCell & isSignificantResponseVsPreviousPeriod(:,4);
isSignificantPreExitFixationVsBaseline = isCell & isSignificantResponseVsBootstrapBaseline(:,6);

isSignificantCueResponseInc = isCell & isSignificantCueResponse & cueResponseVsBootstrapBaselineDirection == 1;
isSignificantCueResponseDec = isCell & isSignificantCueResponse & cueResponseVsBootstrapBaselineDirection == -1;

isSignificantAnySpatialSelectivity = isCell & any(isSignificantSelectivity, 2);
isSignificantEvokedSelectivity = isCell & any(isSignificantSelectivity(:,[1 3 5]), 2);
isSignificantDelaySelectivity = isCell & any(isSignificantSelectivity(:,[2 4]), 2);
isSignificantSelectivityCueTargetDelay = isCell & isSignificantSelectivity(:,2);
isSignificantSelectivityArrayHoldResponse = isCell & isSignificantSelectivity(:,3);
isSignificantSelectivityTargetDimDelay = isCell & isSignificantSelectivity(:,4);

isSignificantSelectivityCueTargetDelayInc = isCell & isSignificantSelectivity(:,2) & diffRates(:,1) > 0;
isSignificantSelectivityCueTargetDelayDec = isCell & isSignificantSelectivity(:,2) & diffRates(:,1) < 0;
isSignificantSelectivityArrayHoldResponseInc = isCell & isSignificantSelectivity(:,3) & diffRates(:,2) > 0;
isSignificantSelectivityArrayHoldResponseDec = isCell & isSignificantSelectivity(:,3) & diffRates(:,2) < 0;
isSignificantSelectivityTargetDimDelayInc = isCell & isSignificantSelectivity(:,4) & diffRates(:,3) > 0;
isSignificantSelectivityTargetDimDelayDec = isCell & isSignificantSelectivity(:,4) & diffRates(:,3) < 0;

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
fprintf('%d/%d = %d%% units show significant target dim activity compared to target-dim delay.\n', ...
        sum(isSignificantTargetDimResponse), nUnitsAll, ...
        round(sum(isSignificantTargetDimResponse)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant pre-saccadic activity compared to previous period.\n', ...
        sum(isSignificantPreExitFixation), nUnitsAll, ...
        round(sum(isSignificantPreExitFixation)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant pre-saccadic activity compared to baseline.\n', ...
        sum(isSignificantPreExitFixationVsBaseline), nUnitsAll, ...
        round(sum(isSignificantPreExitFixationVsBaseline)/nUnitsAll * 100));
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
fprintf('%d/%d = %d%% pulvinar units show significant array response compared to baseline.\n', ...
        sum(isSignificantArrayResponse & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantArrayResponse & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar units show significant pre-saccadic activity compared to baseline.\n', ...
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
        sum(isSignificantSelectivityArrayHoldResponse & isInPulvinar), ...
        sum(isSignificantSelectivityArrayHoldResponseInc & isInPulvinar), ...
        round(sum(isSignificantSelectivityArrayHoldResponseInc & isInPulvinar)/sum(isSignificantSelectivityArrayHoldResponse & isInPulvinar) * 100), ...
        sum(isSignificantSelectivityArrayHoldResponseDec & isInPulvinar), ...
        round(sum(isSignificantSelectivityArrayHoldResponseDec & isInPulvinar)/sum(isSignificantSelectivityArrayHoldResponse & isInPulvinar) * 100));
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
        sum(isCell & isSignificantSelectivityArrayHoldResponse & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantSelectivityArrayHoldResponse & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantSelectivityArrayHoldResponse & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantSelectivityArrayHoldResponse & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
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
fprintf('\t%d (%d%%) show significant selectivity during the array hold response\n', sum(precondition & isSignificantSelectivityArrayHoldResponse), ...
        round(sum(precondition & isSignificantSelectivityArrayHoldResponse)/sum(precondition) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the target-dim delay\n', sum(precondition & isSignificantSelectivityTargetDimDelay), ...
        round(sum(precondition & isSignificantSelectivityTargetDimDelay)/sum(precondition) * 100));
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
    
stop
    
%%
figure_tr_inch(8, 8);
plot(earlyPreExitFixationSlope, latePreExitFixationSlope, '.', 'MarkerSize', 20);
plotFileName = sprintf('%s/allSessions-preSaccadicSlopes-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% 
preSaccadeWindowOffset = [-0.2 0];
preSaccadeWindowIndices = getTimeLogicalWithTolerance(exitFixationT, preSaccadeWindowOffset);

%%
% saccadeTimeIndex = find(preSaccadeWindowIndices, 1, 'last');
% figure_tr_inch(12, 10);
% subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
% hold on;
% % re-align y-axis to have 0 be the time of saccade
% plot(exitFixationT(preSaccadeWindowIndices), ...
%         spdfInfo.exitFixationSpdfInRFNorm(:,preSaccadeWindowIndices) - spdfInfo.exitFixationSpdfInRFNorm(:,saccadeTimeIndex));

%%
tsneVals = tsne(spdfInfo.exitFixationSpdfInRFNorm(:,preSaccadeWindowIndices));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
sh = scatter(tsneVals(:,1), tsneVals(:,2), 100);
sh.MarkerFaceAlpha = 0.9;

set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');

plotFileName = sprintf('%s/allSessions-tSNE-exitFixationSpdfInRFNorm-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca(spdfInfo.exitFixationSpdfInRFNorm(:,preSaccadeWindowIndices));
fprintf('\n');
fprintf('PCA: %d variables, %d observations\n', size(pcaScore, 2), size(pcaScore, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC2 explains %0.1f%% of the variance.\n', pcaPctExplained(2));
fprintf('\tPC3 explains %0.1f%% of the variance.\n', pcaPctExplained(3));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', sum(pcaPctExplained(1:2)));
fprintf('\tPC1 + PC2 + PC3 explain %0.1f%% of the variance.\n', sum(pcaPctExplained(1:3)));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
sh = scatter(pcaScore(:,1), pcaScore(:,2), 100);
sh.MarkerFaceAlpha = 0.9;

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
plot(exitFixationT(preSaccadeWindowIndices), pcaCoeff(:,1), 'LineWidth', 9);
plot(exitFixationT(preSaccadeWindowIndices), pcaCoeff(:,2), 'LineWidth', 5);
plot(exitFixationT(preSaccadeWindowIndices), pcaCoeff(:,3), 'LineWidth', 1);
% note the direction of the coefficient does not hold meaning

plotFileName = sprintf('%s/allSessions-PCA-exitFixationSpdfInRFNorm-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% test diff RT for top third vs bottom third firing rates in delay periods
meanRTRelInRFDiffThirdFiringRateCTDelay = rtFiringRateStruct.meanRTRelInRFTopThirdFiringRateCTDelay - rtFiringRateStruct.meanRTRelInRFBottomThirdFiringRateCTDelay;
meanRTRelExRFDiffThirdFiringRateCTDelay = rtFiringRateStruct.meanRTRelExRFTopThirdFiringRateCTDelay - rtFiringRateStruct.meanRTRelExRFBottomThirdFiringRateCTDelay;
meanRTHoldInRFDiffThirdFiringRateCTDelay = rtFiringRateStruct.meanRTHoldInRFTopThirdFiringRateCTDelay - rtFiringRateStruct.meanRTHoldInRFBottomThirdFiringRateCTDelay;
meanRTHoldExRFDiffThirdFiringRateCTDelay = rtFiringRateStruct.meanRTHoldExRFTopThirdFiringRateCTDelay - rtFiringRateStruct.meanRTHoldExRFBottomThirdFiringRateCTDelay;
meanRTHoldInRFDiffThirdFiringRateTDDelay = rtFiringRateStruct.meanRTHoldInRFTopThirdFiringRateTDDelay - rtFiringRateStruct.meanRTHoldInRFBottomThirdFiringRateTDDelay;
meanRTHoldExRFDiffThirdFiringRateTDDelay = rtFiringRateStruct.meanRTHoldExRFTopThirdFiringRateTDDelay - rtFiringRateStruct.meanRTHoldExRFBottomThirdFiringRateTDDelay;

condition = isSignificantCueResponseInc & isInPulvinar;
fprintf('N = %d\n', sum(condition));
meanRTRelInRFDiffThirdFiringRateCTDelaySub = meanRTRelInRFDiffThirdFiringRateCTDelay(condition);
meanRTRelExRFDiffThirdFiringRateCTDelaySub = meanRTRelExRFDiffThirdFiringRateCTDelay(condition);
meanRTHoldInRFDiffThirdFiringRateCTDelaySub = meanRTHoldInRFDiffThirdFiringRateCTDelay(condition);
meanRTHoldExRFDiffThirdFiringRateCTDelaySub = meanRTHoldExRFDiffThirdFiringRateCTDelay(condition);
meanRTHoldInRFDiffThirdFiringRateTDDelaySub = meanRTHoldInRFDiffThirdFiringRateTDDelay(condition);
meanRTHoldExRFDiffThirdFiringRateTDDelaySub = meanRTHoldExRFDiffThirdFiringRateTDDelay(condition);

[signrank(meanRTRelInRFDiffThirdFiringRateCTDelaySub) ...
        signrank(meanRTRelExRFDiffThirdFiringRateCTDelaySub) ...
        signrank(meanRTHoldInRFDiffThirdFiringRateCTDelaySub) ...
        signrank(meanRTHoldExRFDiffThirdFiringRateCTDelaySub) ...
        signrank(meanRTHoldInRFDiffThirdFiringRateTDDelaySub) ...% **
        signrank(meanRTHoldExRFDiffThirdFiringRateTDDelaySub)]

cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);

maxAbs = max(max(abs([meanRTRelInRFDiffThirdFiringRateCTDelaySub, ...
        meanRTRelExRFDiffThirdFiringRateCTDelaySub, ...
        meanRTHoldInRFDiffThirdFiringRateCTDelaySub, ...
        meanRTHoldExRFDiffThirdFiringRateCTDelaySub, ...
        meanRTHoldInRFDiffThirdFiringRateTDDelaySub, ...
        meanRTHoldExRFDiffThirdFiringRateTDDelaySub])));

binStep = 0.01;
xBounds = [-ceil(maxAbs / binStep) ceil(maxAbs / binStep)] * binStep;
histBinEdges = xBounds(1):binStep:xBounds(2);

figure_tr_inch(9, 6);
set(gcf, 'Color', 'w');
plotHs = nan(6, 1);
plotHs(1) = subaxis(2, 3, 1);
hold on;
hist1 = histogram(meanRTRelInRFDiffThirdFiringRateCTDelaySub, histBinEdges);
hist1.FaceColor = inRFCol;
plotHs(2) = subaxis(2, 3, 2);
hold on;
hist2 = histogram(meanRTHoldInRFDiffThirdFiringRateCTDelaySub, histBinEdges);
hist2.FaceColor = inRFCol;
plotHs(3) = subaxis(2, 3, 3);
hold on;
hist3 = histogram(meanRTHoldInRFDiffThirdFiringRateTDDelaySub, histBinEdges);
hist3.FaceColor = inRFCol;
plotHs(4) = subaxis(2, 3, 4);
hold on;
hist4 = histogram(meanRTRelExRFDiffThirdFiringRateCTDelaySub, histBinEdges);
hist4.FaceColor = exRFCol;
plotHs(5) = subaxis(2, 3, 5);
hold on;
hist5 = histogram(meanRTHoldExRFDiffThirdFiringRateCTDelaySub, histBinEdges);
hist5.FaceColor = exRFCol;
plotHs(6) = subaxis(2, 3, 6);
hold on;
hist6 = histogram(meanRTHoldExRFDiffThirdFiringRateTDDelaySub, histBinEdges);
hist6.FaceColor = exRFCol;

% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);

plotFileName = sprintf('%s/allSessions-rtVsFiringRateThirds-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

% note: interestingly, sessions 1-4 showed more significant
% difference shifted right than sessions 5-6,8-9 which showed more
% significant difference shifted left on hold trials
% both groups showed more sig diff shifted right on release trials

% TODO create model, compare how much variability in RT is explained by CT
% delay and TD delay, using partial determination analysis, which removes
% the correlation between the two. 

%% 
% median/third split may be more appropriate than a regression since 
% firing rates and the correlation of firing rate and RT may be more
% stepwise? there may also be a natural split, e.g. zero vs nonzero firing
% average correlation coefficient is not the way to do this. better to have
% different sessions in the model:
% https://stats.stackexchange.com/questions/8019/averaging-correlation-values
% or transform the r values using Fisher transform

condition = isSignificantCueResponseInc & isInPulvinar;
% isInDPulvinar;
fprintf('N = %d\n', sum(condition));
corrCoefRelInRFCTDelayRTSub = rtFiringRateStruct.spearmanCorrCoefRelInRFCTDelayRT(condition);
corrCoefRelExRFCTDelayRTSub = rtFiringRateStruct.spearmanCorrCoefRelExRFCTDelayRT(condition);
corrCoefHoldInRFCTDelayRTSub = rtFiringRateStruct.spearmanCorrCoefHoldInRFCTDelayRT(condition);
corrCoefHoldExRFCTDelayRTSub = rtFiringRateStruct.spearmanCorrCoefHoldExRFCTDelayRT(condition);
corrCoefHoldInRFTDDelayRTSub = rtFiringRateStruct.spearmanCorrCoefHoldInRFTDDelayRT(condition);
corrCoefHoldExRFTDDelayRTSub = rtFiringRateStruct.spearmanCorrCoefHoldExRFTDDelayRT(condition);

[signrank(atanh(corrCoefRelInRFCTDelayRTSub)) ...
        signrank(atanh(corrCoefRelExRFCTDelayRTSub)) ...
        signrank(atanh(corrCoefHoldInRFCTDelayRTSub)) ...
        signrank(atanh(corrCoefHoldExRFCTDelayRTSub)) ...
        signrank(atanh(corrCoefHoldInRFTDDelayRTSub)) ... % **
        signrank(atanh(corrCoefHoldExRFTDDelayRTSub))]

cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);

maxAbs = max(max(abs([corrCoefHoldInRFCTDelayRTSub, ...
        corrCoefRelInRFCTDelayRTSub, ...
        corrCoefHoldInRFTDDelayRTSub, ...
        corrCoefHoldExRFCTDelayRTSub, ...
        corrCoefRelExRFCTDelayRTSub, ...
        corrCoefHoldExRFTDDelayRTSub])));

binStep = 0.05;
xBounds = [-ceil(maxAbs / binStep) ceil(maxAbs / binStep)] * binStep;
histBinEdges = xBounds(1):binStep:xBounds(2);

figure_tr_inch(9, 6);
set(gcf, 'Color', 'w');
plotHs = nan(6, 1);
plotHs(1) = subaxis(2, 3, 1);
hold on;
hist1 = histogram(corrCoefRelInRFCTDelayRTSub, histBinEdges);
hist1.FaceColor = inRFCol;
xlim(xBounds);
plotHs(2) = subaxis(2, 3, 2);
hold on;
hist2 = histogram(corrCoefHoldInRFCTDelayRTSub, histBinEdges);
hist2.FaceColor = inRFCol;
xlim(xBounds);
plotHs(3) = subaxis(2, 3, 3);
hold on;
hist3 = histogram(corrCoefHoldInRFTDDelayRTSub, histBinEdges);
hist3.FaceColor = inRFCol;
xlim(xBounds);
plotHs(4) = subaxis(2, 3, 4);
hold on;
hist4 = histogram(corrCoefRelExRFCTDelayRTSub, histBinEdges);
hist4.FaceColor = exRFCol;
xlim(xBounds);
plotHs(5) = subaxis(2, 3, 5);
hold on;
hist5 = histogram(corrCoefHoldExRFCTDelayRTSub, histBinEdges);
hist5.FaceColor = exRFCol;
xlim(xBounds);
plotHs(6) = subaxis(2, 3, 6);
hold on;
hist6 = histogram(corrCoefHoldExRFTDDelayRTSub, histBinEdges);
hist6.FaceColor = exRFCol;
xlim(xBounds);


% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);

plotFileName = sprintf('%s/allSessions-rtVsFiringRateCorr-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% mean cue target delay InRF > ExRF
goodUnits = isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
cueTargetDelay = averageFiringRatesBySpdf.cueTargetDelay(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
cueTargetDelayFiringInRF = nan(nGoodUnits, 1);
cueTargetDelayFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    cueTargetDelayFiringInRF(i) = cueTargetDelay(i).byLoc(inRFLocsGoodUnits(i));
    cueTargetDelayFiringExRF(i) = cueTargetDelay(i).byLoc(exRFLocsGoodUnits(i));
end
cueTargetDelayFiringDiff = cueTargetDelayFiringInRF - cueTargetDelayFiringExRF;

fprintf('Mean cue-target delay activity InRF: %0.3f Hz\n', mean(cueTargetDelayFiringInRF));
fprintf('Mean cue-target delay activity ExRF: %0.3f Hz\n', mean(cueTargetDelayFiringExRF));
frBounds = [0 max([max(cueTargetDelayFiringInRF) max(cueTargetDelayFiringExRF)])];
maxAbsDiffFR = max(abs(cueTargetDelayFiringDiff));

binStep = 1;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

% further refining
% minFiringRate = 20;
% goodUnits2 = cueTargetDelayFiringInRF >= minFiringRate;
% cueTargetDelayFiringInRF = cueTargetDelayFiringInRF(goodUnits2);
% cueTargetDelayFiringExRF = cueTargetDelayFiringExRF(goodUnits2);
% cueTargetDelayFiringDiff = cueTargetDelayFiringInRF - cueTargetDelayFiringExRF;
% nGoodUnits = sum(goodUnits2);

figure_tr_inch(10, 5);
subaxis(1, 2, 1);
hold on;
plot(cueTargetDelayFiringInRF, cueTargetDelayFiringExRF, '.', 'MarkerSize', 20);
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(frBounds); 
ylim(frBounds);
xlabel('Mean Cue-Target Delay Firing InRF (Hz)');
ylabel('Mean Cue-Target Delay Firing ExRF (Hz)');
box off;

subaxis(1, 2, 2); 
hold on;
histogram(cueTargetDelayFiringDiff, histBinEdges);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Cue-Target Delay Firing InRF - ExRF (Hz)');
ylabel('Number of Units');
box off;

meanCueTargetDelayFiringDiff = mean(cueTargetDelayFiringDiff);
medianCueTargetDelayFiringDiff = median(cueTargetDelayFiringDiff);
p = signrank(cueTargetDelayFiringDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanCueTargetDelayFiringDiff, medianCueTargetDelayFiringDiff, p, nGoodUnits);

plotFileName = sprintf('%s/allSessions-cueTargetDelayDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% Fano Factor cue target delay InRF ? ExRF
% note that this uses spike counts not spdf rates
% TODO match firing rates (see Mitchell et al 2007)
goodUnits = isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
cueTargetDelay = averageFiringRatesByCount.cueTargetDelay(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
cueTargetDelayFanoFactorInRF = nan(nGoodUnits, 1);
cueTargetDelayFanoFactorExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    cueTargetDelayFanoFactorInRF(i) = cueTargetDelay(i).byLocSD(inRFLocsGoodUnits(i))^2 / cueTargetDelay(i).byLoc(inRFLocsGoodUnits(i));
    cueTargetDelayFanoFactorExRF(i) = cueTargetDelay(i).byLocSD(exRFLocsGoodUnits(i))^2 / cueTargetDelay(i).byLoc(exRFLocsGoodUnits(i));
end
cueTargetDelayFanoFactorDiff = cueTargetDelayFanoFactorInRF - cueTargetDelayFanoFactorExRF;

fprintf('Mean cue-target delay Fano Factor InRF: %0.3f Hz\n', mean(cueTargetDelayFanoFactorInRF));
fprintf('Mean cue-target delay Fano Factor ExRF: %0.3f Hz\n', mean(cueTargetDelayFanoFactorExRF));

frBounds = [0 max([max(cueTargetDelayFanoFactorInRF) max(cueTargetDelayFanoFactorExRF)])];
maxAbsDiffFR = max(abs(cueTargetDelayFanoFactorDiff));

binStep = 1;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

figure_tr_inch(10, 5);
subaxis(1, 2, 1);
hold on;
plot(cueTargetDelayFanoFactorInRF, cueTargetDelayFanoFactorExRF, '.', 'MarkerSize', 20);
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(frBounds); 
ylim(frBounds);
xlabel('Cue-Target Delay Fano Factor InRF (Hz)');
ylabel('Cue-Target Delay Fano Factor ExRF (Hz)');
box off;

subaxis(1, 2, 2); 
hold on;
histogram(cueTargetDelayFanoFactorDiff, histBinEdges);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Cue-Target Delay Fano Factor InRF - ExRF (Hz)');
ylabel('Number of Units');
box off;

meanCueTargetDelayFanoFactorDiff = mean(cueTargetDelayFanoFactorDiff);
medianCueTargetDelayFanoFactorDiff = median(cueTargetDelayFanoFactorDiff);
p = signrank(cueTargetDelayFanoFactorDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanCueTargetDelayFanoFactorDiff, medianCueTargetDelayFanoFactorDiff, p, nGoodUnits);

plotFileName = sprintf('%s/allSessions-cueTargetDelayFanoFactorDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% compare InRF vs ExRF array response latency
maxLatency = 0.125;
goodUnits = ~isnan(arrayHoldResponseLatencyInRF) & ~isnan(arrayHoldResponseLatencyExRF) & ...
        arrayHoldResponseLatencyInRF <= maxLatency & arrayHoldResponseLatencyExRF <= maxLatency & ...
        isInPulvinar & isSignificantCueResponseInc;
arrayOnsetHoldLatencyDiff = arrayHoldResponseLatencyInRF(goodUnits) - arrayHoldResponseLatencyExRF(goodUnits);
fprintf('Mean array onset latency InRF: %0.3f s (median: %0.3f s)\n', mean(arrayHoldResponseLatencyInRF(goodUnits)), median(arrayHoldResponseLatencyInRF(goodUnits)));
fprintf('Mean array onset latency ExRF: %0.3f s (median: %0.3f s)\n', mean(arrayHoldResponseLatencyExRF(goodUnits)), median(arrayHoldResponseLatencyExRF(goodUnits)));
fprintf('Number of units with latency reduction >= 10 ms: %d (%d%%)\n', ...
        sum(arrayOnsetHoldLatencyDiff <= -0.01), round(sum(arrayOnsetHoldLatencyDiff <= -0.01)/numel(arrayOnsetHoldLatencyDiff) * 100));
fprintf('Number of units with latency increase >= 10 ms: %d (%d%%)\n', ...
        sum(arrayOnsetHoldLatencyDiff >= 0.01), round(sum(arrayOnsetHoldLatencyDiff >= 0.01)/numel(arrayOnsetHoldLatencyDiff) * 100));

latBounds = [0 0.13];
maxAbsDiffFR = max(abs(arrayOnsetHoldLatencyDiff));

binStep = 0.005;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

figure_tr_inch(16, 5);
subaxis(1, 3, 1);
hold on;
plot(arrayHoldResponseLatencyInRF(goodUnits), arrayHoldResponseLatencyExRF(goodUnits), '.', 'MarkerSize', 20);
plot([0 1], [0 1], 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(latBounds); 
ylim(latBounds);
xlabel('Array Onset Latency InRF (s)');
ylabel('Array Onset Latency ExRF (s)');
box off;

subaxis(1, 3, 2); 
hold on;
histogram(arrayOnsetHoldLatencyDiff, histBinEdges);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Array Onset Latency InRF - ExRF (s)');
ylabel('Number of Units');
box off;

subaxis(1, 3, 3);
hold on;
c1 = cdfplot(arrayHoldResponseLatencyInRF(goodUnits));
c2 = cdfplot(arrayHoldResponseLatencyExRF(goodUnits));
set(c1, 'LineWidth', 2);
set(c2, 'LineWidth', 2);
xlim(latBounds); 
xlabel('Array Onset Latency (s)');
ylabel('Proportion of Units');
legend({'InRF', 'ExRF'}, 'Location', 'SouthEast');


meanArrayOnsetHoldLatencyDiff = mean(arrayOnsetHoldLatencyDiff);
medianArrayOnsetHoldLatencyDiff = median(arrayOnsetHoldLatencyDiff);
p = signrank(arrayOnsetHoldLatencyDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanArrayOnsetHoldLatencyDiff, medianArrayOnsetHoldLatencyDiff, p, sum(goodUnits));

plotFileName = sprintf('%s/allSessions-arrayHoldResponseLatencyDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% mean array hold response InRF > ExRF
goodUnits = isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
arrayHoldResponse = averageFiringRatesBySpdf.arrayHoldResponse(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
arrayHoldResponseMeanFiringInRF = nan(nGoodUnits, 1);
arrayHoldResponseMeanFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    arrayHoldResponseMeanFiringInRF(i) = arrayHoldResponse(i).byLoc(inRFLocsGoodUnits(i));
    arrayHoldResponseMeanFiringExRF(i) = arrayHoldResponse(i).byLoc(exRFLocsGoodUnits(i));
end
arrayHoldResponseMeanFiringDiff = arrayHoldResponseMeanFiringInRF - arrayHoldResponseMeanFiringExRF;

fprintf('Mean array response InRF: %0.3f Hz\n', mean(arrayHoldResponseMeanFiringInRF));
fprintf('Mean array response ExRF: %0.3f Hz\n', mean(arrayHoldResponseMeanFiringExRF));

frBounds = [0 max([max(arrayHoldResponseMeanFiringInRF) max(arrayHoldResponseMeanFiringExRF)])];
maxAbsDiffFR = max(abs(arrayHoldResponseMeanFiringDiff));

binStep = 1;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

figure_tr_inch(10, 5);
subaxis(1, 2, 1);
hold on;
plot(arrayHoldResponseMeanFiringInRF, arrayHoldResponseMeanFiringExRF, '.', 'MarkerSize', 20);
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(frBounds); 
ylim(frBounds);
xlabel('Mean Array Hold Response InRF (Hz)');
ylabel('Mean Array Hold Response ExRF (Hz)');
box off;

subaxis(1, 2, 2); 
hold on;
histogram(arrayHoldResponseMeanFiringDiff, histBinEdges);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Array Hold Response InRF - ExRF (Hz)');
ylabel('Number of Units');
box off;

meanArrayHoldResponseMeanFiringDiff = mean(arrayHoldResponseMeanFiringDiff);
medianArrayHoldResponseMeanFiringDiff = median(arrayHoldResponseMeanFiringDiff);
p = signrank(arrayHoldResponseMeanFiringDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanArrayHoldResponseMeanFiringDiff, medianArrayHoldResponseMeanDiff, p, nGoodUnits);

plotFileName = sprintf('%s/allSessions-arrayHoldResponseDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% peak array response InRF > ExRF
goodUnits = isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
arrayHoldResponse = averageFiringRatesBySpdf.arrayHoldResponse(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
arrayHoldResponsePeakFiringInRF = nan(nGoodUnits, 1);
arrayHoldResponsePeakFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    arrayHoldResponsePeakFiringInRF(i) = arrayHoldResponse(i).byLocMax(inRFLocsGoodUnits(i));
    arrayHoldResponsePeakFiringExRF(i) = arrayHoldResponse(i).byLocMax(exRFLocsGoodUnits(i));
end
arrayHoldResponsePeakFiringDiff = arrayHoldResponsePeakFiringInRF - arrayHoldResponsePeakFiringExRF;

fprintf('Peak array response InRF: %0.3f Hz\n', mean(arrayHoldResponsePeakFiringInRF));
fprintf('Peak array response ExRF: %0.3f Hz\n', mean(arrayHoldResponsePeakFiringExRF));

frBounds = [0 max([max(arrayHoldResponsePeakFiringInRF) max(arrayHoldResponsePeakFiringExRF)])];
maxAbsDiffFR = max(abs(arrayHoldResponsePeakFiringDiff));

binStep = 2;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

figure_tr_inch(10, 5);
subaxis(1, 2, 1);
hold on;
plot(arrayHoldResponsePeakFiringInRF, arrayHoldResponsePeakFiringExRF, '.', 'MarkerSize', 20);
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(frBounds); 
ylim(frBounds);
xlabel('Mean Array Hold Response InRF (Hz)');
ylabel('Mean Array Hold Response ExRF (Hz)');
box off;

subaxis(1, 2, 2); 
hold on;
histogram(arrayHoldResponsePeakFiringDiff, histBinEdges);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Array Hold Response InRF - ExRF (Hz)');
ylabel('Number of Units');
box off;

meanArrayHoldResponsePeakFiringDiff = mean(arrayHoldResponsePeakFiringDiff);
medianArrayHoldResponsePeakFiringDiff = median(arrayHoldResponsePeakFiringDiff);
p = signrank(arrayHoldResponsePeakFiringDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanArrayHoldResponsePeakFiringDiff, medianArrayHoldResponsePeakFiringDiff, p, nGoodUnits);

plotFileName = sprintf('%s/allSessions-arrayHoldResponsePeakDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% Fano factor array response InRF < ExRF
% note that this uses spike counts not spdf rates
% TODO match firing rates (see Mitchell et al 2007)
goodUnits = isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
arrayHoldResponse = averageFiringRatesByCount.arrayHoldResponse(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
arrayHoldResponseFanoFactorInRF = nan(nGoodUnits, 1);
arrayHoldResponseFanoFactorExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    arrayHoldResponseFanoFactorInRF(i) = arrayHoldResponse(i).byLocSD(inRFLocsGoodUnits(i))^2 / arrayHoldResponse(i).byLoc(inRFLocsGoodUnits(i));
    arrayHoldResponseFanoFactorExRF(i) = arrayHoldResponse(i).byLocSD(exRFLocsGoodUnits(i))^2 / arrayHoldResponse(i).byLoc(exRFLocsGoodUnits(i));
end
arrayHoldResponseFanoFactorDiff = arrayHoldResponseFanoFactorInRF - arrayHoldResponseFanoFactorExRF;

fprintf('Fano Factor array response InRF: %0.3f Hz\n', mean(arrayHoldResponseFanoFactorInRF));
fprintf('Fano Factor array response ExRF: %0.3f Hz\n', mean(arrayHoldResponseFanoFactorExRF));

frBounds = [0 max([max(arrayHoldResponseFanoFactorInRF) max(arrayHoldResponseFanoFactorExRF)])];
maxAbsDiffFR = max(abs(arrayHoldResponseFanoFactorDiff));

binStep = 1;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

figure_tr_inch(10, 5);
subaxis(1, 2, 1);
hold on;
plot(arrayHoldResponseFanoFactorInRF, arrayHoldResponseFanoFactorExRF, '.', 'MarkerSize', 20);
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(frBounds); 
ylim(frBounds);
xlabel('Array Hold Fano Factor InRF (Hz)');
ylabel('Array Hold Fano Factor ExRF (Hz)');
box off;

subaxis(1, 2, 2); 
hold on;
histogram(arrayHoldResponseFanoFactorDiff, histBinEdges);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Array Hold Fano Factor InRF - ExRF (Hz)');
ylabel('Number of Units');
box off;

meanArrayHoldResponseFanoFactorDiff = mean(arrayHoldResponseFanoFactorDiff);
medianArrayHoldResponseFanoFactorDiff = median(arrayHoldResponseFanoFactorDiff);
p = signrank(arrayHoldResponseFanoFactorDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanArrayHoldResponseFanoFactorDiff, medianArrayHoldResponseFanoFactorDiff, p, nGoodUnits);

plotFileName = sprintf('%s/allSessions-arrayHoldResponseFanoFactorDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% does attention-related difference in delay firing correlate with reduction of latency to array onset
goodUnits = ~isnan(arrayHoldResponseLatencyInRF) & ~isnan(arrayHoldResponseLatencyExRF) & ...
        arrayHoldResponseLatencyInRF <= maxLatency & arrayHoldResponseLatencyExRF <= maxLatency & ...
        isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
cueTargetDelay = averageFiringRatesBySpdf.cueTargetDelay(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
cueTargetDelayFiringInRF = nan(nGoodUnits, 1);
cueTargetDelayFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    cueTargetDelayFiringInRF(i) = cueTargetDelay(i).byLoc(inRFLocsGoodUnits(i));
    cueTargetDelayFiringExRF(i) = cueTargetDelay(i).byLoc(exRFLocsGoodUnits(i));
end
cueTargetDelayFiringDiff = cueTargetDelayFiringInRF - cueTargetDelayFiringExRF;
arrayOnsetHoldLatencyDiff = arrayHoldResponseLatencyInRF(goodUnits) - arrayHoldResponseLatencyExRF(goodUnits);

maxAbsDiffFR = max(abs(cueTargetDelayFiringDiff));
frDiffBounds = maxAbsDiffFR * [-1 1];
maxAbsDiffLat = max(abs(arrayOnsetHoldLatencyDiff));
latDiffBounds = maxAbsDiffLat * [-1 1];

figure_tr_inch(6, 6);
hold on;
plot(frDiffBounds, frDiffBounds, 'Color', 0.3*ones(3, 1)); 
plot(frDiffBounds, [0 0], 'Color', zeros(3, 1)); 
plot([0 0], frDiffBounds, 'Color', zeros(3, 1)); 
plot(cueTargetDelayFiringDiff, arrayOnsetHoldLatencyDiff, '.', 'MarkerSize', 20, 'Color', lines(1));
xlim(frDiffBounds);
ylim(latDiffBounds);
xlabel('Cue-Target Delay Firing InRF-ExRF (Hz)');
ylabel('Array Hold Response Latency InRF-ExRF (s)');

[r, p] = corr(cueTargetDelayFiringDiff, arrayOnsetHoldLatencyDiff, 'type', 'Spearman');
fprintf('Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, nGoodUnits);
    
plotFileName = sprintf('%s/allSessions-cueTargetDelayFiringDiffVsArrayHoldResponseLatencyDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
    
%% does attention-related difference in delay firing correlate with attention-related difference in array response
goodUnits = isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
cueTargetDelay = averageFiringRatesBySpdf.cueTargetDelay(goodUnits);
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
cueTargetDelayFiringInRF = nan(nGoodUnits, 1);
cueTargetDelayFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    cueTargetDelayFiringInRF(i) = cueTargetDelay(i).byLoc(inRFLocsGoodUnits(i));
    cueTargetDelayFiringExRF(i) = cueTargetDelay(i).byLoc(exRFLocsGoodUnits(i));
end
cueTargetDelayFiringDiff = cueTargetDelayFiringInRF - cueTargetDelayFiringExRF;

arrayHoldResponse = averageFiringRatesBySpdf.arrayHoldResponse(goodUnits);
arrayHoldResponseMeanFiringInRF = nan(nGoodUnits, 1);
arrayHoldResponseMeanFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    arrayHoldResponseMeanFiringInRF(i) = arrayHoldResponse(i).byLoc(inRFLocsGoodUnits(i));
    arrayHoldResponseMeanFiringExRF(i) = arrayHoldResponse(i).byLoc(exRFLocsGoodUnits(i));
end
arrayHoldResponseMeanFiringDiff = arrayHoldResponseMeanFiringInRF - arrayHoldResponseMeanFiringExRF;

maxAbsDiffFR = max([max(abs(cueTargetDelayFiringDiff)) max(abs(arrayHoldResponseMeanFiringDiff))]);
frDiffBounds = maxAbsDiffFR * [-1 1];

figure_tr_inch(6, 6);
hold on;
plot(frDiffBounds, frDiffBounds, 'Color', 0.3*ones(3, 1)); 
plot(frDiffBounds, [0 0], 'Color', zeros(3, 1)); 
plot([0 0], frDiffBounds, 'Color', zeros(3, 1)); 
plot(cueTargetDelayFiringDiff, arrayHoldResponseMeanFiringDiff, '.', 'MarkerSize', 20, 'Color', lines(1));
xlim(frDiffBounds);
ylim(frDiffBounds);
xlabel('Cue-Target Delay Firing InRF-ExRF (Hz)');
ylabel('Array Hold Response Firing InRF-ExRF (Hz)');

[r, p] = corr(cueTargetDelayFiringDiff, arrayHoldResponseMeanFiringDiff, 'type', 'Spearman');
fprintf('Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, nGoodUnits);
    
plotFileName = sprintf('%s/allSessions-cueTargetDelayFiringDiffVsArrayHoldResponseFiringDiff-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
    
%% does attention-related reduction in latency of array response correlate with attention-related difference in array response
maxLatency = 0.075;
goodUnits = ~isnan(arrayHoldResponseLatencyInRF) & ~isnan(arrayHoldResponseLatencyExRF) & ...
        arrayHoldResponseLatencyInRF <= maxLatency & arrayHoldResponseLatencyExRF <= maxLatency & ...
        isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
arrayHoldResponse = averageFiringRatesBySpdf.arrayHoldResponse(goodUnits);
arrayHoldResponseMeanFiringInRF = nan(nGoodUnits, 1);
arrayHoldResponseMeanFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    arrayHoldResponseMeanFiringInRF(i) = arrayHoldResponse(i).byLoc(inRFLocsGoodUnits(i));
    arrayHoldResponseMeanFiringExRF(i) = arrayHoldResponse(i).byLoc(exRFLocsGoodUnits(i));
end
arrayHoldResponseMeanFiringDiff = arrayHoldResponseMeanFiringInRF - arrayHoldResponseMeanFiringExRF;
arrayOnsetHoldLatencyDiff = arrayHoldResponseLatencyInRF(goodUnits) - arrayHoldResponseLatencyExRF(goodUnits);

maxAbsDiffFR = max(abs(arrayHoldResponseMeanFiringDiff));
frDiffBounds = maxAbsDiffFR * [-1 1];
maxAbsDiffLat = max(abs(arrayOnsetHoldLatencyDiff));
latDiffBounds = maxAbsDiffLat * [-1 1];

figure_tr_inch(6, 6);
hold on;
plot(frDiffBounds, frDiffBounds, 'Color', 0.3*ones(3, 1)); 
plot(frDiffBounds, [0 0], 'Color', zeros(3, 1)); 
plot([0 0], frDiffBounds, 'Color', zeros(3, 1)); 
plot(arrayHoldResponseMeanFiringDiff, arrayOnsetHoldLatencyDiff, '.', 'MarkerSize', 20, 'Color', lines(1));
xlim(frDiffBounds);
ylim(latDiffBounds);
xlabel('Array Hold Response Firing InRF-ExRF (Hz)');
ylabel('Array Hold Response Latency InRF-ExRF (s)');

[r, p] = corr(arrayHoldResponseMeanFiringDiff, arrayOnsetHoldLatencyDiff, 'type', 'Spearman');
fprintf('Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, nGoodUnits);
    
plotFileName = sprintf('%s/allSessions-arrayHoldResponseFiringDiffVsArrayHoldResponseLatencyDiff-maxLat%0.3f-v%d.png', summaryDataDir, maxLatency, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% does latency of array response correlate with strength of array response
goodUnits = ~isnan(arrayHoldResponseLatencyInRF) & ~isnan(arrayHoldResponseLatencyExRF) & ...
        arrayHoldResponseLatencyInRF <= maxLatency & arrayHoldResponseLatencyExRF <= maxLatency & ...
        isInPulvinar & isSignificantCueResponseInc;
nGoodUnits = sum(goodUnits);
% unboxing
inRFLocsGoodUnits = inRFLocs(goodUnits);
exRFLocsGoodUnits = exRFLocs(goodUnits);
arrayHoldResponse = averageFiringRatesBySpdf.arrayHoldResponse(goodUnits);
arrayHoldResponseMeanFiringInRF = nan(nGoodUnits, 1);
arrayHoldResponseMeanFiringExRF = nan(nGoodUnits, 1);
for i = 1:nGoodUnits
    arrayHoldResponseMeanFiringInRF(i) = arrayHoldResponse(i).byLoc(inRFLocsGoodUnits(i));
    arrayHoldResponseMeanFiringExRF(i) = arrayHoldResponse(i).byLoc(exRFLocsGoodUnits(i));
end
arrayHoldResponseLatencyInRFSub = arrayHoldResponseLatencyInRF(goodUnits);
arrayHoldResponseLatencyExRFSub = arrayHoldResponseLatencyExRF(goodUnits);

frBounds = [0 max([max(arrayHoldResponseMeanFiringInRF) max(arrayHoldResponseMeanFiringExRF)])];
latBounds = [0 0.125];
cols = lines(2);

figure_tr_inch(6, 6);
hold on;
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
plot(frBounds, [0 0], 'Color', zeros(3, 1)); 
plot([0 0], latBounds, 'Color', zeros(3, 1)); 
h1 = plot(arrayHoldResponseMeanFiringInRF, arrayHoldResponseLatencyInRFSub, '.', 'MarkerSize', 20, 'Color', cols(1,:));
h2 = plot(arrayHoldResponseMeanFiringExRF, arrayHoldResponseLatencyExRFSub, '.', 'MarkerSize', 20, 'Color', cols(2,:));
xlim(frBounds);
ylim(latBounds);
xlabel('Array Hold Response Firing (Hz)');
ylabel('Array Hold Response Latency (s)');
legend([h1 h2], {'InRF', 'ExRF'}, 'Location', 'NorthEast');

[r, p] = corr(arrayHoldResponseMeanFiringInRF, arrayHoldResponseLatencyInRFSub, 'type', 'Spearman');
fprintf('InRF: Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, nGoodUnits);
[r, p] = corr(arrayHoldResponseMeanFiringExRF, arrayHoldResponseLatencyExRFSub, 'type', 'Spearman');
fprintf('ExRF: Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, nGoodUnits);
    
plotFileName = sprintf('%s/allSessions-arrayHoldResponseFiringVsArrayHoldResponseLatency-v%d.png', summaryDataDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% per-condition baseline-corrected normalized mean
fprintf('\n');
fprintf('Plotting normalized mean SPDFs...\n');
% subdivisions = {'PM', 'PLd', 'PLv', 'PI', 'PUL', 'all'};
subdivisions = {'all', 'PUL', 'dPul', 'vPul'};%, 'notVPul', 'notDPul'};
for j = 1:numel(subdivisions)
    subdivision = subdivisions{j};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(nUnitsAll, 1);
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = isInPulvinar;
    elseif strcmp(subdivision, 'vPul')
        isInSubdivision = isInVPulvinar;
    elseif strcmp(subdivision, 'notVPul')
        isInSubdivision = ~isInVPulvinar;
    elseif strcmp(subdivision, 'dPul')
        isInSubdivision = isInDPulvinar;
    elseif strcmp(subdivision, 'notDPul')
        isInSubdivision = ~isInDPulvinar;
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    condition = isCell & isInSubdivision;
    
    enterFixationSpdfInRFNormSub = (spdfInfo.enterFixationSpdfInRFNorm(condition,:));
    enterFixationSpdfExRFNormSub = (spdfInfo.enterFixationSpdfExRFNorm(condition,:));
    cueOnsetSpdfInRFNormSub = (spdfInfo.cueOnsetSpdfInRFNorm(condition,:));
    cueOnsetSpdfExRFNormSub = (spdfInfo.cueOnsetSpdfExRFNorm(condition,:));
    arrayOnsetHoldSpdfInRFNormSub = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(condition,:));
    arrayOnsetHoldSpdfExRFNormSub = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(condition,:));
    targetDimSpdfInRFNormSub = (spdfInfo.targetDimSpdfInRFNorm(condition,:));
    targetDimSpdfExRFNormSub = (spdfInfo.targetDimSpdfExRFNorm(condition,:));
    exitFixationSpdfInRFNormSub = (spdfInfo.exitFixationSpdfInRFNorm(condition,:));
    exitFixationSpdfExRFNormSub = (spdfInfo.exitFixationSpdfExRFNorm(condition,:));

    fprintf('\t%s: %d cells\n', subdivision, sum(condition));

%     plotFileName = sprintf('%s/allSessions-%s-meanSpdfs3-v%d.png', summaryDataDir, subdivision, v);
%     fprintf('Saving to %s...\n', plotFileName);
%     
%     quickSpdfAllEvents3InARowPopMean(cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
%             arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
%             targetDimSpdfExRFNormSub, cueOnsetT, arrayOnsetT, targetDimT, plotFileName);

    plotFileName = sprintf('%s/allSessions-%s-meanSpdfs5-v%d.png', summaryDataDir, subdivision, v);
    fprintf('Saving to %s...\n', plotFileName);
    
    quickSpdfAllEvents5InARowPopMean(enterFixationSpdfInRFNormSub, enterFixationSpdfExRFNormSub, ...
            cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
            arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
            targetDimSpdfExRFNormSub, exitFixationSpdfInRFNormSub, exitFixationSpdfExRFNormSub, ...
            enterFixationT, cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, plotFileName);    
end

%% mega figure of tiny bc normalized plots per unit by subdivision
fprintf('\n');
fprintf('Plotting mega figure of tiny baseline-corrected, normalized mean SPDFs...\n');
% subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
subdivisions = {'PULCueInc', 'PULCueDec'};%'all', 'vPul', 'dPul'};%, 'notVPul', 'notDPul'};
for j = 1:numel(subdivisions)
    subdivision = subdivisions{j};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(nUnitsAll, 1);
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = isInPulvinar;
    elseif strcmp(subdivision, 'vPul')
        isInSubdivision = isInVPulvinar;
    elseif strcmp(subdivision, 'notVPul')
        isInSubdivision = ~isInVPulvinar;
    elseif strcmp(subdivision, 'dPul')
        isInSubdivision = isInDPulvinar;
    elseif strcmp(subdivision, 'notDPul')
        isInSubdivision = ~isInDPulvinar;
    elseif strcmp(subdivision, 'PULCueInc')
        isInSubdivision = isInPulvinar & isSignificantCueResponseInc;
    elseif strcmp(subdivision, 'PULCueDec')
        isInSubdivision = isInPulvinar & isSignificantCueResponseDec;
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    condition = isCell & isInSubdivision;
    unitNamesSub = unitNames(condition);
    
    enterFixationSpdfInRFNormSub = (spdfInfo.enterFixationSpdfInRFNorm(condition,:));
    enterFixationSpdfExRFNormSub = (spdfInfo.enterFixationSpdfExRFNorm(condition,:));
    cueOnsetSpdfInRFNormSub = (spdfInfo.cueOnsetSpdfInRFNorm(condition,:));
    cueOnsetSpdfExRFNormSub = (spdfInfo.cueOnsetSpdfExRFNorm(condition,:));
    arrayOnsetHoldSpdfInRFNormSub = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(condition,:));
    arrayOnsetHoldSpdfExRFNormSub = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(condition,:));
    targetDimSpdfInRFNormSub = (spdfInfo.targetDimSpdfInRFNorm(condition,:));
    targetDimSpdfExRFNormSub = (spdfInfo.targetDimSpdfExRFNorm(condition,:));
    exitFixationSpdfInRFNormSub = (spdfInfo.exitFixationSpdfInRFNorm(condition,:));
    exitFixationSpdfExRFNormSub = (spdfInfo.exitFixationSpdfExRFNorm(condition,:));
    
    enterFixationSpdfInRFNormErrSub = (spdfInfo.enterFixationSpdfInRFNormErr(condition,:));
    enterFixationSpdfExRFNormErrSub = (spdfInfo.enterFixationSpdfExRFNormErr(condition,:));
    cueOnsetSpdfInRFNormErrSub = (spdfInfo.cueOnsetSpdfInRFNormErr(condition,:));
    cueOnsetSpdfExRFNormErrSub = (spdfInfo.cueOnsetSpdfExRFNormErr(condition,:));
    arrayOnsetHoldSpdfInRFNormErrSub = (spdfInfo.arrayOnsetHoldSpdfInRFNormErr(condition,:));
    arrayOnsetHoldSpdfExRFNormErrSub = (spdfInfo.arrayOnsetHoldSpdfExRFNormErr(condition,:));
    targetDimSpdfInRFNormErrSub = (spdfInfo.targetDimSpdfInRFNormErr(condition,:));
    targetDimSpdfExRFNormErrSub = (spdfInfo.targetDimSpdfExRFNormErr(condition,:));
    exitFixationSpdfInRFNormErrSub = (spdfInfo.exitFixationSpdfInRFNormErr(condition,:));
    exitFixationSpdfExRFNormErrSub = (spdfInfo.exitFixationSpdfExRFNormErr(condition,:));

    fprintf('\t%s: %d cells\n', subdivision, sum(condition));
    
    titleBase = sprintf('%s Cells: Enter Fixation', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf1-enterFixation-v%d', summaryDataDir, subdivision, v);
    makeTinyPlotsOfPopulation(enterFixationSpdfInRFNormSub, enterFixationSpdfInRFNormErrSub, ...
            enterFixationSpdfExRFNormSub, enterFixationSpdfExRFNormErrSub, enterFixationT, unitNamesSub, titleBase, plotFileBaseName);
    
    titleBase = sprintf('%s Cells: Cue Onset', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf2-cueOnset-v%d', summaryDataDir, subdivision, v);
    makeTinyPlotsOfPopulation(cueOnsetSpdfInRFNormSub, cueOnsetSpdfInRFNormErrSub, ...
            cueOnsetSpdfExRFNormSub, cueOnsetSpdfExRFNormErrSub, cueOnsetT, unitNamesSub, titleBase, plotFileBaseName);
    
    titleBase = sprintf('%s Cells: Array Onset Hold', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf3-arrayOnsetHold-v%d', summaryDataDir, subdivision, v);
    makeTinyPlotsOfPopulation(arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfInRFNormErrSub, ...
            arrayOnsetHoldSpdfExRFNormSub, arrayOnsetHoldSpdfExRFNormErrSub, arrayOnsetT, unitNamesSub, titleBase, plotFileBaseName);
    
    titleBase = sprintf('%s Cells: Target Dim', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf4-targetDim-v%d', summaryDataDir, subdivision, v);
    makeTinyPlotsOfPopulation(targetDimSpdfInRFNormSub, targetDimSpdfInRFNormErrSub, ...
            targetDimSpdfExRFNormSub, targetDimSpdfExRFNormErrSub, targetDimT, unitNamesSub, titleBase, plotFileBaseName);
        
    titleBase = sprintf('%s Cells: Exit Fixation', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf5-exitFixation-v%d', summaryDataDir, subdivision, v);
    makeTinyPlotsOfPopulation(exitFixationSpdfInRFNormSub, exitFixationSpdfInRFNormErrSub, ...
            exitFixationSpdfExRFNormSub, exitFixationSpdfExRFNormErrSub, exitFixationT, unitNamesSub, titleBase, plotFileBaseName);
end

stop
%% Pie Chart localizing attentional modulation in the pulvinar
% We found x pulvinar cells that show significant attentional modulation. 
% In which subdivision of the pulvinar are they located?
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.1, 'MB', 0.1, 'MT', 0.1, 'MR', 0.1);
set(gcf, 'Color', 'white');
hold on;

g2 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g2all = isCell & strcmp(localization, 'PLd');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PLv');
g3all = isCell & strcmp(localization, 'PLv');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g4all = isCell & strcmp(localization, 'PM');
g5 = isSignificantDelaySelectivity & strcmp(localization, 'PI');
g5all = isCell & strcmp(localization, 'PI');

x = [sum(g2)/sum(g2all) sum(g3)/sum(g3all) sum(g4)/sum(g4all) sum(g5)/sum(g5all)];
h = pie(x / sum(x));
h(1).FaceColor = [75 172 198]/255; % use colors from the ppt
h(3).FaceColor = [255 192 0]/255;
h(5).FaceColor = [255 0 0]/255;
h(7).FaceColor = [100 50 190]/255;
h(2).String = sprintf('PLd: %s', h(2).String);
h(4).String = sprintf('PLv: %s', h(4).String);
h(6).String = sprintf('PM: %s', h(6).String);
h(8).String = sprintf('PI: %s', h(8).String);
set(h(2), 'FontSize', 20);
set(h(4), 'FontSize', 20);
set(h(6), 'FontSize', 20);
set(h(8), 'FontSize', 20);
axis off;
% adjust the locations of the text by hand

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by Sig in Which 
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity;
g2 = isCell & isInPulvinar; % & any(isSignificantCueTargetDelay, 2) & ~any(isSignificantTargetDimDelay, 2);
% g3 = isCell & pulLocalization; % & ~any(isSignificantCueTargetDelay, 2) & any(isSignificantTargetDimDelay, 2);
% g4 = isCell & pulLocalization; % & any(isSignificantCueTargetDelay, 2) & any(isSignificantTargetDimDelay, 2);


plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
% plot(infoRates(g3, 2), infoRates(g3, 4), ...
%         '.', 'MarkerSize', 20, 'Color', [1 0.5 1]);
% plot(infoRates(g4, 2), infoRates(g4, 4), ...
%         '.', 'MarkerSize', 20, 'Color', [0.5 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = max([currXLim(2) currYLim(2)]);
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by PLd vs PM vs PI
figure_tr_inch(6.5, 6.5);
subaxis(1, 1, 1, 'ML', 0.22, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PLv');
g5 = isSignificantDelaySelectivity & strcmp(localization, 'PI');

cols = lines(5);
cols = [cols(1,:); cols(3,:); cols(5,:); cols(2,:)];

maxLim = 0.55;
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

% plot(infoRates(g1,2), infoRates(g1,4), ...
%         '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
% plot(infoRates(g2, 2), infoRates(g2, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(1,:));
% plot(infoRates(g3, 2), infoRates(g3, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(2,:));
% plot(infoRates(g4, 2), infoRates(g4, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(3,:));
% plot(infoRates(g5, 2), infoRates(g5, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(4,:));
sh = scatter(infoRates(g1,2), infoRates(g1,4), 100, 0.85*ones(1, 3), 'MarkerFaceColor', 0.85*ones(3, 1));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g2,2), infoRates(g2,4), 100, cols(1,:), 'MarkerFaceColor', cols(1,:));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g3,2), infoRates(g3,4), 100, cols(2,:), 'MarkerFaceColor', cols(2,:));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g4,2), infoRates(g4,4), 100, cols(3,:), 'MarkerFaceColor', cols(3,:));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g5,2), infoRates(g5,4), 100, cols(4,:), 'MarkerFaceColor', cols(4,:));
sh.MarkerFaceAlpha = 0.9;

currXLim = xlim();
currYLim = ylim();
%max([currXLim(2) currYLim(2)]);
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
% xlabel({'Selectivity in', 'Cue-Target Delay (bits/s)'});
% ylabel({'Selectivity in', 'Target-Dim Delay (bits/s)'});
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-PLdvsPMvsPI-v%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by PLd vs PLv vs PM vs PI
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity;% & isInPulvinar;
g2 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PLv');
g5 = isSignificantDelaySelectivity & strcmp(localization, 'PI');

plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(infoRates(g3, 2), infoRates(g3, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
plot(infoRates(g4, 2), infoRates(g4, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
plot(infoRates(g5, 2), infoRates(g5, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0.7 0]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.2;%0.5;
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-PLdvsPMvsPI-zoom-v%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by PLd vs PLv vs PM vs PI
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantEvokedSelectivity & isInPulvinar;
g2 = isSignificantEvokedSelectivity & strcmp(localization, 'PM');
g3 = isSignificantEvokedSelectivity & strcmp(localization, 'PLd');
g4 = isSignificantEvokedSelectivity & strcmp(localization, 'PLv');
g5 = isSignificantEvokedSelectivity & strcmp(localization, 'PI');

plot(meanNormSpdfInRFAllWindowsAll(g1, 4), meanNormSpdfInRFAllWindowsAll(g1, 6), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(meanNormSpdfInRFAllWindowsAll(g2, 4), meanNormSpdfInRFAllWindowsAll(g2, 6), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(meanNormSpdfInRFAllWindowsAll(g3, 4), meanNormSpdfInRFAllWindowsAll(g3, 6), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
plot(meanNormSpdfInRFAllWindowsAll(g4, 4), meanNormSpdfInRFAllWindowsAll(g4, 6), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
plot(meanNormSpdfInRFAllWindowsAll(g5, 4), meanNormSpdfInRFAllWindowsAll(g5, 6), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0.7 0]);
currXLim = xlim();
currYLim = ylim();
maxLim = max([currXLim currYLim]);
plot([-maxLim maxLim], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([-maxLim maxLim]);
ylim([-maxLim maxLim]);
axis square;
box off;
xlabel({'Normalized Firing Rate After Cue Onset'});
ylabel({'Normalized Firing Rate After Array Onset'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-cueResponseArrayResponse-PLdvsPMvsPI-zoom-v%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');


%% Cue-Target Delay Info vs Target-Dim Delay AI, Color by PLd vs PLv vs PM vs PI
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
% TODO for now use significance based on info rate BUT should use AI
% shuffle
g1 = isCell & ~isSignificantDelaySelectivity;% & isInPulvinar;
g2 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PI');

plot(attnIndices(g1, 1), attnIndices(g1, 2), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(attnIndices(g2, 1), attnIndices(g2, 2), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(attnIndices(g3, 1), attnIndices(g3, 2), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
plot(attnIndices(g4, 1), attnIndices(g4, 2), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.3; % set manually
assert(~any(any(abs(attnIndices(g1 | g2 | g3 | g4,:) > maxLim))));
plot([0 0], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([-maxLim maxLim]);
ylim([-maxLim maxLim]);
axis square;
box off;
xlabel({'Attention Index during', 'Cue-Target Delay'});
ylabel({'Attention Index during', 'Target-Dim Delay'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-AI1vsAI2-PLdvsPMvsPI-v%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by BS vs NS
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & isBroadSpiking & isInPulvinar;
g3 = isSignificantDelaySelectivity & isNarrowSpiking & isInPulvinar;

plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(infoRates(g3, 2), infoRates(g3, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = max([currXLim(2) currYLim(2)]);
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-BSvsNS-%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by BS vs NS - ZOOM
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & isBroadSpiking & isInPulvinar;
g3 = isSignificantDelaySelectivity & isNarrowSpiking & isInPulvinar;

plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(infoRates(g3, 2), infoRates(g3, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.2;
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-BSvsNS-zoom-%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by BS vs NS
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & isBroadSpiking & isInPulvinar;
g3 = isSignificantDelaySelectivity & isNarrowSpiking & isInPulvinar;

plot(attnIndices(g1, 1), attnIndices(g1, 2), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(attnIndices(g2, 1), attnIndices(g2, 2), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(attnIndices(g3, 1), attnIndices(g3, 2), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.7; % set manually
assert(~any(any(abs(attnIndices(g1 | g2 | g3 | g4,:) > maxLim))));
plot([0 0], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([-maxLim maxLim]);
ylim([-maxLim maxLim]);
axis square;
box off;
xlabel({'Attention Index during', 'Cue-Target Delay'});
ylabel({'Attention Index during', 'Target-Dim Delay'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-AI1vsAI2-BSvsNS-%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');






%% PCA on activity space by cell
superdivision = isCell & isInPulvinar;
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca([meanNormSpdfInRFAllWindowsAll(superdivision,[1:2 4:end])]);
fprintf('\n');
fprintf('PCA: %d variables, %d observations\n', size(pcaScore, 2), size(pcaScore, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', pcaPctExplained(1) + pcaPctExplained(2));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;

cols = lines(5);
cols = [cols(1,:); cols(3,:); cols(5,:); cols(2,:)];

subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
localizationSuper = localization(superdivision);
for j = 1:numel(subdivisions)
    subdivision = subdivisions{j};
    isInSubdivision = strcmp(localizationSuper, subdivision);
    fprintf('\t%s: %d cells\n', subdivision, sum(isInSubdivision));

    sh = scatter(pcaScore(isInSubdivision,1), pcaScore(isInSubdivision,2), 100, cols(j,:), 'MarkerFaceColor', cols(j,:));
    sh.MarkerFaceAlpha = 0.9;
end

% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% xlim([-25 120]);
% ylim([-30 60]);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');
% set(gca, 'FontWeight', 'bold');

plotFileName = sprintf('%s/allSessions-pcaByActivity-splitBySubdivision-v%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% PCA on activity space by cell
superdivision = isCell & isInPulvinar;
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca([meanNormSpdfInRFAllWindowsAll(superdivision,[1:2 4:end])]);
fprintf('\n');
fprintf('PCA: %d variables, %d observations\n', size(pcaScore, 2), size(pcaScore, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', pcaPctExplained(1) + pcaPctExplained(2));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;

cols = lines(5);
cols = [cols(1,:); cols(3,:); cols(5,:); cols(2,:)];

subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
localizationSuper = localization(superdivision);
isSignificantDelaySelectivitySuper = isSignificantDelaySelectivity(superdivision);
for j = 1:numel(subdivisions)
    subdivision = subdivisions{j};
    isInSubdivision = strcmp(localizationSuper, subdivision);
    fprintf('\t%s: %d cells\n', subdivision, sum(isInSubdivision));

    sh = scatter(pcaScore(isInSubdivision,1), pcaScore(isInSubdivision,2), 100, cols(j,:), 'MarkerFaceColor', cols(j,:));
    sh.MarkerFaceAlpha = 0.9;
end

% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% xlim([-25 120]);
% ylim([-30 60]);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');
% set(gca, 'FontWeight', 'bold');

plotFileName = sprintf('%s/allSessions-pcaByActivity-splitBySubdivision-delaySigHighlight-v%d.png', summaryDataDir, v);
export_fig(plotFileName, '-nocrop');

%% use t-sne (random position each run)

tsneVals = tsne([meanNormSpdfInRFAllWindowsAll(superdivision,[1:2 4:end])]);

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
for j = 1:numel(subdivisions)
    subdivision = subdivisions{j};
    isInSubdivision = strcmp(localizationSuper, subdivision);
    fprintf('\t%s: %d cells\n', subdivision, sum(isInSubdivision));

    sh = scatter(tsneVals(isInSubdivision,1), tsneVals(isInSubdivision,2), 100, cols(j,:), 'MarkerFaceColor', cols(j,:));
    sh.MarkerFaceAlpha = 0.9;
end


% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% xlim([-25 120]);
% ylim([-30 60]);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');
% set(gca, 'FontWeight', 'bold');
