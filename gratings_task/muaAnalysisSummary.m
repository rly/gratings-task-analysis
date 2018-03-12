function muaAnalysisSummary(processedDataRootDir, recordingInfoFileName)

v = 11;

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);

summaryDataDir = sprintf('%s/%s/', processedDataRootDir, 'MUA_GRATINGS_SUMMARY');
if ~exist(summaryDataDir, 'dir')
    error('No directory %s\n', summaryDataDir);
end

nUnitsApprox = 16;%numel(recordingInfo) * 16; % should be equal or an underestimate

unitNames = cell(nUnitsApprox, 1);
isSignificantResponseVsBaseline = false(nUnitsApprox, 6); % 6 periods > baseline, 5 periods info rate
isSignificantSelectivity = false(nUnitsApprox, 5); % 6 periods > baseline, 5 periods info rate
infoRates = nan(nUnitsApprox, 5); % 5 periods
attnIndices = nan(nUnitsApprox, 2); % 2 delay periods
localization = cell(nUnitsApprox, 1);
isInVPulvinar = false(nUnitsApprox, 1);
isInDPulvinar = false(nUnitsApprox, 1);
spdfInfo = struct();

meanRTHoldInRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldInRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelInRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelInRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldInRFTopThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
meanRTHoldInRFBottomThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelExRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelExRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFTopThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFBottomThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
corrCoefHoldInRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefRelInRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefHoldInRFTDDelayRTAll = nan(nUnitsApprox, 3);
corrCoefHoldExRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefRelExRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefHoldExRFTDDelayRTAll = nan(nUnitsApprox, 3);

unitCount = 0;
% should also be running a lot of shuffle tests given the number of trials

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');

%% session loop
for sessionInd = 7%1:numel(recordingInfo)
    sessionName = recordingInfo(sessionInd).sessionName;
    saveFileName = sprintf('%s/%s-sessionInd%d-muaAnalysisSummaryData-v%d.mat', summaryDataDir, sessionName, sessionInd, v);
    fprintf('Loading file %s ...\n', saveFileName);
    S = load(saveFileName);
    
    fprintf('Found %d units...\n', numel(S.unitNames));
    currentUnitInds = (unitCount + 1):(unitCount + numel(S.unitNames));
    unitCount = unitCount + numel(S.unitNames);
    unitNames(currentUnitInds) = S.unitNames;
    isSignificantResponseVsBaseline(currentUnitInds,:) = S.isSignificantResponseVsBaseline;
    isSignificantSelectivity(currentUnitInds,:) = S.isSignificantSelectivity;
    infoRates(currentUnitInds,:) = S.infoRates;
    attnIndices(currentUnitInds,:) = S.attnIndices;
    localization(currentUnitInds) = S.localization;
    isInVPulvinar(currentUnitInds) = S.isInVPulvinar;
    isInDPulvinar(currentUnitInds) = S.isInDPulvinar;
    
    fn = fieldnames(S.spdfInfo);
    for j = 1:numel(fn)
        if isfield(spdfInfo, fn{j})
            spdfInfo.(fn{j})(currentUnitInds,:) = S.spdfInfo.(fn{j});
        else
            spdfInfo.(fn{j}) = S.spdfInfo.(fn{j}); % no pre-allocation
        end
    end
    
    % overwrite each session but that's ok. they should all be the same
    enterFixationT = S.enterFixationT;
    cueOnsetT = S.cueOnsetT;
    arrayOnsetT = S.arrayOnsetT;
    targetDimT = S.targetDimT;
    exitFixationT = S.exitFixationT;
end
clear S;

%% summarize
nUnitsAll = unitCount;
isCell = true(unitCount, 1); % for MUA, cannot distinguish between cell and not cell

% TODO compute these per unit above so that the right alpha is used
isSignificantAnyTaskMod = isCell & any(isSignificantResponseVsBaseline, 2);
isSignificantCueResponse = isCell & isSignificantResponseVsBaseline(:,1);
isSignificantCueTargetDelay = isCell & isSignificantSelectivity(:,2);
isSignificantTargetDimDelay = isCell & isSignificantSelectivity(:,4);
isSignificantPreExitFixation = isCell & isSignificantResponseVsBaseline(:,6);

isSignificantAnySpatialSelectivity = isCell & any(isSignificantSelectivity, 2);
isSignificantEvokedSelectivity = isCell & any(isSignificantSelectivity(:,[1 3 5]), 2);
isSignificantDelaySelectivity = isCell & any(isSignificantSelectivity(:,[2 4]), 2);
% isSigSelectEvokedNotDelay = isCell & any(isSignificantSelectivity(:,[1 3 5]), 2) & any(isSignificantSelectivity(:,[2 4]), 2);
% isSigSelectOnlyEvoked = isCell & any(isSignificantSelectivity(:,[1 3 5]), 2) & ~any(isSignificantSelectivity(:,[2 4]), 2);
% isSigSelectOnlyDelay = isCell & ~any(isSignificantSelectivity(:,[1 3 5]), 2) & any(isSignificantSelectivity(:,[2 4]), 2);
isInPulvinar = strcmp(localization, 'vPul') | strcmp(localization, 'dPul');
% strcmp(localization, 'PLd') | strcmp(localization, 'PLv') | strcmp(localization, 'PM') | strcmp(localization, 'PI');

fprintf('-------------------------------\n');
fprintf('%d/%d = %d%% units were localized to the pulvinar.\n', ...
        sum(isCell & isInPulvinar), nUnitsAll, ...
        round(sum(isCell & isInPulvinar)/nUnitsAll * 100));
fprintf('\n');

fprintf('%d/%d = %d%% units show significant task modulation compared to baseline (not corrected for multiple comparisons).\n', ...
        sum(isSignificantAnyTaskMod), nUnitsAll, ...
        round(sum(isSignificantAnyTaskMod)/nUnitsAll * 100));
fprintf('%d/%d = %d%% units show significant cue response compared to baseline.\n', ...
        sum(isSignificantCueResponse), nUnitsAll, ...
        round(sum(isSignificantCueResponse)/nUnitsAll * 100));
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

fprintf('Of the %d units that show significant spatial selectivity during some task period,\n', ...
        sum(isSignificantAnySpatialSelectivity));
% fprintf('\t%d (%d%%) show significant spatial selectivity after a visual change AND during attention,\n', ...
%         sum(isSigSelectEvokedNotDelay), ...
%         round(sum(isSigSelectEvokedNotDelay)/sum(isSignificantAnySpatialSelectivity) * 100));
% fprintf('\t%d (%d%%) show significant spatial selectivity ONLY after a visual change, and\n', ...
%         sum(isSigSelectOnlyEvoked), ...
%         round(sum(isSigSelectOnlyEvoked)/sum(isSignificantAnySpatialSelectivity) * 100));
% fprintf('\t%d (%d%%) show significant spatial selectivity ONLY during attention.\n\n', ...
%         sum(isSigSelectOnlyDelay), ...
%         round(sum(isSigSelectOnlyDelay)/sum(isSignificantAnySpatialSelectivity) * 100));
% fprintf('Of the %d units in the pulvinar that show spatial selectivity during attention, \n\t%d (%d%%) are in PLd, %d (%d%%) are in PLv, %d (%d%%) are in PM, and %d (%d%%) are in PI\n', ...
%         sum(isSignificantDelaySelectivity & isInPulvinar), ...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PM')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PM'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PI')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PI'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100));
fprintf('Of the %d units in the pulvinar that show significant cue response vs baseline, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantCueResponse & isInPulvinar), ...
        sum(isSignificantCueResponse & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantCueResponse & strcmp(localization, 'vPul'))/sum(isSignificantCueResponse & isInPulvinar) * 100), ...
        sum(isSignificantCueResponse & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantCueResponse & strcmp(localization, 'dPul'))/sum(isSignificantCueResponse & isInPulvinar) * 100));
fprintf('\n');
fprintf('Of the %d units in the pulvinar that show significant pre-saccade activity vs baseline, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100));
fprintf('\n');
fprintf('Of the %d units in the pulvinar that show spatial selectivity during attention, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100));
fprintf('\n');


fprintf('Of the %d units in the pulvinar that show significant cue response, \n', sum(isSignificantCueResponse & isInPulvinar));
fprintf('\t%d (%d%%) show significant selectivity during the cue-target delay\n', sum(isSignificantCueResponse & isInPulvinar & isSignificantCueTargetDelay), ...
        round(sum(isSignificantCueResponse & isInPulvinar & isSignificantCueTargetDelay)/sum(isSignificantCueResponse & isInPulvinar) * 100));
fprintf('\t%d (%d%%) show significant selectivity during the target-dim delay\n', sum(isSignificantCueResponse & isInPulvinar & isSignificantTargetDimDelay), ...
        round(sum(isSignificantCueResponse & isInPulvinar & isSignificantTargetDimDelay)/sum(isSignificantCueResponse & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d units in the pulvinar that show spatial selectivity during attention, \n\t%d (%d%%) are in vPul, %d (%d%%) are in dPul\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100));
fprintf('\n');


% fprintf('\t%d/%d = %d%% PLd units\n', ...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd')), ...
%         sum(isCell & strcmp(localization, 'PLd')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
% fprintf('\t%d/%d = %d%% PLv units\n', ...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv')), ...
%         sum(isCell & strcmp(localization, 'PLv')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
% fprintf('\t%d/%d = %d%% PM units\n',...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PM')), ...
%         sum(isCell & strcmp(localization, 'PM')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
% fprintf('\t%d/%d = %d%% PI units\n',...
%         sum(isSignificantDelaySelectivity & strcmp(localization, 'PI')), ...
%         sum(isCell & strcmp(localization, 'PI')), ...
%         round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in cue-target delay:\n');
% fprintf('\t%d/%d = %d%% PLd units\n', ...
%         sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PLd')), ...
%         sum(isCell & strcmp(localization, 'PLd')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
% fprintf('\t%d/%d = %d%% PLv units\n', ...
%         sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PLv')), ...
%         sum(isCell & strcmp(localization, 'PLv')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
% fprintf('\t%d/%d = %d%% PM units\n',...
%         sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PM')), ...
%         sum(isCell & strcmp(localization, 'PM')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
% fprintf('\t%d/%d = %d%% PI units\n',...
%         sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PI')), ...
%         sum(isCell & strcmp(localization, 'PI')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantCueTargetDelay & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in target-dim delay:\n');
% fprintf('\t%d/%d = %d%% PLd units\n', ...
%         sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PLd')), ...
%         sum(isCell & strcmp(localization, 'PLd')), ...
%         round(sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
% fprintf('\t%d/%d = %d%% PLv units\n', ...
%         sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PLv')), ...
%         sum(isCell & strcmp(localization, 'PLv')), ...
%         round(sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
% fprintf('\t%d/%d = %d%% PM units\n',...
%         sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PM')), ...
%         sum(isCell & strcmp(localization, 'PM')), ...
%         round(sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
% fprintf('\t%d/%d = %d%% PI units\n',...
%         sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PI')), ...
%         sum(isCell & strcmp(localization, 'PI')), ...
%         round(sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantTargetDimDelay & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

fprintf('Significant selectivity in both delay periods:\n');
% fprintf('\t%d/%d = %d%% PLd units\n', ...
%         sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PLd')), ...
%         sum(isCell & strcmp(localization, 'PLd')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
% fprintf('\t%d/%d = %d%% PLv units\n', ...
%         sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PLv')), ...
%         sum(isCell & strcmp(localization, 'PLv')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
% fprintf('\t%d/%d = %d%% PM units\n',...
%         sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PM')), ...
%         sum(isCell & strcmp(localization, 'PM')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
% fprintf('\t%d/%d = %d%% PI units\n',...
%         sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PI')), ...
%         sum(isCell & strcmp(localization, 'PI')), ...
%         round(sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\t%d/%d = %d%% vPul units\n',...
        sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'vPul')), ...
        sum(isCell & strcmp(localization, 'vPul')), ...
        round(sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'vPul'))/sum(isCell & strcmp(localization, 'vPul')) * 100));
fprintf('\t%d/%d = %d%% dPul units\n',...
        sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'dPul')), ...
        sum(isCell & strcmp(localization, 'dPul')), ...
        round(sum(isCell & isSignificantCueTargetDelay & isSignificantTargetDimDelay & strcmp(localization, 'dPul'))/sum(isCell & strcmp(localization, 'dPul')) * 100));
fprintf('\n');

% stop

%% per-condition baseline-corrected normalized mean
fprintf('\n');
% subdivisions = {'PM', 'PLd', 'PLv', 'PI', 'PUL', 'all'};
subdivisions = {'all', 'vPul', 'notVPul', 'dPul', 'notDPul'};
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

    plotFileName = sprintf('%s/allSessions-%s-meanSpdfs3-v%d.png', summaryDataDir, subdivision, v);
    fprintf('Saving to %s...\n', plotFileName);
    
    quickSpdfAllEvents3InARowPopMean(cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
            arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
            targetDimSpdfExRFNormSub, cueOnsetT, arrayOnsetT, targetDimT, plotFileName);

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
% subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
subdivisions = {'all', 'vPul', 'dPul'};%, 'notVPul', 'notDPul'};
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



%% test diff RT for top third vs bottom third firing rates in delay periods
meanRTRelInRFDiffThirdFiringRateCTDelayAll = meanRTRelInRFTopThirdFiringRateCTDelayAll - meanRTRelInRFBottomThirdFiringRateCTDelayAll;
meanRTRelExRFDiffThirdFiringRateCTDelayAll = meanRTRelExRFTopThirdFiringRateCTDelayAll - meanRTRelExRFBottomThirdFiringRateCTDelayAll;
meanRTHoldInRFDiffThirdFiringRateCTDelayAll = meanRTHoldInRFTopThirdFiringRateCTDelayAll - meanRTHoldInRFBottomThirdFiringRateCTDelayAll;
meanRTHoldExRFDiffThirdFiringRateCTDelayAll = meanRTHoldExRFTopThirdFiringRateCTDelayAll - meanRTHoldExRFBottomThirdFiringRateCTDelayAll;
meanRTHoldInRFDiffThirdFiringRateTDDelayAll = meanRTHoldInRFTopThirdFiringRateTDDelayAll - meanRTHoldInRFBottomThirdFiringRateTDDelayAll;
meanRTHoldExRFDiffThirdFiringRateTDDelayAll = meanRTHoldExRFTopThirdFiringRateTDDelayAll - meanRTHoldExRFBottomThirdFiringRateTDDelayAll;

meanRTRelInRFDiffThirdFiringRateCTDelayPul = meanRTRelInRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTRelExRFDiffThirdFiringRateCTDelayPul = meanRTRelExRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTHoldInRFDiffThirdFiringRateCTDelayPul = meanRTHoldInRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTHoldExRFDiffThirdFiringRateCTDelayPul = meanRTHoldExRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTHoldInRFDiffThirdFiringRateTDDelayPul = meanRTHoldInRFDiffThirdFiringRateTDDelayAll(isCell & isInPulvinar);
meanRTHoldExRFDiffThirdFiringRateTDDelayPul = meanRTHoldExRFDiffThirdFiringRateTDDelayAll(isCell & isInPulvinar);

[~,p] = ttest(meanRTRelInRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTRelExRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTHoldInRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTHoldExRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTHoldInRFDiffThirdFiringRateTDDelayPul)
[~,p] = ttest(meanRTHoldExRFDiffThirdFiringRateTDDelayPul)


cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);
binEdges = -0.1:0.01:0.1;

figure_tr_inch(9, 6);
set(gcf, 'Color', 'w');
plotHs = nan(6, 1);
plotHs(1) = subaxis(2, 3, 1);
hold on;
hist1 = histogram(meanRTRelInRFDiffThirdFiringRateCTDelayPul, binEdges);
hist1.FaceColor = inRFCol;
plotHs(2) = subaxis(2, 3, 2);
hold on;
hist2 = histogram(meanRTHoldInRFDiffThirdFiringRateCTDelayPul, binEdges);
hist2.FaceColor = inRFCol;
plotHs(3) = subaxis(2, 3, 3);
hold on;
hist3 = histogram(meanRTHoldInRFDiffThirdFiringRateTDDelayPul, binEdges);
hist3.FaceColor = inRFCol;
plotHs(4) = subaxis(2, 3, 4);
hold on;
hist4 = histogram(meanRTRelExRFDiffThirdFiringRateCTDelayPul, binEdges);
hist4.FaceColor = exRFCol;
plotHs(5) = subaxis(2, 3, 5);
hold on;
hist5 = histogram(meanRTHoldExRFDiffThirdFiringRateCTDelayPul, binEdges);
hist5.FaceColor = exRFCol;
plotHs(6) = subaxis(2, 3, 6);
hold on;
hist6 = histogram(meanRTHoldExRFDiffThirdFiringRateTDDelayPul, binEdges);
hist6.FaceColor = exRFCol;

% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);

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

corrCoefHoldInRFCTDelayRTPul = corrCoefHoldInRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefRelInRFCTDelayRTPul = corrCoefRelInRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefHoldInRFTDDelayRTPul = corrCoefHoldInRFTDDelayRTAll(isCell & isInPulvinar,3);
corrCoefHoldExRFCTDelayRTPul = corrCoefHoldExRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefRelExRFCTDelayRTPul = corrCoefRelExRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefHoldExRFTDDelayRTPul = corrCoefHoldExRFTDDelayRTAll(isCell & isInPulvinar,3);

[~,p] = ttest(corrCoefRelInRFCTDelayRTPul)
[~,p] = ttest(corrCoefRelExRFCTDelayRTPul)
[~,p] = ttest(corrCoefHoldInRFCTDelayRTPul)
[~,p] = ttest(corrCoefHoldExRFCTDelayRTPul)
[~,p] = ttest(corrCoefHoldInRFTDDelayRTPul)
[~,p] = ttest(corrCoefHoldExRFTDDelayRTPul)


cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);
binEdges = -0.3:0.05:0.3;

figure_tr_inch(9, 6);
set(gcf, 'Color', 'w');
plotHs = nan(6, 1);
plotHs(1) = subaxis(2, 3, 1);
hold on;
hist1 = histogram(corrCoefRelInRFCTDelayRTPul, binEdges);
hist1.FaceColor = inRFCol;
plotHs(2) = subaxis(2, 3, 2);
hold on;
hist2 = histogram(corrCoefHoldInRFCTDelayRTPul, binEdges);
hist2.FaceColor = inRFCol;
plotHs(3) = subaxis(2, 3, 3);
hold on;
hist3 = histogram(corrCoefHoldInRFTDDelayRTPul, binEdges);
hist3.FaceColor = inRFCol;
plotHs(4) = subaxis(2, 3, 4);
hold on;
hist4 = histogram(corrCoefRelExRFCTDelayRTPul, binEdges);
hist4.FaceColor = exRFCol;
plotHs(5) = subaxis(2, 3, 5);
hold on;
hist5 = histogram(corrCoefHoldExRFCTDelayRTPul, binEdges);
hist5.FaceColor = exRFCol;
plotHs(6) = subaxis(2, 3, 6);
hold on;
hist6 = histogram(corrCoefHoldExRFTDDelayRTPul, binEdges);
hist6.FaceColor = exRFCol;

% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);


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
