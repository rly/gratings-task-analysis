
clear;
readDataLocally;
sessionInds = 1:23;

v = 12;
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\GRC\';

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);

summaryDataDir = sprintf('%s/MUA_GRATINGS_SUMMARY/', processedDataRootDir);
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
% isSignificantResponseVsBootstrapBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
% isSignificantResponseVsBootstrapPreviousPeriod = false(nUnitsApprox, 4);
isSignificantSelectivity = false(nUnitsApprox, 5); % 5 periods info rate
cueResponseVsBaselineDirection = zeros(nUnitsApprox, 1);
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
    
    unitNames(currentUnitInds) = S.unitNames; % TODO include sessionInd
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
isSignificantSelectivityCueResponse = isCell & isSignificantSelectivity(:,1);
isSignificantSelectivityCueTargetDelay = isCell & isSignificantSelectivity(:,2);
isSignificantSelectivityArrayHoldResponse = isCell & isSignificantSelectivity(:,3);
isSignificantSelectivityTargetDimDelay = isCell & isSignificantSelectivity(:,4);
isSignificantSelectivityTargetDimResponse = isCell & isSignificantSelectivity(:,5);

isSignificantSelectivityCueTargetDelayInc = isCell & isSignificantSelectivity(:,2) & diffRates(:,1) > 0;
isSignificantSelectivityCueTargetDelayDec = isCell & isSignificantSelectivity(:,2) & diffRates(:,1) < 0;
isSignificantSelectivityArrayHoldResponseInc = isCell & isSignificantSelectivity(:,3) & diffRates(:,2) > 0;
isSignificantSelectivityArrayHoldResponseDec = isCell & isSignificantSelectivity(:,3) & diffRates(:,2) < 0;
isSignificantSelectivityTargetDimDelayInc = isCell & isSignificantSelectivity(:,4) & diffRates(:,3) > 0;
isSignificantSelectivityTargetDimDelayDec = isCell & isSignificantSelectivity(:,4) & diffRates(:,3) < 0;

isInRFP1 = isCell & inRFLocs == 1;
isInRFP3 = isCell & inRFLocs == 3;

isInPulvinar = strcmp(localization, 'vPul') | strcmp(localization, 'dPul');

cols = lines(6);

%% firing rate & fano factor loop across subdivisions
subdivisions = {'dPul2', 'vPul2'};
for i = 1:numel(subdivisions)
    if strcmp(subdivisions{i}, 'dPul')
        goodUnits = isInDPulvinar & isSignificantCueResponseInc;
    elseif strcmp(subdivisions{i}, 'vPul')
        goodUnits = isInVPulvinar & isSignificantCueResponseInc;
    elseif strcmp(subdivisions{i}, 'dPul2')
        goodUnits = isInDPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse;
    elseif strcmp(subdivisions{i}, 'vPul2')
        goodUnits = isInVPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse;
    elseif strcmp(subdivisions{i}, 'dPul3')
        goodUnits = isInDPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse & isInRFP3;
    elseif strcmp(subdivisions{i}, 'vPul3')
        goodUnits = isInVPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse & isInRFP3;
    end

    %% cue target delay mean firing rate InRF vs ExRF
    [~,ax1] = grcPlotRateDiff(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
            inRFLocs(goodUnits), exRFLocs(goodUnits), ...
            isInDPulvinar(goodUnits), isInVPulvinar(goodUnits), ...
            isSignificantSelectivityCueTargetDelay(goodUnits));
    xlabel(ax1, 'Attend-RF Firing Rate (Hz)');
    ylabel(ax1, 'Attend-Away Firing Rate (Hz)');

    pctSig = round(sum(isSignificantSelectivityCueTargetDelay(goodUnits)) / sum(goodUnits) * 100);

    legend(ax1.Children(1:2), {sprintf('Sig. Modulation (%d%%)', pctSig), 'N.S. Modulation'}, ...
            'Location', 'NorthWest', 'FontSize', 14, 'box', 'off');

    plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayMeanFRDiff-v%d.png', outputDir, subdivisions{i}, v);
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');

    %% cue target delay fano factor InRF vs ExRF
    [~,ax1] = grcPlotFanoFactorDiff(averageFiringRatesByCount.cueTargetDelay(goodUnits), ...
            inRFLocs(goodUnits), exRFLocs(goodUnits), ...
            isInDPulvinar(goodUnits), isInVPulvinar(goodUnits), ...
            false(sum(goodUnits), 1));
    xlabel(ax1, 'Attend-RF Fano Factor');
    ylabel(ax1, 'Attend-Away Fano Factor');

    plotFileName = sprintf('%s/allSessions-%s-cueTargetDelayFanoFactorDiff-v%d.png', outputDir, subdivisions{i}, v);
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');

    %% array hold response mean firing rate InRF vs ExRF
    [~,ax1] = grcPlotRateDiff(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
            inRFLocs(goodUnits), exRFLocs(goodUnits), ...
            isInDPulvinar(goodUnits), isInVPulvinar(goodUnits), ...
            isSignificantSelectivityArrayHoldResponse(goodUnits));
    xlabel(ax1, 'Attend-RF Firing Rate (Hz)');
    ylabel(ax1, 'Attend-Away Firing Rate (Hz)');
    
    pctSig = round(sum(isSignificantSelectivityArrayHoldResponse(goodUnits)) / sum(goodUnits) * 100);
    
    legend(ax1.Children(1:2), {sprintf('Sig. Modulation (%d%%)', pctSig), 'N.S. Modulation'}, ...
            'Location', 'NorthWest', 'FontSize', 14, 'box', 'off');

    plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldMeanFRDiff-v%d.png', outputDir, subdivisions{i}, v);
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
    
    %% array hold response fano factor InRF vs ExRF
    [~,ax1] = grcPlotFanoFactorDiff(averageFiringRatesByCount.arrayResponseHoldBal(goodUnits), ...
            inRFLocs(goodUnits), exRFLocs(goodUnits), ...
            isInDPulvinar(goodUnits), isInVPulvinar(goodUnits), ...
            false(sum(goodUnits), 1));
    xlabel(ax1, 'Attend-RF Fano Factor');
    ylabel(ax1, 'Attend-Away Fano Factor');

    plotFileName = sprintf('%s/allSessions-%s-arrayResponseHoldFanoFactorDiff-v%d.png', outputDir, subdivisions{i}, v);
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');

end


%% cue target delay noise correlations P3 vs P1
goodUnitsDPul = isInDPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse & isInRFP3 & [averageFiringRatesByCount.cueTargetDelay.all]' >= 0;
goodUnitsVPul = isInVPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse & isInRFP3 & [averageFiringRatesByCount.cueTargetDelay.all]' >= 0;

[~,ax1] = grcPlotNoiseCorrDiff(cueTargetDelayNoiseCorr, goodUnitsDPul, false(size(goodUnitsVPul)));
xlabel(ax1, 'Noise Correlation Attend-RF');
ylabel(ax1, 'Noise Correlation Attend-Away');

plotFileName = sprintf('%s/allSessions-dPul2-cueTargetDelayNoiseCorrDiff-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

[~,ax1] = grcPlotNoiseCorrDiff(cueTargetDelayNoiseCorr, false(size(goodUnitsDPul)), goodUnitsVPul);
xlabel(ax1, 'Noise Correlation Attend-RF');
ylabel(ax1, 'Noise Correlation Attend-Away');

plotFileName = sprintf('%s/allSessions-vPul2-cueTargetDelayNoiseCorrDiff-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% array hold response mid noise correlations P3 vs P1
goodUnitsDPul = isInDPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse & isInRFP3 & [averageFiringRatesByCount.arrayResponseHoldBal.all]' >= 0;
goodUnitsVPul = isInVPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse & isInRFP3 & [averageFiringRatesByCount.arrayResponseHoldBal.all]' >= 0;

[~,ax1] = grcPlotNoiseCorrDiff(arrayResponseHoldMidNoiseCorr, goodUnitsDPul, false(size(goodUnitsVPul)));
xlabel(ax1, 'Noise Correlation Attend-RF');
ylabel(ax1, 'Noise Correlation Attend-Away');

plotFileName = sprintf('%s/allSessions-dPul2-arrayResponseHoldMidNoiseCorrDiff-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

[~,ax1] = grcPlotNoiseCorrDiff(arrayResponseHoldMidNoiseCorr, false(size(goodUnitsDPul)), goodUnitsVPul);
xlabel(ax1, 'Noise Correlation Attend-RF');
ylabel(ax1, 'Noise Correlation Attend-Away');

plotFileName = sprintf('%s/allSessions-vPul2-arrayResponseHoldMidNoiseCorrDiff-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% population mean plots
grcPopMeanMuaSpdfCueOnsetArrayOnsetHoldExitFixation;

%% where are the example cells in the FR scatterplots
match1 = strcmp(unitNames, 'M20170311_PUL_48M');
match2 = strcmp(unitNames, 'M20170331_PUL_46M');

matches = [match1 match2];
for i = 1:size(matches, 2)
    m = matches(:,i);
    if sum(m) > 1 % two different blocks
        m = false(size(m));
        m(find(matches(:,i), 1)) = true;
    end
    inRFLoc = inRFLocs(m);
    exRFLoc = exRFLocs(m);
    fprintf('Unit: %s\n', unitNames{m});
    fr = averageFiringRatesByCount.cueTargetDelay(m);
    fprintf('CT Delay Firing InRF: %0.1f, ExRF: %0.1f\n', fr.byLoc(inRFLoc), fr.byLoc(exRFLoc));
    fr = averageFiringRatesByCount.arrayResponseHoldBal(m);
    fprintf('Array Response Hold Firing InRF: %0.1f, ExRF: %0.1f\n\n', fr.byLoc(inRFLoc), fr.byLoc(exRFLoc));
end