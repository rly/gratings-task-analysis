function usefulEvents = getUsefulEvents2(logDir, logIndices, nLoc, D, blockName)

holdDurMid = 850; % ms; for splitting short vs long hold
logRoot = 'gratings_task_eye_tracking_';

% cueP3Event = D.events{2}(:,1);
cueOnsetEvent = D.events{4}(:,1);
arrayOnsetHoldEvent = D.events{6}(:,1); % TODO
arrayOnsetReleaseEvent = D.events{7}(:,1);
% arrayOnsetHoldP1Event = D.events{14}(:,1);
% arrayOnsetReleaseP1Event = D.events{15}(:,1);
targetDimEvent = D.events{16}(:,1);
juiceEvent = D.events{8}(:,1);

%%
trialParamsAllCorrect = readPresentationLogCorrectTrial(logDir, logRoot, logIndices);

%% load more detailed trial info from json -> mat
% use saveGratingsTaskResultsJsonToMatRunner.m to generate mat files
load(sprintf('%s/gratingsTaskResultsJson-%s.mat', logDir, blockName), 'trialStructs', 'arrayShapes');
trialResults = cellfun(@(x) x.response.type, trialStructs, 'UniformOutput', false);
isCorrect = strcmp(trialResults, 'correct-response');
trialStructsCorrect = trialStructs(isCorrect);
verifyGratingsTaskLogsAndJsonCorrect(trialParamsAllCorrect, trialStructsCorrect);
arrayShapesCorrect = arrayShapes(isCorrect);
isArrayRelBalanced = cellfun(@(x) x(1) == 'R' && x(3) == 'R', arrayShapesCorrect)';
isArrayHoldBalanced = cellfun(@(x) x(1) == 'H' && x(3) == 'H', arrayShapesCorrect)';

isCorrectRelTrial = cellfun(@(x,y) y(x.cueLoc) == 'R', trialStructsCorrect, arrayShapesCorrect)';
isCorrectHoldTrial = cellfun(@(x,y) y(x.cueLoc) == 'H', trialStructsCorrect, arrayShapesCorrect)';

%% clean up juice events
% store only the first juice event
firstJuiceEvent = nan(1, 1);
count = 0;
diffJuiceEvent = diff(juiceEvent);
maxDiffJuiceEventSameTrial = 1; % seconds
for i = 1:numel(juiceEvent)
    if i >= 2 && diffJuiceEvent(i-1) < maxDiffJuiceEventSameTrial
        continue;
    end
    count = count + 1;
    firstJuiceEvent(count) = juiceEvent(i);
end

firstJuiceEvent = firstJuiceEvent';
fprintf('There are %d correct trials in the log and %d correct trials based on juice events.\n', ...
        size(trialParamsAllCorrect, 1), numel(firstJuiceEvent));
assert(size(trialParamsAllCorrect, 1) == numel(firstJuiceEvent));

%% get cue onsets corresponding to correct trials only
maxCueOnsetToJuiceTime = 5; % seconds
cueOnset = nan(numel(firstJuiceEvent), 1);
for i = 1:numel(firstJuiceEvent)  
    prevCueOnsetInd = find((firstJuiceEvent(i) - cueOnsetEvent < maxCueOnsetToJuiceTime) & ...
            (firstJuiceEvent(i) - cueOnsetEvent > 0), 1, 'last');
    if ~isempty(prevCueOnsetInd)
        cueOnset(i) = cueOnsetEvent(prevCueOnsetInd);
    else
        warning('Juice event without a corresponding cue onset event: %d, %f', i, firstJuiceEvent(i));
    end
end

%% remove trials that do not have an earlier lever press
[cueOnset,trialsToKeep] = removeHandStartedTrials(D, cueOnset, firstJuiceEvent);
firstJuiceEvent = firstJuiceEvent(trialsToKeep);
trialParamsAllCorrect = trialParamsAllCorrect(trialsToKeep,:);
fprintf('Removed %d/%d putative hand-started trials.\n', sum(~trialsToKeep), numel(trialsToKeep));

%%
% TODO get only trials that are not repeats
% notRepeatLogical = trialParamsAllCorrect(:,3) == 1;
% juiceEvent = juiceEvent(notRepeatLogical,:);
trialParams = trialParamsAllCorrect;%(notRepeatLogical,:);
% numNonRptCorrectTrials = size(trialParams, 1);

cueLoc = trialParams(:,4);%(trialParams(:,4) + 3) / 2;
% here, P1 is btm right, P2 is top right, P3 is top left, P4 is btm left
cueTargetDelayDur = trialParams(:,6);
targetDimDelayDur = trialParams(:,7);

%% get cue onset events
cueOnsetByLoc = cell(nLoc, 1);
for i = 1:nLoc
    cueOnsetByLoc{i} = cueOnset(cueLoc == i);
end    

%% get array onsets corresponding to correct trials only
maxArrayOnsetToJuiceTime = 2.5; % seconds
arrayOnset = nan(numel(firstJuiceEvent), 1);
allArrayOnsetEvents = [arrayOnsetHoldEvent; arrayOnsetReleaseEvent];
for i = 1:numel(firstJuiceEvent) 
    prevArrayOnsetInd = find((firstJuiceEvent(i) - allArrayOnsetEvents < maxArrayOnsetToJuiceTime) & ...
            (firstJuiceEvent(i) - allArrayOnsetEvents > 0), 1, 'last');
    if ~isempty(prevArrayOnsetInd)
        arrayOnset(i) = allArrayOnsetEvents(prevArrayOnsetInd);
    else
        warning('Juice event without a corresponding array onset event: %d, %f', i, firstJuiceEvent(i));
    end
end

% split release and hold shapes - not the best way, but works for now
maxArrayOnsetToJuiceTimeReleaseShape = 0.925;
isHoldTrial = firstJuiceEvent - arrayOnset >= maxArrayOnsetToJuiceTimeReleaseShape;
assert(all(isHoldTrial == (trialParamsAllCorrect(:,7) ~= -1)));
% assert(all(isHoldTrial == isCorrectHoldTrial)); % FIXME

arrayOnsetRel = arrayOnset(~isHoldTrial);
arrayOnsetHold = arrayOnset(isHoldTrial);

arrayOnsetRelBal = arrayOnset(~isHoldTrial & isArrayRelBalanced);
arrayOnsetHoldBal = arrayOnset(isHoldTrial & isArrayHoldBalanced);

arrayOnsetByLoc = cell(nLoc, 1);
arrayOnsetRelByLoc = cell(nLoc, 1);
arrayOnsetHoldByLoc = cell(nLoc, 1);
arrayOnsetBalByLoc = cell(nLoc, 1);
arrayOnsetRelBalByLoc = cell(nLoc, 1);
arrayOnsetHoldBalByLoc = cell(nLoc, 1);
arrayOnsetHoldBalShortHoldByLoc = cell(nLoc, 1);
arrayOnsetHoldBalLongHoldByLoc = cell(nLoc, 1);

for i = 1:nLoc
    arrayOnsetByLoc{i} = arrayOnset(cueLoc == i);
    arrayOnsetRelByLoc{i} = arrayOnset(~isHoldTrial & cueLoc == i);
    arrayOnsetHoldByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i);
    arrayOnsetBalByLoc{i} = arrayOnset(cueLoc == i & (isArrayRelBalanced | isArrayHoldBalanced));
    arrayOnsetRelBalByLoc{i} = arrayOnset(~isHoldTrial & cueLoc == i & isArrayRelBalanced);
    arrayOnsetHoldBalByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & isArrayHoldBalanced);
    arrayOnsetHoldBalShortHoldByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & targetDimDelayDur < holdDurMid & isArrayHoldBalanced);
    arrayOnsetHoldBalLongHoldByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & targetDimDelayDur >= holdDurMid & isArrayHoldBalanced);
end

cueLocHold = cueLoc(isHoldTrial);
cueLocRel = cueLoc(~isHoldTrial);

cueOnsetRel = cueOnset(~isHoldTrial);
cueOnsetHold = cueOnset(isHoldTrial);

cueOnsetRelByLoc = cell(nLoc, 1);
cueOnsetHoldByLoc = cell(nLoc, 1);
for i = 1:nLoc
    cueOnsetRelByLoc{i} = cueOnset(~isHoldTrial & cueLoc == i);
    cueOnsetHoldByLoc{i} = cueOnset(~isHoldTrial & cueLoc == i);
end 

%% get target dim events corresponding to correct trials only
% TODO use time of lever response instead of juice onset
maxTargetDimToJuiceTime = 1; % seconds
targetDimMatch = nan(numel(firstJuiceEvent), 1);
rt = nan(numel(firstJuiceEvent), 1);
for i = 1:numel(firstJuiceEvent)  
    prevEvent7Ind = find((firstJuiceEvent(i) - targetDimEvent < maxTargetDimToJuiceTime) & ...
            (firstJuiceEvent(i) - targetDimEvent > 0), 1, 'last');
    if ~isempty(prevEvent7Ind)
        targetDimMatch(i) = targetDimEvent(prevEvent7Ind);
        rt(i) = firstJuiceEvent(i) - targetDimMatch(i);
    else
        rt(i) = firstJuiceEvent(i) - arrayOnset(i);
    end
end

targetDim = targetDimMatch(~isnan(targetDimMatch));
targetDimBal = targetDimMatch(~isnan(targetDimMatch) & isArrayHoldBalanced);
targetDimBalShortHold = targetDimMatch(~isnan(targetDimMatch) & targetDimDelayDur < holdDurMid & isArrayHoldBalanced);
targetDimBalLongHold = targetDimMatch(~isnan(targetDimMatch) & targetDimDelayDur >= holdDurMid & isArrayHoldBalanced);
assert(numel(unique(targetDim)) == numel(targetDim));

targetDimByLoc = cell(nLoc, 1);
targetDimBalByLoc = cell(nLoc, 1);
targetDimBalShortHoldByLoc = cell(nLoc, 1);
targetDimBalLongHoldByLoc = cell(nLoc, 1);
nTrialShortHold = sum(~isnan(targetDimMatch) & targetDimDelayDur < holdDurMid & isArrayHoldBalanced);
nTrialLongHold = sum(~isnan(targetDimMatch) & targetDimDelayDur >= holdDurMid & isArrayHoldBalanced);
for i = 1:nLoc
    targetDimByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & cueLoc == i);
    targetDimBalByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & cueLoc == i & isArrayHoldBalanced);
    targetDimBalShortHoldByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & ...
            cueLoc == i & targetDimDelayDur < holdDurMid & isArrayHoldBalanced);
    targetDimBalLongHoldByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & ...
            cueLoc == i & targetDimDelayDur >= holdDurMid & isArrayHoldBalanced);
end

%%
fixationAndLeverTimes = [];
if isfield(D, 'adjDirects') && ~isempty(D.adjDirects)
    fprintf('Determining fixation and lever event times around each trial.\n');
    fixationAndLeverTimes = getFixationAndLeverTimes(D, cueOnset, firstJuiceEvent, cueLoc, isHoldTrial, nLoc);
end

%%
usefulEvents = var2struct(cueOnset, cueOnsetByLoc, ...
        cueOnsetRel, cueOnsetHold, cueOnsetRelByLoc, cueOnsetHoldByLoc, ...
        arrayOnset, arrayOnsetRel, arrayOnsetHold, arrayOnsetByLoc, ...
        arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
        arrayOnsetHoldBalShortHoldByLoc, arrayOnsetHoldBalLongHoldByLoc, ...
        targetDim, targetDimByLoc, ...
        targetDimBalShortHold, targetDimBalLongHold, ...
        targetDimBalShortHoldByLoc, targetDimBalLongHoldByLoc, ...
        nTrialShortHold, nTrialLongHold, cueLoc, cueLocHold, cueLocRel, ...
        cueTargetDelayDur, targetDimDelayDur, isHoldTrial, rt, ...
        fixationAndLeverTimes, firstJuiceEvent, targetDimMatch, ...
        arrayOnsetRelBal, arrayOnsetHoldBal, arrayOnsetRelBalByLoc, ...
        arrayOnsetHoldBalByLoc, targetDimBal, targetDimBalByLoc, ...
        isArrayRelBalanced, isArrayHoldBalanced);


