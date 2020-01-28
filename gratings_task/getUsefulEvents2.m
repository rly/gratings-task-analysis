function usefulEvents = getUsefulEvents2(logDir, logIndices, nLoc, D, blockName)

holdDurMid = 850; % ms; for splitting short vs long hold

% cueP3Event = D.events{2};
cueOnsetEvent = D.events{4};
arrayOnsetHoldEvent = D.events{6}; % TODO
arrayOnsetReleaseEvent = D.events{7};
% firstFixationEvent = D.events{13}; % TODO use for assertions
% arrayOnsetHoldP1Event = D.events{14};
% arrayOnsetReleaseP1Event = D.events{15};
targetDimEvent = D.events{16};
juiceEvent = D.events{8};

%% load presentation log information on correct trials only
logRoot = 'gratings_task_eye_tracking_';
trialParamsAllCorrect = readPresentationLogCorrectTrial(logDir, logRoot, logIndices);

%% load more detailed trial info from json -> mat
% use saveGratingsTaskResultsJsonToMatRunner.m to generate .mat files
load(sprintf('%s/gratingsTaskResultsJson-%s.mat', logDir, blockName), 'trialStructs', 'arrayShapes');
trialStructs = trialStructs'; % make column
arrayShapes = arrayShapes'; % make column
trialResults = cellfun(@(x) x.response.type, trialStructs, 'UniformOutput', false);
isCorrect = strcmp(trialResults, 'correct-response');
trialStructsCorrect = trialStructs(isCorrect);
verifyGratingsTaskLogsAndJsonCorrect(trialParamsAllCorrect, trialStructsCorrect);
arrayShapesCorrect = arrayShapes(isCorrect);

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
trialStructsCorrect = trialStructsCorrect(trialsToKeep);
arrayShapesCorrect = arrayShapesCorrect(trialsToKeep);
fprintf('Removed %d/%d putative hand-started trials.\n', sum(~trialsToKeep), numel(trialsToKeep));

%%
cueLoc = trialParamsAllCorrect(:,4);
isHoldTrial = trialParamsAllCorrect(:,7) ~= -1;

% fixationAndLeverTimes = [];
% if isfield(D, 'adjDirects') && ~isempty(D.adjDirects)
    fprintf('Determining fixation and lever event times around each trial.\n');
    [fixationAndLeverTimes,isMissingData] = getFixationAndLeverTimes(D, cueOnset, firstJuiceEvent, cueLoc, isHoldTrial, nLoc);
% end
cueOnset = cueOnset(~isMissingData);
firstJuiceEvent = firstJuiceEvent(~isMissingData);
trialParamsAllCorrect = trialParamsAllCorrect(~isMissingData,:);
trialStructsCorrect = trialStructsCorrect(~isMissingData);
arrayShapesCorrect = arrayShapesCorrect(~isMissingData);
fprintf('Removed %d/%d trials with missing matching lever/fixation data.\n', sum(isMissingData), numel(isMissingData));

%%
isRelBal = cellfun(@(x) x(1) == 'R' && x(3) == 'R', arrayShapesCorrect);
isHoldBal = cellfun(@(x) x(1) == 'H' && x(3) == 'H', arrayShapesCorrect);

isCorrectRelTrialLog = cellfun(@(x,y) y(1) == 'R', trialStructsCorrect, arrayShapesCorrect);
isCorrectHoldTrialLog = cellfun(@(x,y) y(1) == 'H', trialStructsCorrect, arrayShapesCorrect);

%%
% TODO get only trials that are not repeats
% notRepeatLogical = trialParamsAllCorrect(:,3) == 1;
% juiceEvent = juiceEvent(notRepeatLogical,:);
trialParams = trialParamsAllCorrect;
% numNonRptCorrectTrials = size(trialParams, 1);
cueLoc = trialParams(:,4);
cueTargetDelayDur = trialParams(:,6);
targetDimDelayDur = trialParams(:,7);
% these are the reported delay durations, but may not be the actual delay
% durations
assert(all(cellfun(@(x) x.cueLoc, trialStructsCorrect) == cueLoc));
assert(all(cellfun(@(x) x.cueArrayDelayDuration, trialStructsCorrect) == cueTargetDelayDur));
assert(all(cellfun(@(x) x.holdShapeDuration .* (x.holdShapeDuration ~= 15000) + -1 * (x.holdShapeDuration == 15000), ...
        trialStructsCorrect) == targetDimDelayDur));

%% get cue onset events corresponding to correct trials only
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
% maxArrayOnsetToJuiceTimeReleaseShape = 0.925;
% isHoldTrial = firstJuiceEvent - arrayOnset >= maxArrayOnsetToJuiceTimeReleaseShape;
isHoldTrial = trialParams(:,7) ~= -1;
% assert(all(isHoldTrial == (trialParams(:,7) ~= -1)));
assert(all(isHoldTrial == isCorrectHoldTrialLog));
assert(all(~isHoldTrial == isCorrectRelTrialLog));
% arrayOnsetRel = arrayOnset(~isHoldTrial);
% arrayOnsetHold = arrayOnset(isHoldTrial);

arrayOnsetRelBal = arrayOnset(~isHoldTrial & isRelBal);
arrayOnsetHoldBal = arrayOnset(isHoldTrial & isHoldBal);

arrayOnsetByLoc = cell(nLoc, 1);
% arrayOnsetRelByLoc = cell(nLoc, 1);
% arrayOnsetHoldByLoc = cell(nLoc, 1);
% arrayOnsetBalByLoc = cell(nLoc, 1);
arrayOnsetRelBalByLoc = cell(nLoc, 1);
arrayOnsetHoldBalByLoc = cell(nLoc, 1);
% arrayOnsetHoldShortHoldBalByLoc = cell(nLoc, 1);
% arrayOnsetHoldLongHoldBalByLoc = cell(nLoc, 1);

for i = 1:nLoc
    arrayOnsetByLoc{i} = arrayOnset(cueLoc == i);
%     arrayOnsetRelByLoc{i} = arrayOnset(~isHoldTrial & cueLoc == i);
%     arrayOnsetHoldByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i);
%     arrayOnsetBalByLoc{i} = arrayOnset(cueLoc == i & (isArrayRelBalanced | isArrayHoldBalanced));
    arrayOnsetRelBalByLoc{i} = arrayOnset(~isHoldTrial & cueLoc == i & isRelBal);
    arrayOnsetHoldBalByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & isHoldBal);
%     arrayOnsetHoldShortHoldBalByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & targetDimDelayDur < holdDurMid & isArrayHoldBalanced);
%     arrayOnsetHoldLongHoldBalByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & targetDimDelayDur >= holdDurMid & isArrayHoldBalanced);
end

%% split cue events and info by whether the trial is a hold trial or not
cueLocRelBal = cueLoc(~isHoldTrial & isRelBal);
cueLocHoldBal = cueLoc(isHoldTrial & isHoldBal);

cueOnsetRelBal = cueOnset(~isHoldTrial & isRelBal);
cueOnsetHoldBal = cueOnset(isHoldTrial & isHoldBal);

cueOnsetRelBalByLoc = cell(nLoc, 1);
cueOnsetHoldBalByLoc = cell(nLoc, 1);
for i = 1:nLoc
    cueOnsetRelBalByLoc{i} = cueOnset(~isHoldTrial & cueLoc == i & isRelBal);
    cueOnsetHoldBalByLoc{i} = cueOnset(~isHoldTrial & cueLoc == i & isHoldBal);
end 

%% get target dim events corresponding to correct trials only
maxTargetDimToJuiceTime = 1; % seconds
targetDimMatch = nan(numel(firstJuiceEvent), 1);
for i = 1:numel(firstJuiceEvent)  
    prevEvent7Ind = find((firstJuiceEvent(i) - targetDimEvent < maxTargetDimToJuiceTime) & ...
            (firstJuiceEvent(i) - targetDimEvent > 0), 1, 'last');
    if ~isempty(prevEvent7Ind)
        targetDimMatch(i) = targetDimEvent(prevEvent7Ind);
    end
end

targetDim = targetDimMatch(~isnan(targetDimMatch));
targetDimBal = targetDimMatch(~isnan(targetDimMatch) & isHoldBal);
targetDimShortHoldBal = targetDimMatch(~isnan(targetDimMatch) & targetDimDelayDur < holdDurMid & isHoldBal);
targetDimLongHoldBal = targetDimMatch(~isnan(targetDimMatch) & targetDimDelayDur >= holdDurMid & isHoldBal);
nTrialShortHold = sum(~isnan(targetDimMatch) & targetDimDelayDur < holdDurMid & isHoldBal);
nTrialLongHold = sum(~isnan(targetDimMatch) & targetDimDelayDur >= holdDurMid & isHoldBal);
assert(numel(unique(targetDimBal)) == numel(targetDimBal));

targetDimByLoc = cell(nLoc, 1);
targetDimBalByLoc = cell(nLoc, 1);
targetDimShortHoldBalByLoc = cell(nLoc, 1);
targetDimLongHoldBalByLoc = cell(nLoc, 1);
for i = 1:nLoc
    targetDimByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & cueLoc == i);
    targetDimBalByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & cueLoc == i & isHoldBal);
    targetDimShortHoldBalByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & ...
            cueLoc == i & targetDimDelayDur < holdDurMid & isHoldBal);
    targetDimLongHoldBalByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & ...
            cueLoc == i & targetDimDelayDur >= holdDurMid & isHoldBal);
end

%% check rt match across logs and event timing
% four ways to compute RT
% 1) rtPresLogs: time between correct event and array onset or target
% dimming event in presentation logs
% 2) rtJuice: juice digital event minus array onset or target dimming event
% -- this is discretized in 12 ms bins
% 3) rtJsonLogs: reported RT from array onset in json file (using hold
% shape duration reported in json file)
% 4) rtLever: lever release minus array onset or target dimming event
% 5) rtLogEvent: time from array onset (where hold shape duration is
% defined as time between target dimming and array onset events)
rtPresLogs = trialParams(:,8);
rtJuice = nan(numel(trialStructsCorrect), 1);
rtJsonLogs = nan(numel(trialStructsCorrect), 1);
rtLever = nan(numel(trialStructsCorrect), 1);
rtLogEvent = nan(numel(trialStructsCorrect), 1);
for i = 1:numel(trialStructsCorrect)
    if isCorrectRelTrialLog(i)
        rtJuice(i) = firstJuiceEvent(i) - arrayOnset(i);
        rtJsonLogs(i) = trialStructsCorrect{i}.response.responseTime / 1000;
        rtLever(i) = fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(i) - arrayOnset(i);
        rtLogEvent(i) = trialStructsCorrect{i}.response.responseTime / 1000;
        % fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(i) is
        % based on analog signal
        % arrayOnset is based on digital event signal
    else
        rtJuice(i) = firstJuiceEvent(i) - targetDimMatch(i);
        rtJsonLogs(i) = (trialStructsCorrect{i}.response.responseTime - trialStructsCorrect{i}.holdShapeDuration) / 1000;
        rtLever(i) = fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(i) - targetDimMatch(i);
        rtLogEvent(i) = trialStructsCorrect{i}.response.responseTime / 1000 - (targetDimMatch(i) - arrayOnset(i));
    end
end
rt = rtPresLogs;
 
% mean(rtJuice - rtPresLogs) % offset by mean 8 ms, sd 4 ms
% mean(rtJuice - rtJsonLogs) % offset by mean 9 ms, sd 5 ms
% mean(rtJuice - rtLever) % offset by mean 37 ms, sd 4 ms
% mean(rtJuice - rtLogEvent) % offset by mean 8 ms, sd 5 ms
% mean(rtPresLogs - rtLever) % offset by mean 30 ms, sd 3 ms

%% process error events
trialResponseState = cellfun(@(x) x.response.trialState, trialStructs);
isEarlyReleaseHoldAll = strcmp(trialResults, 'array-on-early-release');
assert(sum(isEarlyReleaseHoldAll) == sum(isEarlyReleaseHoldAll & trialResponseState == -108)); % ARRAY_ON state

isLateReleaseAll = strcmp(trialResults, 'late-response');
assert(sum(isLateReleaseAll) == sum((isLateReleaseAll & trialResponseState == -113) | ... % PAST_RESPONSE_PERIOD state
        (isLateReleaseAll & trialResponseState == -118) | ... % PAST_RESPONSE_PERIOD_RELEASE_SHAPE state
        (isLateReleaseAll & trialResponseState == -114))); % RELEASE_SHAPE_DIMMED
isLateReleaseRelAll = isLateReleaseAll & (trialResponseState == -118 | trialResponseState == -114);

isEyeErrorRel = strcmp(trialResults, 'release-shape-on-eye-error');
assert(sum(isEyeErrorRel) == sum((isEyeErrorRel & trialResponseState == -109) | ... % ARRAY_ON_RELEASE_SHAPE state
        (isEyeErrorRel & trialResponseState == -117) | ... % WAITING_FOR_RESPONSE_RELEASE_SHAPE state
        (isEyeErrorRel & trialResponseState == -118) | ... % PAST_RESPONSE_PERIOD_RELEASE_SHAPE state
        (isEyeErrorRel & trialResponseState == -114))); % RELEASE_SHAPE_DIMMED
isLateSaccadeRelAll = isEyeErrorRel & (trialResponseState == -118 | trialResponseState == -114);

isEyeErrorHold = strcmp(trialResults, 'hold-shape-on-eye-error');
assert(sum(isEyeErrorHold) == sum((isEyeErrorHold & trialResponseState == -108) | ... % ARRAY_ON state
        (isEyeErrorHold & trialResponseState == -110) | ... % TARGET_DIMMED state
        (isEyeErrorHold & trialResponseState == -112) | ... % WAITING_FOR_RESPONSE state
        (isEyeErrorHold & trialResponseState == -113))); % PAST_RESPONSE_PERIOD
isEarlySaccadeHoldAll = isEyeErrorHold & trialResponseState == -104;

% only process these errors for now
isProcessedError = isEarlyReleaseHoldAll | isLateReleaseRelAll | isLateSaccadeRelAll | isEarlySaccadeHoldAll;

% get new indices
temp = zeros(size(isProcessedError));
temp(isEarlyReleaseHoldAll | isEarlySaccadeHoldAll) = 1;
temp(isLateReleaseRelAll | isLateSaccadeRelAll) = 2;
temp(temp == 0) = [];
isErrorEarlyResponseHold = temp == 1;
isErrorLateResponseRel = temp == 2;
clear temp;
fprintf('Found %d early response errors (release and eye) on hold trials.\n', sum(isErrorEarlyResponseHold));
fprintf('Found %d late response errors (release and eye) on release trials.\n', sum(isErrorLateResponseRel));

% int EYE_ERROR = 1; # event 9 to port 1
% int RELEASE_ERROR_PRE_ARRAY_ON = 2; # event 10 to port 1
% int RELEASE_ERROR_ARRAY_ON = 3; # event 9,10 to port 1
% int RELEASE_ERROR_ARRAY_ON_AFTER_DCHANGE = 4; # event 11 to port 1
% int RELEASE_ERROR_JUST_EARLY = 5; # event 9,11 to port 1
% int RELEASE_ERROR_LATE = 6; # event 10,11 to port 1
eventTol = 0.001;
errorTimesAll = sort([D.events{9}; D.events{10}; D.events{11}; D.events{12}]);
toRemove = diff(errorTimesAll) < eventTol;
toRemove = [false; toRemove];
errorTimesAll(toRemove) = [];
eventCode = zeros(size(errorTimesAll));
for i = 1:numel(errorTimesAll)
    if any(abs(D.events{9} - errorTimesAll(i)) < eventTol)
        eventCode(i) = eventCode(i) + 1;
    end
    if any(abs(D.events{10} - errorTimesAll(i)) < eventTol)
        eventCode(i) = eventCode(i) + 2;
    end
    if any(abs(D.events{11} - errorTimesAll(i)) < eventTol)
        eventCode(i) = eventCode(i) + 4;
    end
    if any(abs(D.events{12} - errorTimesAll(i)) < eventTol)
        eventCode(i) = eventCode(i) + 8;
    end
end
assert(all(eventCode > 0))
errorTimesAll(eventCode >= 8) = [];
eventCode(eventCode >= 8) = [];
assert(sum(~isCorrect) == numel(errorTimesAll));
assert(sum(isEarlyReleaseHoldAll) == sum(eventCode == 3));
assert(sum(isLateReleaseAll) == sum(eventCode == 6));
errorTimesProcessed = errorTimesAll(isProcessedError(~isCorrect));

trialStructsError = trialStructs(isProcessedError);
arrayShapesError = arrayShapes(isProcessedError);
cueLocError = cellfun(@(x) x.cueLoc, trialStructsError);
cueTargetDelayDurError = cellfun(@(x) x.cueArrayDelayDuration, trialStructsError);

%% get cue onsets corresponding to error trials only
maxCueOnsetToErrorTime = 5; % seconds
cueOnsetError = nan(numel(errorTimesProcessed), 1);
for i = 1:numel(errorTimesProcessed)  
    prevCueOnsetInd = find((errorTimesProcessed(i) - cueOnsetEvent < maxCueOnsetToErrorTime) & ...
            (errorTimesProcessed(i) - cueOnsetEvent > 0), 1, 'last');
    if ~isempty(prevCueOnsetInd)
        cueOnsetError(i) = cueOnsetEvent(prevCueOnsetInd);
    else
        warning('Error event without a corresponding cue onset event: %d, %f', i, errorTimesProcessed(i));
    end
end

%% remove trials that do not have an earlier lever press
% [cueOnset,trialsToKeep] = removeHandStartedTrials(D, cueOnset, firstJuiceEvent);
% firstJuiceEvent = firstJuiceEvent(trialsToKeep);
% trialParamsAllCorrect = trialParamsAllCorrect(trialsToKeep,:);
% fprintf('Removed %d/%d putative hand-started trials.\n', sum(~trialsToKeep), numel(trialsToKeep));
% 
% trialStructsCorrect = trialStructsCorrect(trialsToKeep);
% arrayShapesCorrect = arrayShapesCorrect(trialsToKeep);

%%
isRelBalError = cellfun(@(x) x(1) == 'R' && x(3) == 'R', arrayShapesError);
isHoldBalError = cellfun(@(x) x(1) == 'H' && x(3) == 'H', arrayShapesError);

% isErrorRelTrialLog = cellfun(@(x,y) y(1) == 'R', trialStructsError, arrayShapesError);
isErrorHoldTrialLog = cellfun(@(x,y) y(1) == 'H', trialStructsError, arrayShapesError);
assert(all(isErrorHoldTrialLog == isErrorEarlyResponseHold));

%% get cue onset events corresponding to correct trials only
cueOnsetErrorByLoc = cell(nLoc, 1);
for i = 1:nLoc
    cueOnsetErrorByLoc{i} = cueOnsetError(cueLocError == i);
end    

%% get array onsets corresponding to correct trials only
maxArrayOnsetToJuiceTime = 2.5; % seconds
arrayOnsetError = nan(numel(errorTimesProcessed), 1);
allArrayOnsetEvents = [arrayOnsetHoldEvent; arrayOnsetReleaseEvent];
for i = 1:numel(errorTimesProcessed) 
    prevArrayOnsetInd = find((errorTimesProcessed(i) - allArrayOnsetEvents < maxArrayOnsetToJuiceTime) & ...
            (errorTimesProcessed(i) - allArrayOnsetEvents > 0), 1, 'last');
    if ~isempty(prevArrayOnsetInd)
        arrayOnsetError(i) = allArrayOnsetEvents(prevArrayOnsetInd);
    else
        warning('Error event without a corresponding array onset event: %d, %f', i, errorTimesProcessed(i));
    end
end
isHoldTrialError = isErrorHoldTrialLog;

% arrayOnsetRel = arrayOnset(~isHoldTrial);
% arrayOnsetHold = arrayOnset(isHoldTrial);

arrayOnsetRelBalError = arrayOnsetError(~isHoldTrialError & isRelBalError);
arrayOnsetHoldBalError = arrayOnsetError(isHoldTrialError & isHoldBalError);

arrayOnsetErrorByLoc = cell(nLoc, 1);
% arrayOnsetRelByLoc = cell(nLoc, 1);
% arrayOnsetHoldByLoc = cell(nLoc, 1);
% arrayOnsetBalByLoc = cell(nLoc, 1);
arrayOnsetRelBalErrorByLoc = cell(nLoc, 1);
arrayOnsetHoldBalErrorByLoc = cell(nLoc, 1);
% arrayOnsetHoldShortHoldBalByLoc = cell(nLoc, 1);
% arrayOnsetHoldLongHoldBalByLoc = cell(nLoc, 1);

for i = 1:nLoc
    arrayOnsetErrorByLoc{i} = arrayOnsetError(cueLocError == i);
%     arrayOnsetRelByLoc{i} = arrayOnset(~isHoldTrial & cueLoc == i);
%     arrayOnsetHoldByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i);
%     arrayOnsetBalByLoc{i} = arrayOnset(cueLoc == i & (isArrayRelBalanced | isArrayHoldBalanced));
    arrayOnsetRelBalErrorByLoc{i} = arrayOnsetError(~isHoldTrialError & cueLocError == i & isRelBalError);
    arrayOnsetHoldBalErrorByLoc{i} = arrayOnsetError(isHoldTrialError & cueLocError == i & isHoldBalError);
%     arrayOnsetHoldShortHoldBalByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & targetDimDelayDur < holdDurMid & isArrayHoldBalanced);
%     arrayOnsetHoldLongHoldBalByLoc{i} = arrayOnset(isHoldTrial & cueLoc == i & targetDimDelayDur >= holdDurMid & isArrayHoldBalanced);
end

%% split cue events and info by whether the trial is a hold trial or not
cueLocRelBalError = cueLocError(~isHoldTrialError & isRelBalError);
cueLocHoldBalError = cueLocError(isHoldTrialError & isHoldBalError);

cueOnsetRelBalError = cueOnsetError(~isHoldTrialError & isRelBalError);
cueOnsetHoldBalError = cueOnsetError(isHoldTrialError & isHoldBalError);

cueOnsetRelBalErrorByLoc = cell(nLoc, 1);
cueOnsetHoldBalErrorByLoc = cell(nLoc, 1);
for i = 1:nLoc
    cueOnsetRelBalErrorByLoc{i} = cueOnsetError(~isHoldTrialError & cueLocError == i & isRelBalError);
    cueOnsetHoldBalErrorByLoc{i} = cueOnsetError(~isHoldTrialError & cueLocError == i & isHoldBalError);
end 

%% TODO fixation and lever times and RT for errors


%%
usefulEvents = var2struct(...
        cueOnset, ...
        cueOnsetByLoc, ...
        cueOnsetRelBal, cueOnsetHoldBal, ...
        cueOnsetRelBalByLoc, cueOnsetHoldBalByLoc, ...
        arrayOnset, ...
        arrayOnsetByLoc, ...
        arrayOnsetRelBal, arrayOnsetHoldBal, ...
        arrayOnsetRelBalByLoc, arrayOnsetHoldBalByLoc, ...
        ...arrayOnsetHoldShortHoldBalByLoc, arrayOnsetHoldLongHoldBalByLoc, ...
        targetDim, targetDimByLoc, ...
        targetDimBal, ...
        targetDimBalByLoc, ...
        targetDimShortHoldBal, targetDimLongHoldBal, ...
        targetDimShortHoldBalByLoc, targetDimLongHoldBalByLoc, ...
        nTrialShortHold, nTrialLongHold, ...
        cueLoc, cueLocHoldBal, cueLocRelBal, ...
        cueTargetDelayDur, targetDimDelayDur, rt, ...
        fixationAndLeverTimes, firstJuiceEvent, targetDimMatch, ...
        isHoldTrial, isRelBal, isHoldBal, ...
        cueOnsetError, ...
        cueOnsetErrorByLoc, ...
        cueOnsetRelBalError, cueOnsetHoldBalError, ...
        cueOnsetRelBalErrorByLoc, cueOnsetHoldBalErrorByLoc, ...
        arrayOnsetError, ...
        arrayOnsetErrorByLoc, ...
        arrayOnsetRelBalError, arrayOnsetHoldBalError, ...
        arrayOnsetRelBalErrorByLoc, arrayOnsetHoldBalErrorByLoc, ...
        cueLocError, cueLocHoldBalError, cueLocRelBalError, ...
        cueTargetDelayDurError, ...
        isHoldTrialError, isRelBalError, isHoldBalError, ...
        isErrorEarlyResponseHold, isErrorLateResponseRel, ...
        arrayShapes, arrayShapesCorrect, arrayShapesError);


