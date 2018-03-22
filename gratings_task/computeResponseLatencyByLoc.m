function eventStruct = computeResponseLatencyByLoc(eventStruct, isLocUsed)

eventPsthT = eventStruct.t;
baselineWindowOffset = [-0.3 0];
baselineWindow = eventStruct.window(1) + baselineWindowOffset; % from array onset
baselineIndices = getTimeLogicalWithTolerance(eventPsthT, baselineWindow);
latencyWindowOffset = [0.025 0.125];
latencyWindow = eventStruct.window(1) + latencyWindowOffset;
latencyIndices = getTimeLogicalWithTolerance(eventPsthT, latencyWindow);
minPeakForLatency = 3; % SDs from mean
maxTroughForLatency = -minPeakForLatency;
nLoc = numel(isLocUsed);

%% use time to half peak method
eventStruct.latencyInfoByLoc = cell(nLoc, 1);
eventStruct.latencyByLoc = nan(nLoc, 1);
for k = 1:nLoc
    if ~isLocUsed(k)
        continue;
    end
    psthResponse = eventStruct.spdfByLoc(k,:);
    baselineResponseByPsth = mean(psthResponse(baselineIndices));
    baselineResponseSDOverTimeByPsth = std(psthResponse(baselineIndices));
    normPsthResponse = (psthResponse - baselineResponseByPsth) / baselineResponseSDOverTimeByPsth;
    % check peak
    peakLatencyInfo = computeLatencyPeakMethod(normPsthResponse, eventPsthT, ...
            eventStruct.window, latencyIndices, 0, minPeakForLatency, 0); 
    % check trough
    troughLatencyInfo = computeLatencyPeakMethod(normPsthResponse, eventPsthT, ...
            eventStruct.window, latencyIndices, 0, maxTroughForLatency, 1);

    % get the earlier one
    eventStruct.latencyInfoByLoc{k} = getEarlierPeakTroughLatencyInfo(peakLatencyInfo, troughLatencyInfo);
    eventStruct.latencyByLoc(k) = eventStruct.latencyInfoByLoc{k}.latency;
end

%% use time to significant deviation from resampled baseline method

alpha = 0.05;
kernelSigma = 0.01;
numRandomizations = 25; % few randomizations are needed because we look at each sample in the 300 ms baseline independently. also it is very slow.
minNumSigSamples = 10;
eventStruct.latencyBootInfoByLoc = cell(nLoc, 1);
eventStruct.latencyBootByLoc = nan(nLoc, 1);
for k = 1:nLoc
    if ~isLocUsed(k)
        continue;
    end
    nTrials = numel(eventStruct.spikeTimesByLoc{k});
    resampledBaselines = zeros(numRandomizations, sum(baselineIndices));
    for m = 1:numRandomizations
        bootRandInd = randi(nTrials, 1, nTrials);
        bootstrapPsth = fixedPsth(eventStruct.spikeTimesByLoc{k}(bootRandInd), kernelSigma, 0, eventPsthT);
        resampledBaselines(m,:) = bootstrapPsth(baselineIndices);
    end
    fprintf('.');

    upperThresh = prctile(resampledBaselines(:), (1 - alpha/2) * 100);
    lowerThresh = prctile(resampledBaselines(:), alpha/2 * 100);

    isSignificant = (eventStruct.spdfByLoc(k,:) > upperThresh | eventStruct.spdfByLoc(k,:) < lowerThresh);
    crossings = find(isSignificant & latencyIndices);
    latency = NaN;
    latencyRate = NaN;
    latencyTInd = NaN;
    for i = 1:numel(crossings)
        if all(ismember(crossings(i):(crossings(i) + minNumSigSamples - 1), crossings))
            latencyTInd = crossings(i);
            latency = eventPsthT(latencyTInd) - eventStruct.window(1);
            latencyRate = eventStruct.spdfByLoc(k,latencyTInd);
            break;
        end
    end
    eventStruct.latencyBootInfoByLoc{k} = var2struct(latency, latencyTInd, latencyRate);
    eventStruct.latencyBootByLoc(k) = latency;
end

