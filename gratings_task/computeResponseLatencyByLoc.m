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