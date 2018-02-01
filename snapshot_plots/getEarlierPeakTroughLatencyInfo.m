function latencyInfo = getEarlierPeakTroughLatencyInfo(peakLatencyInfo, troughLatencyInfo)

if isnan(peakLatencyInfo.latency)
    % peak doesn't exist, use trough latency even if it doesn't exist
    latencyInfo = troughLatencyInfo;
elseif ~isnan(troughLatencyInfo.latency)
    % both peak and trough exist
    if peakLatencyInfo.latency < troughLatencyInfo.latency
        latencyInfo = peakLatencyInfo;
    else
        latencyInfo = troughLatencyInfo;
    end
else
    % peak exists but trough doesn't
    latencyInfo = peakLatencyInfo;
end