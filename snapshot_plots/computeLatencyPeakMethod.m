function latencyInfo = computeLatencyPeakMethod(spdf, eventT, ...
        eventWindow, eventTAnalysisLogical, meanBaseline, peakMin, ...
        isFindTrough)
% compute latency using the peak-based method min 1.5 SDs above baseline
% (computing SD based on variability of mean firing per time point)  (note
% that 2 SDs is standard. 1 SD has been used but with additional
% bootstrapping to check for robustness) where baseline for cue period is
% the pre-cue period and baseline for array period is the baseline for cue
% period or the pre-array period, whichever is smaller, to account for
% responses that are obvious after suppressed activity during the delay
% period min: 25% of greatest peak in peak range cue peak range: 25ms to
% 150ms array peak range: 25ms to 150ms minimum latency: 25ms, which
% doesn't matter much because I am reporting the median and comparing
% latency distributions in a non-parametric way. latency is defined as the
% time for the response to reach halfway between the peak firing rate and
% the baseline firing rate . the firing rate during the time period between
% the latency and the time of its associated peak cannot drop below 25% of
% that peak.

fprintf('10\n');
spdf
numel(eventT)
numel(eventTAnalysisLogical)
assert(all(~isnan(spdf)));
assert(numel(eventT) == numel(eventTAnalysisLogical));
fprintf('11\n');

if isFindTrough
    % reverse the sign
    spdf = -1*spdf;
    meanBaseline = -1*meanBaseline;
    peakMin = -1*peakMin;
end

fracPeak = 0.5;
minFracOfLargestPeak = 0.5;
minFracOfLargestPeakDrop = 0.25;

% peak is where the difference goes from + to -, with a time lag of 1
% expand eventTAnalysisLogical so that the first point is considered for
% having positive slope
indexBeforeEventTAnalysis = find(eventTAnalysisLogical, 1) - 1;
augEventTAnalysisLogical = eventTAnalysisLogical;
augEventTAnalysisLogical(indexBeforeEventTAnalysis) = 1;
absDiffSpdf = sign(diff(spdf(augEventTAnalysisLogical)));
% peak is where the difference of sign diff is -2, with an
% overall time lag of 2
% include last time point of window as possible peak
diffAbsDiffSpdf = [diff(absDiffSpdf) -2];
allPeakIndices = find(diffAbsDiffSpdf == -2 & ...
        spdf(eventTAnalysisLogical) >= meanBaseline) + ...
        find(eventTAnalysisLogical, 1) - 1;
allPeakRates = spdf(allPeakIndices);

% TODO discard peaks based on <= 3 spikes or based on <= 1 trial
% drop any peaks that are <15% magnitude of the largest peak
peakCutoff = (max(allPeakRates) - meanBaseline)*minFracOfLargestPeak + meanBaseline;
allPeakIndices(allPeakRates < peakCutoff) = [];
allPeakRates(allPeakRates < peakCutoff) = [];

% drop any peaks that are less than 0.5 SDs more than baseline
allPeakIndices(allPeakRates < peakMin) = [];
allPeakRates(allPeakRates < peakMin) = [];

if ~isempty(allPeakIndices)
    firstPeakIndex = allPeakIndices(1);
else
    firstPeakIndex = find(eventTAnalysisLogical, 1, 'last'); % temp
end

peakRate = spdf(firstPeakIndex);

fracPeakThreshold = (peakRate - meanBaseline)*fracPeak + meanBaseline;
indicesFracPeakUpward = find(spdf(eventTAnalysisLogical) >= fracPeakThreshold & ...
        absDiffSpdf == 1) + find(eventTAnalysisLogical, 1) - 1;
indicesFracPeakUpward(indicesFracPeakUpward > firstPeakIndex) = [];

% get earliest half-peak time where the rate between half-peak
% and peak doesn't drop below minFracOfLargestPeakDrop of peak
fracPeakMinDrop = (peakRate - meanBaseline)*minFracOfLargestPeakDrop + meanBaseline;
if ~isempty(indicesFracPeakUpward)
    firstIndexFracPeakUpward = indicesFracPeakUpward(1);
    while any(spdf(firstIndexFracPeakUpward:firstPeakIndex) < fracPeakMinDrop)
        indicesFracPeakUpward(1) = [];
        if isempty(indicesFracPeakUpward)
            firstIndexFracPeakUpward = [];
            break;
        else
            firstIndexFracPeakUpward = indicesFracPeakUpward(1);
        end
    end
end
if isempty(indicesFracPeakUpward) || isempty(firstIndexFracPeakUpward)
    firstIndexFracPeakUpward = find(eventTAnalysisLogical, 1, 'last'); % temp;
end

timeEventToPeak = eventT(firstPeakIndex) - eventWindow(1);

if peakRate <= peakMin
    latency = NaN;
    latencyRate = NaN;
    timeEventToPeak = NaN;
    peakRate = NaN;
    latencyTInd = NaN;
    peakTInd = NaN;
else
    latency = eventT(firstIndexFracPeakUpward) - eventWindow(1);
    latencyRate = spdf(firstIndexFracPeakUpward);
    latencyTInd = firstIndexFracPeakUpward;
    peakTInd = firstPeakIndex;
end

if isFindTrough
    % reverse the sign
    latencyRate = -1*latencyRate;
    peakRate = -1*peakRate;
end

latencyInfo = var2struct(latency, latencyTInd, latencyRate, timeEventToPeak, ...
        peakTInd, peakRate, isFindTrough);

