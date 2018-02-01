function [aligned,t,sampleInfo,offset] = sliceLfpAroundTrials(data, cueOnset, firstJuice, Fs, ...
        timeFromCueOnset, timeFromFirstJuice)
% data: nChannel x nSamples

assert(numel(cueOnset) == numel(firstJuice));
nTrial = numel(cueOnset);

assert(ismatrix(data));

startTimes = round((cueOnset + timeFromCueOnset) * Fs);
endTimes = round((firstJuice + timeFromFirstJuice) * Fs) - 1;
offset = round(timeFromCueOnset * Fs);

aligned = cell(1, nTrial); % 1 x nTrial
t = cell(1, nTrial); % 1 x nTrial
tInd = cell(1, nTrial); % 1 x nTrial

for i = 1:nTrial
    tInd{i} = startTimes(i):endTimes(i);
    t{i} = ((1:numel(tInd{i})) - 1 + round(timeFromCueOnset * Fs)) / Fs; %tInd{i} / Fs;
    aligned{i} = data(:,tInd{i});
    assert(~any(any(isnan(aligned{i}))));
end

sampleInfo = [startTimes endTimes];
