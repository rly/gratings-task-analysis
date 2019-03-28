function adjLfpsClean = interpolateLfpOverSpikeTimes(adjLfps, channelIDs, lfpFs, allSpikeStructs)

% remove 2 ms window around spike time (i.e. the data point before and 
% after since we're sampling at 1000 Hz) and linearly interpolate
fprintf('Interpolating over spike times...\n');
adjLfpsClean = adjLfps;
nChannels = size(adjLfps, 1);
nTime = size(adjLfps, 2);
assert(nChannels == numel(channelIDs));
for j = 1:nChannels
    isSpikeTimeThisCh = false(1, nTime);
    unitsThisCh = findAllUnitsOnCh(allSpikeStructs, channelIDs(j));
    for k = 1:numel(unitsThisCh)
        if floor(allSpikeStructs{unitsThisCh(k)}.ts(end) * lfpFs) > nTime
            error('Spikes extend past LFP array time');
        end
        isSpikeTimeThisCh(max(1, floor(allSpikeStructs{unitsThisCh(k)}.ts * lfpFs))) = true;
        isSpikeTimeThisCh(min(nTime, ceil(allSpikeStructs{unitsThisCh(k)}.ts * lfpFs))) = true;
    end
    ix = 1:size(adjLfps, 2);
    adjLfpsClean(j,isSpikeTimeThisCh) = interp1(ix(~isSpikeTimeThisCh), ...
            adjLfps(j,~isSpikeTimeThisCh), ix(isSpikeTimeThisCh), 'linear');
end