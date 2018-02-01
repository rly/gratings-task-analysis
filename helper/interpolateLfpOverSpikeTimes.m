function adjLfpsClean = interpolateLfpOverSpikeTimes(adjLfps, lfpFs, allSpikeStructs)

% % remove 2 ms window around spike time (i.e. the data point before and 
% % after since we're sampling at 1000 Hz) and interpolate
fprintf('Interpolating over spike times...\n');
adjLfpsClean = adjLfps;
nChannels = size(adjLfps, 1);
for j = 1:nChannels
    isSpikeTimeThisCh = false(size(adjLfps(j,:)));
    unitsThisCh = findAllUnitsOnCh(allSpikeStructs, j);
    for k = 1:numel(unitsThisCh)
        isSpikeTimeThisCh(floor(allSpikeStructs{k}.ts * lfpFs)) = true;
        isSpikeTimeThisCh(ceil(allSpikeStructs{k}.ts * lfpFs)) = true;
    end
    ix = 1:size(adjLfps, 2);
    adjLfpsClean(j,isSpikeTimeThisCh) = interp1(ix(~isSpikeTimeThisCh), adjLfps(j,~isSpikeTimeThisCh), ix(isSpikeTimeThisCh));
end