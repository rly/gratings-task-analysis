function meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, alignWindowOffset)

assert(isrow(alignedSpikeTs) || iscolumn(alignedSpikeTs));
assert(isfield(alignedSpikeTs, 'times'));
nTrial = numel(alignedSpikeTs);

nTime = diff(alignWindowOffset); % in seconds
assert(nTime > 0);

nSpike = sum(cellfun(@numel, {alignedSpikeTs.times}));

meanFR = nSpike / nTrial / nTime; % in Hz