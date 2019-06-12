function [sfcLF, fAxisLF, sfcHF, fAxisHF] = computeSfcWrapper(spikeTs, eventTimes, ...
        eventOffset, alignedLfps, paramsLF, paramsHF, adjCohNormRate)

nTrial = numel(eventTimes);
assert(nTrial == size(alignedLfps, 2));

fAxisLF = NaN;
fAxisHF = NaN;
sfcLF = NaN;
sfcHF = NaN;

alignedSpikeTs = createnonemptydatamatpt(spikeTs, eventTimes, eventOffset .* [-1 1]);
meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, eventOffset);
fprintf('\tMean FR: %0.2f Hz\n', meanFR);

if meanFR < 0.01
    fprintf('\tMean spike rate (%0.2f Hz) is too low; skipping...\n', meanFR); 
    return; % skip this channel
end

% high frequency range uses different number of tapers
[C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(alignedLfps, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
if any(C(:) > 0.8) % coherence above 0.8 is unusual and also atanh() values close to 1 goes to Inf
    fprintf('\tAbnormally high coherence; skipping...\n'); 
    return; % skip this channel
end
sfcLF = atanh(C)-(1/((2*paramsLF.tapers(2)*nTrial)-2)); % adjust for num trials

% high frequency range uses different number of tapers
[C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(alignedLfps, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
if any(C(:) > 0.8) % coherence above 0.8 is unusual and also atanh() values close to 1 goes to Inf
    fprintf('\tAbnormally high coherence; skipping...\n'); 
    return; % skip this channel
end
sfcHF = atanh(C)-(1/((2*paramsHF.tapers(2)*nTrial)-2)); % adjust for num trials
