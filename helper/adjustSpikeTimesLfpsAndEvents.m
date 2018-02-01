function D = adjustSpikeTimesLfpsAndEvents(D, blockInds)

firstBlockStartTime = D.blockStartTimes(blockInds(1));
lastBlockStopTime = D.blockStopTimes(blockInds(end));

for i = 1:numel(D.events)
    D.events{i} = D.events{i} - firstBlockStartTime;
end

if isfield(D, 'allSpikeStructs')
    for i = 1:numel(D.allSpikeStructs)
        D.allSpikeStructs{i}.ts = D.allSpikeStructs{i}.ts - firstBlockStartTime;
    end
end

if isfield(D, 'adjLfps')
    D.adjLfps(:,ceil(lastBlockStopTime * D.lfpFs):end) = [];
    D.adjLfps(:,1:floor(firstBlockStartTime * D.lfpFs)) = [];
end
