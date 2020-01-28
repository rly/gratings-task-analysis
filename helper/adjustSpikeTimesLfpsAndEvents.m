function D = adjustSpikeTimesLfpsAndEvents(D, blockInds)
% D events, spike structs, and LFPs have already been trimmed

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

if isfield(D, 'allUnitStructs')
    for i = 1:numel(D.allUnitStructs)
        D.allUnitStructs{i}.ts = D.allUnitStructs{i}.ts - firstBlockStartTime;
    end
end

if isfield(D, 'allMUAStructs')
    for i = 1:numel(D.allMUAStructs)
        D.allMUAStructs{i}.ts = D.allMUAStructs{i}.ts - firstBlockStartTime;
    end
end

if isfield(D, 'allUnitStructs')
    for i = 1:numel(D.allUnitStructs)
        D.allUnitStructs{i}.ts = D.allUnitStructs{i}.ts - firstBlockStartTime;
    end
end

if isfield(D, 'adjLfps')
    D.adjLfps(:,ceil(lastBlockStopTime * D.lfpFs):end) = [];
    D.adjLfps(:,1:floor(firstBlockStartTime * D.lfpFs)) = [];
end

if isfield(D, 'adjDirects')
    D.adjDirects(:,ceil(lastBlockStopTime * D.directFs):end) = [];
    D.adjDirects(:,1:floor(firstBlockStartTime * D.directFs)) = [];
end
