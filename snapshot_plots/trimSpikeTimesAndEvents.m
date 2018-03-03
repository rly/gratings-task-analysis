function D = trimSpikeTimesAndEvents(D, blockInds)

% determine which spike times to keep from each block, and logical OR them
% together
eventsToKeep = cell(size(D.events));
for i = 1:numel(D.events)
    eventsToKeep{i} = zeros(size(D.events{i}));
end
if isfield(D, 'allSpikeStructs')
    spikesToKeep = cell(size(D.allSpikeStructs));
    for i = 1:numel(D.allSpikeStructs)
        spikesToKeep{i} = zeros(size(D.allSpikeStructs{i}.ts));
    end
end
if isfield(D, 'allMUAStructs')
    muaToKeep = cell(size(D.allMUAStructs));
    for i = 1:numel(D.allMUAStructs)
        muaToKeep{i} = zeros(size(D.allMUAStructs{i}.ts));
    end
end
if isfield(D, 'adjLfps')
    lfpIndicesToKeep = false(1, size(D.adjLfps, 2));
end
for j = 1:numel(blockInds)
    blockStartTime = D.blockStartTimes(blockInds(j));
    blockStopTime = D.blockStopTimes(blockInds(j));
    fprintf('Block %d: Keeping event and spike times between %0.1f s and %0.1f s\n', blockInds(j), blockStartTime, blockStopTime);

    for i = 1:numel(D.events)
        eventsToKeep{i} = eventsToKeep{i} | (D.events{i} >= blockStartTime & D.events{i} <= blockStopTime);
    end
    if isfield(D, 'allSpikeStructs')
        for i = 1:numel(D.allSpikeStructs)
            spikesToKeep{i} = spikesToKeep{i} | (D.allSpikeStructs{i}.ts >= blockStartTime & D.allSpikeStructs{i}.ts <= blockStopTime);
        end
    end
    if isfield(D, 'allMUAStructs')
        for i = 1:numel(D.allMUAStructs)
            muaToKeep{i} = muaToKeep{i} | (D.allMUAStructs{i}.ts >= blockStartTime & D.allMUAStructs{i}.ts <= blockStopTime);
        end
    end
    if isfield(D, 'adjLfps')
        blockLfpIndices = max(1, floor(blockStartTime * D.lfpFs)):ceil(blockStopTime * D.lfpFs);
        lfpIndicesToKeep(blockLfpIndices) = true;
    end
end

% remove the non-marked events and spikes
for i = 1:numel(D.events)
    D.events{i}(~eventsToKeep{i},:) = [];
%     fprintf('EVT%02d: Removing %d events\n', i, sum(~eventsToKeep{i}));
end
if isfield(D, 'allSpikeStructs')
    for i = 1:numel(D.allSpikeStructs)
        D.allSpikeStructs{i}.wf(~spikesToKeep{i},:) = [];
        D.allSpikeStructs{i}.ts(~spikesToKeep{i}) = [];
        % TODO recompute mean and sd?
    end
end
if isfield(D, 'allMUAStructs')
    for i = 1:numel(D.allMUAStructs)
        D.allMUAStructs{i}.wf(~muaToKeep{i},:) = [];
        D.allMUAStructs{i}.ts(~muaToKeep{i}) = [];
        % TODO recompute mean and sd?
    end
end
if isfield(D, 'adjLfps')
    D.adjLfps(:,~lfpIndicesToKeep) = NaN;
    % TODO special case for adjLfps -- make the matrix smaller and adjust
    % the event times accordingly
end
