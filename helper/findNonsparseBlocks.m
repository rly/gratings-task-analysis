function isFiringNonsparseByBlock = findNonsparseBlocks(D, allSpikeStructs, blockInds)
% compute whether firing is sparse for each block x unit
% 0 = sparse firing for that block, 1 = active firing
% allSpikeStructs can be single units or multi units

minFiringRate = 0.1; % Hz

overallFiringRateByBlock = nan(numel(blockInds), numel(allSpikeStructs));
for j = 1:numel(blockInds)
    blockStartTime = D.blockStartTimes(blockInds(j));
    blockStopTime = D.blockStopTimes(blockInds(j));
    
    for i = 1:numel(allSpikeStructs)
        overallFiringRateByBlock(j,i) = sum(allSpikeStructs{i}.ts >= blockStartTime & allSpikeStructs{i}.ts <= blockStopTime) / ...
                (blockStopTime - blockStartTime);
    end
end

isFiringNonsparseByBlock = overallFiringRateByBlock >= minFiringRate;
