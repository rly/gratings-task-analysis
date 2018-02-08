function isFiringNonsparseByBlock = findNonsparseBlocks(D, blockInds)
% block x unit
% 0 = sparse firing for that block, 1 = active firing

minFiringRate = 0.1; % Hz

overallFiringRateByBlock = nan(numel(blockInds), numel(D.allSpikeStructs));
for j = 1:numel(blockInds)
    blockStartTime = D.blockStartTimes(blockInds(j));
    blockStopTime = D.blockStopTimes(blockInds(j));
    
    for i = 1:numel(D.allSpikeStructs)
        overallFiringRateByBlock(j,i) = sum(D.allSpikeStructs{i}.ts >= blockStartTime & D.allSpikeStructs{i}.ts <= blockStopTime) / ...
                (blockStopTime - blockStartTime);
    end
end

isFiringNonsparseByBlock = overallFiringRateByBlock >= minFiringRate;
