function allUnitsThisChannel = findAllUnitsOnCh(allSpikeStructs, channelID)
allUnitsThisChannel = [];
for j = 1:numel(allSpikeStructs)
    if channelID == allSpikeStructs{j}.channelID
        allUnitsThisChannel = [allUnitsThisChannel j]; %#ok<AGROW>
    end
end