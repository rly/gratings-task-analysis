function allUnitsThisChannel = findAllUnitsSameCh(allSpikeStructs, index)
allUnitsThisChannel = [];
for j = 1:numel(allSpikeStructs)
    if strcmp(allSpikeStructs{index}.sessionName, allSpikeStructs{j}.sessionName) && ...
            strcmp(allSpikeStructs{index}.areaName, allSpikeStructs{j}.areaName) && ...
            allSpikeStructs{index}.channelID == allSpikeStructs{j}.channelID && ...
            ~allSpikeStructs{j}.isMUA
        allUnitsThisChannel = [allUnitsThisChannel j]; %#ok<AGROW>
    end
end