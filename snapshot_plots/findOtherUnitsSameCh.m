function otherUnitsThisChannel = findOtherUnitsSameCh(allSpikeStructs, index)
otherUnitsThisChannel = [];
for j = 1:numel(allSpikeStructs)
    if index ~= j && strcmp(allSpikeStructs{index}.animalName, allSpikeStructs{j}.animalName) && ...
            strcmp(allSpikeStructs{index}.sessionName, allSpikeStructs{j}.sessionName) && ...
            strcmp(allSpikeStructs{index}.areaName, allSpikeStructs{j}.areaName) && ...
            allSpikeStructs{index}.channelID == allSpikeStructs{j}.channelID
        otherUnitsThisChannel = [otherUnitsThisChannel j]; %#ok<AGROW>
    end
end