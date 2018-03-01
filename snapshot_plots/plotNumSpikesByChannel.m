function plotNumSpikesByChannel(channelIDs, allSpikeStructs, processedDataDir, blockName, v)

nChannels = numel(channelIDs);
nSpikesByChannel = zeros(nChannels, 1);
for j = 1:nChannels
    unitsThisCh = findAllUnitsOnCh(allSpikeStructs, channelIDs(j));
    for k = 1:numel(unitsThisCh)
        nSpikesByChannel(j) = nSpikesByChannel(j) + numel(allSpikeStructs{unitsThisCh(k)}.ts);
    end
end

figure_tr_inch(6, 6);
bar(channelIDs, nSpikesByChannel);
title('Number of Spikes by Channel');
xlim(channelIDs([1 end]) + [-0.5 0.5]);
xlabel('Channel Number');
ylabel('Number of Spikes');
plotFileName = sprintf('%s/allSpikes-ch%d-ch%d-%s-numSpikesByChannel_v%d.png', ...
        processedDataDir, channelIDs([1 end]), blockName, v);
fprintf('Saving plot to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
close;