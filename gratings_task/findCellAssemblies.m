function findCellAssemblies(D, processedDataDir, blockName, totalTimeOverall)
%% currently ignores timing during the task

minFiringRateOverall = 0.2;

% find max number of spikes
maxNSpikes = 0;

nUnits = numel(D.allSpikeStructs);
unitLabel = cell(nUnits, 1);
isIncludeUnit = false(nUnits, 1);
for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitLabel{i} = sprintf('%s (%s)', spikeStruct.unitIDChar, spikeStruct.physClass);
    spikeTimes = spikeStruct.ts;
    firingRateOverall = numel(spikeTimes) / totalTimeOverall;
    if firingRateOverall >= minFiringRateOverall
        isIncludeUnit(i) = true;
        if numel(spikeTimes) > maxNSpikes
            maxNSpikes = numel(spikeTimes);
        end
    end
end

% create nan-padded data matrix of spike times
allSpikesMat = nan(nUnits, maxNSpikes);
for i = 1:nUnits
    if isIncludeUnit(i)
        allSpikesMat(i,1:numel(D.allSpikeStructs{i}.ts)) = D.allSpikeStructs{i}.ts;
    end
end
allSpikesMat = trimNanRows(allSpikesMat);
unitLabel = unitLabel(isIncludeUnit);
nUnitsIncluded = size(allSpikesMat, 1);

%% run cell assembly detection script
addpath(genpath('C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\Cell_assembly_detection'));

% assembly search settings
binSizes = [0.005 0.010 0.025];
maxLags = round(0.05 ./ binSizes); % consider at most +/- 50 ms away
filePrefix = sprintf('%s/%s-assembly', processedDataDir, blockName);

tic;
assembly = Main_assemblies_detection(allSpikesMat, maxLags, binSizes, [], [], [], [], [], filePrefix);
toc;

%%
[As_across_bins,As_across_bins_index] = assemblies_across_bins(assembly, binSizes);

saveFileName = sprintf('%s/%s-assembly.mat', processedDataDir, blockName);
savefast(saveFileName, 'allSpikesMat', 'maxLags', 'binSizes', 'assembly', ...
        'As_across_bins', 'As_across_bins_index', 'nUnitsIncluded');

figure_tr_inch(18, 10);
display = 'ordunit';
plotTitle = sprintf('%s - Assembly Assignment', D.sessionName);
[Amatrix,Binvector,Unit_order,As_order] = assembly_assignment_matrix(As_across_bins, ...
        nUnitsIncluded, binSizes, display, plotTitle, unitLabel);
% note colorbar was removed from the created figure temporarily because the
% bars below did not align with them. TODO fix this

plotFileName = sprintf('%s/%s-assembly_assignment.png', processedDataDir, blockName);
export_fig(plotFileName, '-nocrop');
plotFileName = sprintf('%s/%s-assembly_assignment.fig', processedDataDir, blockName);
saveas(gca, plotFileName);
% close;

%%
pruningCriteria = 'biggest';
[As_across_bins_pr,As_across_bins_index_pr] = pruning_across_bins(As_across_bins, As_across_bins_index, nUnitsIncluded, pruningCriteria);

saveFileName = sprintf('%s/%s-assembly-pruned.mat', processedDataDir, blockName);
savefast(saveFileName, ...
        'As_across_bins_pr', 'As_across_bins_index_pr', 'pruningCriteria', 'nUnitsIncluded');

figure_tr_inch(18, 10);
display = 'ordunit';
plotTitle = sprintf('%s - Pruned Assembly Assignment', D.sessionName);
[Amatrix,Binvector,Unit_order,As_order] = assembly_assignment_matrix(As_across_bins_pr, nUnitsIncluded, binSizes, display, plotTitle, unitLabel);

plotFileName = sprintf('%s/%s-assembly_assignment_pruned.png', processedDataDir, blockName);
export_fig(plotFileName, '-nocrop');
plotFileName = sprintf('%s/%s-assembly_assignment_pruned.fig', processedDataDir, blockName);
saveas(gca, plotFileName);
% close;


%% get assembly time courses
lagChoice = 'duration';
activityCount = 'partial';
assemblyActivity = Assembly_activity_function(As_across_bins_pr, assembly, ...
        allSpikesMat, binSizes, lagChoice, activityCount);
nAssembly = size(assemblyActivity, 1);

%% plot assembly time courses
figure_tr_inch(18, 10);
hold on;
for i = 1:nAssembly
    activeIndex = find(assemblyActivity{i}(:,2) == 1);
    plot(assemblyActivity{i}(activeIndex,1), ones(size(activeIndex)) * i, '.', 'MarkerSize', 10)
end

