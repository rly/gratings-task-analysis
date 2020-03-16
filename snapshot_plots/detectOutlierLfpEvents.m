function [eventOnsetClean,isEventOutlier] = detectOutlierLfpEvents(channelData, Fs, ...
        channelNames, eventOnset, outlierCheckWindowOffset, ...
        processedDataDir, plotFileNamePrefix, rfMappingMode, v)

doOutlierCheckPlot = 1;
outlierMaxSDStep2 = 8; % second-pass threshold after LPF

nEvents = numel(eventOnset);
nChannels = size(channelData, 1);

%% detect events with outlier activity (e.g. NaN, amp saturated)
outlierCheckStartIndices = round((eventOnset + outlierCheckWindowOffset(1)) * Fs); 
outlierCheckEndIndices = outlierCheckStartIndices + round(diff(outlierCheckWindowOffset) * Fs) - 1;
outlierCheckT = outlierCheckWindowOffset(1):1/Fs:outlierCheckWindowOffset(2)-1/Fs;
nOutlierCheckTime = numel(outlierCheckT);
assert(all(nOutlierCheckTime == (outlierCheckEndIndices - outlierCheckStartIndices + 1)));

isEventOutlier = false(nEvents, 1);

% check for nans around each event in any channel 
fprintf('Checking for events with NaN...\n');
isEventNan = false(nEvents, 1);
for n = 1:nEvents
    dataInOutlierCheckPeriod = channelData(:,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
    if any(any(isnan(dataInOutlierCheckPeriod)))
        isEventNan(n) = 1;
        fprintf('\t\tFlash %d: Detected %d NaN time points because in trial\n', ...
                n, sum(any(isnan(dataInOutlierCheckPeriod))));
    end
end
isEventOutlier = isEventOutlier | isEventNan;

if doOutlierCheckPlot && any(isEventNan)
    figure_tr_inch(6, 6);
    suptitle('NaN Check - Each Line is Separate Channel');
    hold on;
    dataInOutlierCheckPeriod = channelData(:,outlierCheckStartIndices(isEventNan):outlierCheckEndIndices(isEventNan));
    plot(outlierCheckT, dataInOutlierCheckPeriod);

    title('Original Responses');
    xlabel('Time from Event Onset');
    ylabel('SDs from Overall Mean');
    plot([0 0], [-1000 1000], 'k-');
    xlim(outlierCheckWindowOffset);
    ylim([-8 8]);

    plotFileName = sprintf('%s/%s-rfm_mode%d-nanCheck-v%d.png', ...
            processedDataDir, plotFileNamePrefix, rfMappingMode, v);
    fprintf('\t\tSaving NaN check plot to %s...\n', plotFileName);
    %export_fig(plotFileName, '-nocrop');
    close;
end


fprintf('Checking for events with outlier activity (e.g. amp saturated)...\n');
for j = 1:nChannels
    fprintf('\tProcessing %s...\n', channelNames{j});
    
    isEventOutlierThisChannel = false(nEvents, 1);
    for n = 1:nEvents
        dataInOutlierCheckPeriod = channelData(j,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
        % mark outlier flashes
        if any(abs(dataInOutlierCheckPeriod) > outlierMaxSDStep2)
            isEventOutlierThisChannel(n) = 1;
            fprintf('\t\tFlash %d: Detected %d outlier time points because voltage in trial is too extreme (> %0.1f SDs)\n', ...
                    n, sum(abs(dataInOutlierCheckPeriod) > outlierMaxSDStep2), outlierMaxSDStep2);
        end
    end
    
    isEventOutlier = isEventOutlier | isEventOutlierThisChannel;
    
    if doOutlierCheckPlot && any(isEventOutlierThisChannel)
        figure_tr_inch(6, 6);
        suptitle(sprintf('Outlier Check on Channel %s', channelNames{j}));
        hold on;
        % plot all responses with thin line
        for n = 1:nEvents
            dataInOutlierCheckPeriod = channelData(j,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
            plot(outlierCheckT, dataInOutlierCheckPeriod);
        end
        for n = 1:nEvents
            if isEventOutlierThisChannel(n)
                dataInOutlierCheckPeriod = channelData(j,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
                % plot outlier responses over other lines, with thicker line
                plot(outlierCheckT, dataInOutlierCheckPeriod, 'LineWidth', 3);
            end
        end
        
        title('Original Responses');
        xlabel('Time from Event Onset');
        ylabel('SDs from Overall Mean');
        plot([0 0], [-1000 1000], 'k-');
        xlim(outlierCheckWindowOffset);
        ylim([-8 8]);
        
        plotFileName = sprintf('%s/%s-rfm_mode%d-%s-outlierCheck-v%d.png', ...
                processedDataDir, plotFileNamePrefix, rfMappingMode, channelNames{j}, v);
        fprintf('\t\tSaving outlier check plot to %s...\n', plotFileName);
        %export_fig(plotFileName, '-nocrop');
        close;
    end
end

% remove outlier flashes
fprintf('Detected %d outlier events out of %d (%d%%).\n', sum(isEventOutlier), ...
        nEvents, round(sum(isEventOutlier)/nEvents * 100));
eventOnsetClean = eventOnset(~isEventOutlier);
if sum(isEventOutlier) > 0.1 * nEvents
    error('More than 10%% of flashes have outlier responses.');
end
