function [channelDataNorm,eventOnsetClean,isEventOutlier] = preprocessLfps(adjLfps, Fs, channelNames, eventOnset)
% remove events with outlier data and apply CAR and low-pass filtering to adjLfps

doOutlierCheckPlot = 0;

maxAbsYNonOutlier = 1000;%0.25;
outlierMaxSDStep1 = 6; % first-pass threshold before LPF
outlierMaxSDStep2 = 4.5; % second-pass threshold after LPF
outlierCheckWindowOffset = [-0.25 0.3]; % seconds around event

isNanOrig = isnan(adjLfps); % track nans in original data

nEvents = numel(eventOnset);
nChannels = size(adjLfps, 1);
%assert(nChannels == 32);

%% do common average reference first
adjLfpsCAR = adjLfps - nanmean(adjLfps, 1);
clear adjLfps;

%% remove obvious outliers and spike/lever artifacts before low-pass 
% filtering further
fprintf('Removing major outliers and spike/lever artifacts...\n');
for j = 1:nChannels
    fprintf('\tProcessing %s...\n', channelNames{j});
    meanAdjLfp = nanmean(adjLfpsCAR(j,:));
    sdAdjLfp = nanstd(adjLfpsCAR(j,:));
    extremeIndices = adjLfpsCAR(j,:) > meanAdjLfp + outlierMaxSDStep1 * sdAdjLfp | ...
            adjLfpsCAR(j,:) < meanAdjLfp - outlierMaxSDStep1 * sdAdjLfp;
    if sum(extremeIndices) > 0
        fprintf('\tDetected %d outlier time points because voltage in trial is too extreme (> %0.1f SDs). Zeroing them out.\n', ...
                sum(extremeIndices), outlierMaxSDStep1);
    end
    % TODO probably bad just to zero it out instead of smoothly dealing
    % with it or NaNing it. fix this with filtering or case by case
    % replace outlier time points with the average for this channel
    adjLfpsCAR(j,extremeIndices) = meanAdjLfp; 
end

%% low-pass filter data even further
% ideally don't low-pass at 300 Hz and then again at another freq
% use FIR1 filter
fprintf('Low-pass filtering...\n');
hiCutoffFreqCustom = 100; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
adjLfpsCAR(isnan(adjLfpsCAR)) = 0; % zero out nans
adjLfpsLP = filtfilt(bFirLowPassCustom, 1, adjLfpsCAR')';
clear adjLfpsCAR;

% set missing vals back to nan
adjLfpsLP(isNanOrig) = NaN;

%% temp
% for j = 1:nChannels, figure; plot(adjLfpsLP(j,:)); end

%% redo common average reference
adjLfpsLpCAR = adjLfpsLP - nanmean(adjLfpsLP, 1);
clear adjLfpsLP;

%% normalize
% normalize each channel by its mean and standard deviation
% units are now standard deviations away from the mean
channelDataNorm = (adjLfpsLpCAR - nanmean(adjLfpsLpCAR, 2)) ./ nanstd(adjLfpsLpCAR, 0, 2);
assert(all(size(adjLfpsLpCAR) == size(channelDataNorm)));
clear adjLfpsLpCAR;

%% detect events with outlier activity (e.g. amp saturated)
outlierCheckStartIndices = round((eventOnset + outlierCheckWindowOffset(1)) * Fs); 
outlierCheckEndIndices = outlierCheckStartIndices + round(diff(outlierCheckWindowOffset) * Fs) - 1;
outlierCheckT = outlierCheckWindowOffset(1):1/Fs:outlierCheckWindowOffset(2)-1/Fs;
nOutlierCheckTime = numel(outlierCheckT);
assert(all(nOutlierCheckTime == (outlierCheckEndIndices - outlierCheckStartIndices + 1)));

isEventOutlier = false(nEvents, 1);
fprintf('Checking for events with outlier activity (e.g. amp saturated)...\n');
for j = 1:nChannels
    fprintf('\tProcessing %s...\n', channelNames{j});
        
    % plot per channel responses before outlier removal
    if doOutlierCheckPlot
        fHandles(j) = figure_tr_inch(12,6);
        suptitle(sprintf('Outlier Check on Channel %d', j));
        subaxis(1,2,1);
        hold on;
    end
    
    for i = 1:nEvents
        dataInOutlierCheckPeriod = channelDataNorm(j,outlierCheckStartIndices(i):outlierCheckEndIndices(i));
        % mark outlier flashes
        if any(isnan(dataInOutlierCheckPeriod))
            isEventOutlier(i) = 1;
            fprintf('\t\tFlash %d: Detected %d outlier time points because of NaN in trial\n', ...
                    i, sum(isnan(dataInOutlierCheckPeriod)));
        elseif any(abs(dataInOutlierCheckPeriod) > maxAbsYNonOutlier)
            isEventOutlier(i) = 1;
            fprintf('\t\tFlash %d: Detected %d outlier time points because absolute value of voltage in trial is too high (> %0.1f)\n', ...
                    i, sum(abs(dataInOutlierCheckPeriod) > maxAbsYNonOutlier), maxAbsYNonOutlier);
        elseif any(abs(dataInOutlierCheckPeriod) > outlierMaxSDStep2)
            isEventOutlier(i) = 1;
            fprintf('\t\tFlash %d: Detected %d outlier time points because voltage in trial is too extreme (> %0.1f SDs)\n', ...
                    i, sum(abs(dataInOutlierCheckPeriod) > outlierMaxSDStep2), outlierMaxSDStep2);
        end
    
        if doOutlierCheckPlot
            plot(outlierCheckT, dataInOutlierCheckPeriod);
            xlim(outlierCheckWindowOffset);
        end
    end
    
    if doOutlierCheckPlot
        title('Original Responses');
        xlabel('Time from Event Onset');
        ylabel('SDs from Overall Mean');
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        ylim([-8 8]);
    end
end

% remove outlier flashes
fprintf('Detected %d outlier events out of %d.\n', sum(isEventOutlier), nEvents);
eventOnsetClean = eventOnset(~isEventOutlier);

%% plot per channel responses after outlier removal
% TODO recompute baseline (?)
% plot for sanity check
if doOutlierCheckPlot
    nEventsClean = numel(eventOnsetClean);
    startIndices = round((eventOnsetClean + outlierCheckWindowOffset(1)) * Fs); 
    endIndices = startIndices + round(diff(outlierCheckWindowOffset) * Fs) - 1;
    t = outlierCheckWindowOffset(1):1/Fs:outlierCheckWindowOffset(2)-1/Fs;
    nTime = numel(t);
    assert(all(nTime == (endIndices - startIndices + 1)));

    fprintf('Collecting flash responses...\n');
    responses = nan(nChannels, nTime, nEventsClean);
    for j = 1:nChannels
        % normalize each channel by mean baseline across trials and time and
        % std across time of mean baseline across trials
        % TODO: bar plot of average baseline response and std by channel
        fprintf('\tProcessing %s...\n', channelNames{j});

    %     channelDataBaselined{j} = (channelData{j} - nanmean(averageBaselineResponses(j,:))) / nanstd(averageBaselineResponses(j,:));

        % can vectorize???
        for i = 1:nEventsClean
            responses(j,:,i) = channelDataNorm(j,startIndices(i):endIndices(i));
        end

        figure(fHandles(j));
        subaxis(1,2,2);
        hold on;
        plot(t, squeeze(responses(j,:,:)));
        xlim(outlierCheckWindowOffset);
        title('After Outlier Removal');
        xlabel('Time from Flash Onset');
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        ylim([-8 8]);
    end
end