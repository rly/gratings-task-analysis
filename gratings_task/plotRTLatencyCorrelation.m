function plotRTLatencyCorrelation(eventStruct, isLocUsed, rtByLoc, spikeTimesByLoc, inRFLocs, plotFileName)

% currently referencing a unit's ES though it shouldn't matter
kernelSigma = 0.01;
eventPsthT = eventStruct.t;
baselineWindowOffset = [-0.3 0];
baselineWindow = eventStruct.window(1) + baselineWindowOffset; % from array onset
baselineIndices = getTimeLogicalWithTolerance(eventPsthT, baselineWindow);
latencyWindowOffset = [0.025 0.125];
latencyWindow = eventStruct.window(1) + latencyWindowOffset;
latencyIndices = getTimeLogicalWithTolerance(eventPsthT, latencyWindow);
minPeakForLatency = 3; % SDs from mean
maxTroughForLatency = -minPeakForLatency;
nLoc = numel(isLocUsed);

xBounds = [0 0.13];
figure_tr_inch(15, 5);

for k = 1:nLoc
    if isLocUsed(k)
        nTrials = size(spikeTimesByLoc{k}, 2);
        assert(numel(rtByLoc{k}) == nTrials);
        latencies = nan(nTrials, 1);
%         f1 = figure;
%         hold on;
        unitCondition = inRFLocs == k;
        for i = 1:nTrials
            psthResponse = fixedPsth(spikeTimesByLoc{k}(unitCondition,i)', kernelSigma, 0, eventPsthT);
            baselineResponseByPsth = mean(psthResponse(baselineIndices));
            baselineResponseSDOverTimeByPsth = std(psthResponse(baselineIndices));
            normPsthResponse = (psthResponse - baselineResponseByPsth) / baselineResponseSDOverTimeByPsth;
            % check peak
            peakLatencyInfo = computeLatencyPeakMethod(normPsthResponse, eventPsthT, ...
                    eventStruct.window, latencyIndices, 0, minPeakForLatency, 0); 
            % check trough
            troughLatencyInfo = computeLatencyPeakMethod(normPsthResponse, eventPsthT, ...
                    eventStruct.window, latencyIndices, 0, maxTroughForLatency, 1);

            % get the earlier one
            latencyInfo = getEarlierPeakTroughLatencyInfo(peakLatencyInfo, troughLatencyInfo);
            latencies(i) = latencyInfo.latency;
            
%             figure;
%             hold on;
%             plot(eventPsthT - eventStruct.window(1), psthResponse);
%             if ~isnan(latencies(i))
%                 plot(latencies(i), psthResponse(latencyInfo.latencyTInd), 'o', ...
%                         'Color', 'k', 'MarkerSize', 8, 'LineWidth', 2);
%                 plot(latencyInfo.timeEventToPeak, psthResponse(latencyInfo.peakTInd), 'o', ...
%                         'Color', 'k', 'MarkerSize', 8, 'LineWidth', 2);
%                 figure(f1);
%                 plot(latencyInfo.timeEventToPeak, rtByLoc{k}(i), '.', 'MarkerSize', 20);
%             end
%             pause;
        end

        subaxis(1, nLoc, k);
        hold on;
        plot(latencies, rtByLoc{k}, '.', 'MarkerSize', 20);
        LM = fitlm(latencies, rtByLoc{k});
        predY = predict(LM, xBounds');
        plot(xBounds, predY, 'LineWidth', 2);
        xlim(xBounds);
        xlabel('Response Latency (s)');
        ylabel('Response Time (s)');
        [r,p] = corr(latencies(~isnan(latencies)), rtByLoc{k}(~isnan(latencies)), 'type', 'Spearman');
        title(sprintf('Loc %d, NUnits = %d, NTrials = %d, r = %0.2f, p=%0.3f', k, sum(unitCondition), sum(~isnan(latencies)), r, p));
    end
end

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end