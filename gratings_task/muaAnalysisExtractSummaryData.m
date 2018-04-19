function muaAnalysisExtractSummaryData(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)

v = 11;
tic;

% for preallocation. make sure this is an underestimate or equal to actual
% number of units saved
nUnitsApprox = 1; 

unitNames = cell(nUnitsApprox, 1);
isSignificantResponseVsBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
isSignificantResponseVsPreviousPeriod = false(nUnitsApprox, 4);
isSignificantResponseVsBootstrapBaseline = false(nUnitsApprox, 6); % 6 periods > baseline
isSignificantResponseVsBootstrapPreviousPeriod = false(nUnitsApprox, 4);
isSignificantSelectivity = false(nUnitsApprox, 5); % 5 periods info rate
cueResponseVsBaselineDirection = zeros(nUnitsApprox, 1);
infoRates = nan(nUnitsApprox, 5); % 5 periods
diffRates = nan(nUnitsApprox, 3); % 2 delay periods + array response
attnIndices = nan(nUnitsApprox, 3); % 2 delay periods + array response
localization = cell(nUnitsApprox, 1);
isInVPulvinar = false(nUnitsApprox, 1);
isInDPulvinar = false(nUnitsApprox, 1);
preCueBaselineExpFit = cell(nUnitsApprox, 1);
earlyPreExitFixationSlope = nan(nUnitsApprox, 1);
latePreExitFixationSlope = nan(nUnitsApprox, 1);
earlyPreExitFixationWindowOffset = [-0.2 -0.05];
latePreExitFixationWindowOffset = [-0.05 0];

diffTargetDimLatencySplitThirdsRT = nan(nUnitsApprox, 1);
rtFiringRateStruct = struct();

nLoc = 4;
arrayOnsetRelBalSpikeTimesByLoc = cell(nLoc, 1);
arrayOnsetHoldBalSpikeTimesByLoc = cell(nLoc, 1);
targetDimBalSpikeTimesByLoc = cell(nLoc, 1);
inRFLocs = nan(nUnitsApprox, 1);
exRFLocs = nan(nUnitsApprox, 1);

arrayHoldBalLatencyInRF = nan(nUnitsApprox, 1);
arrayHoldBalLatencyExRF = nan(nUnitsApprox, 1);
targetDimBalLatencyInRF = nan(nUnitsApprox, 1);
targetDimBalLatencyExRF = nan(nUnitsApprox, 1);

averageFiringRatesBySpdf = struct();
averageFiringRatesByCount = struct();

targetDimDelayHoldInRFRateAll = cell(nUnitsApprox, 1);
rtHoldInRFAll = cell(nUnitsApprox, 1);

unitCount = 0;
minFiringRate = 2; % use only cells with a time-locked response > 1 Hz in any window
statAlpha = 0.05; % account for multiple comparisons later
% should also be running a lot of shuffle tests given the number of trials

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');

outputDir = sprintf('%s/%s/', processedDataRootDir, 'MUA_GRATINGS_SUMMARY');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%%
taskName = 'GRATINGS';
scriptName = 'MUA_GRATINGS';
isZeroDistractors = 0;
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, muaChannelsToLoad, taskName, scriptName, 1, 0);
sessionName = R.sessionName;
areaName = R.areaName;

%% summarize across population of cells    
totalTimeOverall = sum(D.blockStopTimes(R.blockIndices) - D.blockStartTimes(R.blockIndices));
minFiringRateOverall = 0.2;

nUnits = numel(muaChannelsToLoad);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing Session %s %s, Channels %d-%d\n', sessionName, areaName, muaChannelsToLoad([1 end]));
fprintf('Processing %d units...\n', nUnits);

for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    firingRateOverall = numel(spikeTimes) / totalTimeOverall;

    if firingRateOverall >= minFiringRateOverall
        saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                processedDataDir, unitName, blockName, v);
        fprintf('Processing %s...\n', saveFileName);

        ES = load(saveFileName);

        unitCount = unitCount + 1;

        if (any(ES.averageFiringRatesBySpdf.preEnterFixation.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.postEnterFixation.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.postEnterFixationLate.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.preCueBaseline.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.cueResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.arrayResponseHold.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.arrayResponseRel.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.targetDimDelay.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.targetDimResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.preExitFixation.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.postExitFixation.byLoc >= minFiringRate))

            unitNames{unitCount} = unitName;

%             if ismember(spikeStruct.channelID, R.pldChannels)
%                 localization{unitCount} = 'PLd';
%             elseif ismember(spikeStruct.channelID, R.plvChannels)
%                 localization{unitCount} = 'PLv';
%             elseif ismember(spikeStruct.channelID, R.pmChannels)
%                 localization{unitCount} = 'PM';
%             elseif ismember(spikeStruct.channelID, R.piChannels)
%                 localization{unitCount} = 'PI';
%             else
                localization{unitCount} = '';
%             end

            isInVPulvinar(unitCount) = 0;
            isInDPulvinar(unitCount) = 0;
            if ismember(spikeStruct.channelID, R.vPulChannels)
                isInVPulvinar(unitCount) = 1;
                localization{unitCount} = 'vPul'; % TEMP
            end
            if ismember(spikeStruct.channelID, R.dPulChannels)
                isInDPulvinar(unitCount) = 1;
                localization{unitCount} = 'dPul'; % TEMP
            end
            if ismember(spikeStruct.channelID, R.vPulChannels) && ...
                    ismember(spikeStruct.channelID, R.dPulChannels)
                error('Channel %d cannot be in both vPul and dPul', spikeStruct.channelID);
            end

            nLocUsed = sum(ES.isLocUsed);

            % significance test using rank sum test. correct for multiple
            % comparisons by comparing p < statAlpha / nLocUsed
            isSignificantResponseVsBaseline(unitCount,:) = [...
                    getMinPValueOverLoc([ES.cueResponseVsBaselineStatsByLoc.rankSum]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.cueTargetDelayVsBaselineStatsByLoc.rankSum]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.arrayResponseHoldVsBaselineStatsByLoc.rankSum]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.targetDimDelayVsBaselineStatsByLoc.rankSum]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.targetDimResponseVsBaselineStatsByLoc.rankSum]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.preExitFixationVsBaselineStatsByLoc.rankSum]) < statAlpha / nLocUsed];
                
            % significance test using rank sum test. correct for multiple
            % comparisons by comparing p < statAlpha / nLocUsed
            isSignificantResponseVsPreviousPeriod(unitCount,:) = [...
                    getMinPValueOverLoc([ES.cueResponseVsBaselineStatsByLoc.signRank]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.arrayResponseHoldVsPrevStatsByLoc.signRank]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.targetDimResponseVsPrevStatsByLoc.signRank]) < statAlpha / nLocUsed ...
                    getMinPValueOverLoc([ES.preExitFixationVsPrevStatsByLoc.signRank]) < statAlpha / nLocUsed];
            
            % significance test using non-parametric test vs bootstrapped
            % baseline. correct for multiple comparisons by comparing p <
            % statAlpha / nLocUsed
%             isSignificantResponseVsBootstrapBaseline(unitCount,:) = [...
%                     getMinPValueOverLoc([ES.cueResponseVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.cueTargetDelayVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.arrayResponseHoldVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.targetDimDelayVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.targetDimResponseVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.preExitFixationVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed];
%                 
%             isSignificantResponseVsBootstrapPreviousPeriod(unitCount,:) = [...
%                     getMinPValueOverLoc([ES.cueResponseVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.arrayResponseHoldVsPrevStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.targetDimResponseVsPrevStatsByLoc.bootstrap]) < statAlpha / nLocUsed ...
%                     getMinPValueOverLoc([ES.preExitFixationVsPrevStatsByLoc.bootstrap]) < statAlpha / nLocUsed];
%             
            % for significant responses, look at largest ABSOLUTE 
            % difference from baseline to determine direction of response
            if isSignificantResponseVsPreviousPeriod(unitCount,1)
                respDiffFromBaseline = ES.averageFiringRatesBySpdf.cueResponse.byLoc - ES.averageFiringRatesBySpdf.preCueBaseline.byLoc;
%                 isSig = arrayfun(@(x) x.p, [ES.cueResponseVsBaselineStatsByLoc.bootstrap]) < statAlpha / nLocUsed;
%                 assert(~isempty(union(ES.inRFLocByExtreme, find(isSig)))); % check InRF Extreme loc has significant response
                cueResponseVsBaselineDirection(unitCount) = sign(respDiffFromBaseline(ES.inRFLocByExtreme));
                clear respDiffFromBaseline;% isSig;
            else
                cueResponseVsBaselineDirection(unitCount) = 0;
            end
            
            % determine slope of firing rate within baseline period
            % not used currently
%             preCueBaselineExpFit{unitCount} = computeExpFiringRateBySpdf([-0.25 0], ES.cueOnset);
            
            % determine slope of firing rate prior to exit fixation saccade
            earlyPreExitFixationSlopeStruct = computeSlopeFiringRateBySpdf(earlyPreExitFixationWindowOffset, ES.exitFixation);
            latePreExitFixationSlopeStruct = computeSlopeFiringRateBySpdf(latePreExitFixationWindowOffset, ES.exitFixation);
            earlyPreExitFixationSlope(unitCount) = earlyPreExitFixationSlopeStruct.all;
            latePreExitFixationSlope(unitCount) = latePreExitFixationSlopeStruct.all;
            
            isSignificantSelectivity(unitCount,:) = [...
                    ES.cueResponseInfoRate.permutation.p < statAlpha ...
                    ES.cueTargetDelayInfoRate.permutation.p < statAlpha ...
                    ES.arrayResponseHoldInfoRate.permutation.p < statAlpha ...
                    ES.targetDimDelayInfoRate.permutation.p < statAlpha ...
                    ES.targetDimResponseInfoRate.permutation.p < statAlpha];

            infoRates(unitCount,:) = [...
                    ES.cueResponseInfoRate.infoRate ...
                    ES.cueTargetDelayInfoRate.infoRate ...
                    ES.arrayResponseHoldInfoRate.infoRate ...
                    ES.targetDimDelayInfoRate.infoRate ...
                    ES.targetDimResponseInfoRate.infoRate];

            diffRates(unitCount,:) = [...
                    ES.cueTargetDelayAttnStats.diff.actual ...
                    ES.arrayResponseHoldAttnStats.diff.actual ...
                    ES.targetDimDelayAttnStats.diff.actual];    
            
            attnIndices(unitCount,:) = [...
                    ES.cueTargetDelayAttnStats.ai.actual ...
                    ES.arrayResponseHoldAttnStats.ai.actual ...
                    ES.targetDimDelayAttnStats.ai.actual];  

            inRFLoc = ES.inRFLoc;
            exRFLoc = ES.exRFLoc;
            
            %%
            % extract firing rates (count-based) at InRF and ExRF locations
            % note: cueTargetDelayLongWindowOffset = [-0.4 0];
            % note: targetDimDelayLongWindowOffset = [-0.4 0];
            
            rtRelInRF = ES.UE.rt(ES.UE.cueLoc == inRFLoc & ES.UE.isRelBal);
            rtRelExRF = ES.UE.rt(ES.UE.cueLoc == exRFLoc & ES.UE.isRelBal);
            rtHoldInRF = ES.UE.rt(ES.UE.cueLoc == inRFLoc & ES.UE.isHoldBal);
            rtHoldExRF = ES.UE.rt(ES.UE.cueLoc == exRFLoc & ES.UE.isHoldBal);
            targetDimInRF = ES.UE.targetDimMatch(ES.UE.cueLoc == inRFLoc & ES.UE.isHoldBal);
            targetDimExRF = ES.UE.targetDimMatch(ES.UE.cueLoc == exRFLoc & ES.UE.isHoldBal);
            
            % use long balanced trials only
            cueTargetDelayRelInRFRate = ES.averageFiringRatesByCount.cueTargetDelayRelBalLong.trialRateByLoc{inRFLoc};
            cueTargetDelayRelExRFRate = ES.averageFiringRatesByCount.cueTargetDelayRelBalLong.trialRateByLoc{exRFLoc};
            
            cueTargetDelayHoldInRFRate = ES.averageFiringRatesByCount.cueTargetDelayHoldBalLong.trialRateByLoc{inRFLoc};
            cueTargetDelayHoldExRFRate = ES.averageFiringRatesByCount.cueTargetDelayHoldBalLong.trialRateByLoc{exRFLoc};

            targetDimDelayHoldInRFRate = ES.averageFiringRatesByCount.targetDimDelayBalLong.trialRateByLoc{inRFLoc};
            targetDimDelayHoldExRFRate = ES.averageFiringRatesByCount.targetDimDelayBalLong.trialRateByLoc{exRFLoc};

            nTrialsRelInRF = numel(rtRelInRF);
            topThirdIndicesRelInRF = round(nTrialsRelInRF * 2/3)+1:nTrialsRelInRF;
            bottomThirdIndicesRelInRF = 1:round(nTrialsRelInRF * 1/3);
            
            nTrialsRelExRF = numel(rtRelExRF);
            topThirdIndicesRelExRF = round(nTrialsRelExRF * 2/3)+1:nTrialsRelExRF;
            bottomThirdIndicesRelExRF = 1:round(nTrialsRelExRF * 1/3);
            
            nTrialsHoldInRF = numel(rtHoldInRF);
            topThirdIndicesHoldInRF = round(nTrialsHoldInRF * 2/3)+1:nTrialsHoldInRF;
            bottomThirdIndicesHoldInRF = 1:round(nTrialsHoldInRF * 1/3);

            nTrialsHoldExRF = numel(rtHoldExRF);
            topThirdIndicesHoldExRF = round(nTrialsHoldExRF * 2/3)+1:nTrialsHoldExRF;
            bottomThirdIndicesHoldExRF = 1:round(nTrialsHoldExRF * 1/3);
            
            assert(nTrialsRelInRF == numel(cueTargetDelayRelInRFRate));
            assert(nTrialsRelExRF == numel(cueTargetDelayRelExRFRate));
            assert(all(nTrialsHoldInRF == [numel(cueTargetDelayHoldInRFRate) numel(targetDimDelayHoldInRFRate)]));
            assert(all(nTrialsHoldExRF == [numel(cueTargetDelayHoldExRFRate) numel(targetDimDelayHoldExRFRate)]));
            
            [~,sortRTHoldInRFInd] = sortBreakOrder(rtHoldInRF);
            [~,sortRTHoldExRFInd] = sortBreakOrder(rtHoldExRF);
            
            % slowest RTs
            targetDimByLocSlowThirdRT = cell(nLoc, 1);
            targetDimByLocSlowThirdRT{inRFLoc} = targetDimInRF(sortRTHoldInRFInd(topThirdIndicesHoldInRF));
            targetDimByLocSlowThirdRT{exRFLoc} = targetDimExRF(sortRTHoldExRFInd(topThirdIndicesHoldExRF));
            targetDimSlowThirdRT.window = ES.targetDimBal.window;
            targetDimSlowThirdRT.spdfWindowOffset = ES.targetDimBal.spdfWindowOffset;
            targetDimSlowThirdRT = createTimeLockedSpdf(ES.spikeTs, [], targetDimByLocSlowThirdRT, targetDimSlowThirdRT, ES.kernelSigma);
            targetDimSlowThirdRT = computeResponseLatencyByLoc(targetDimSlowThirdRT, ES.isLocUsed);
            
            % fastest RTs
            targetDimByLocFastThirdRT = cell(nLoc, 1);
            targetDimByLocFastThirdRT{inRFLoc} = targetDimInRF(sortRTHoldInRFInd(bottomThirdIndicesHoldInRF));
            targetDimByLocFastThirdRT{exRFLoc} = targetDimExRF(sortRTHoldExRFInd(bottomThirdIndicesHoldExRF));
            targetDimFastThirdRT.window = ES.targetDimBal.window;
            targetDimFastThirdRT.spdfWindowOffset = ES.targetDimBal.spdfWindowOffset;
            targetDimFastThirdRT = createTimeLockedSpdf(ES.spikeTs, [], targetDimByLocFastThirdRT, targetDimFastThirdRT, ES.kernelSigma);
            targetDimFastThirdRT = computeResponseLatencyByLoc(targetDimFastThirdRT, ES.isLocUsed);

            % InRF only
            diffTargetDimLatencySplitThirdsRT(unitCount) = targetDimSlowThirdRT.latencyInfoByLoc{inRFLoc}.latency - ...
                    targetDimFastThirdRT.latencyInfoByLoc{inRFLoc}.latency;
            if ~isnan(diffTargetDimLatencySplitThirdsRT(unitCount))
                figure_tr_inch(6, 5);
                hold on;
                targetDimT = targetDimSlowThirdRT.t - targetDimSlowThirdRT.window(1);
                plot(targetDimT, targetDimSlowThirdRT.spdfByLoc(inRFLoc,:), 'LineWidth', 2);
                plot(targetDimT, targetDimFastThirdRT.spdfByLoc(inRFLoc,:), 'LineWidth', 2);
                xlabel('Time from Target Dimming (s)');
                ylabel('Estimated Firing Rate (Hz)');
                title(sprintf('%s Target Dim InRF - Slow vs Fast RTs', unitName), 'Interpreter', 'none');
                legend({'Slow Third of RTs', 'Fast Third of RTs'});
                plotFileName = sprintf('%s/%s-%s-targetDimSlowVsFastRT-v%d.png', processedDataDir, unitName, blockName, v);
                fprintf('\tSaving figure to file %s...\n', plotFileName);
                export_fig(plotFileName, '-nocrop');
            end
            
            % algorithm:
            % sort trials according to firing rate during each delay period
            % sort the RTs based on the above sorting
            % (alternatively sort based on RT and apply to FR)

            % WARNING: many of these firing rates are 0 because of the short
            % time window and the sparse firing. even with using
            % sortBreakOrder, this could lead to problems
            % this plot is commented out because of the above WARNING
            % repeated runs of this code will yield vastly different p
            % values. TODO need to account for this.
            [~,sortCTDelayRelInRFInd] = sortBreakOrder(cueTargetDelayRelInRFRate);
            rtRelInRFSortedByCTDelay = rtRelInRF(sortCTDelayRelInRFInd);
            
            [~,sortCTDelayRelExRFInd] = sortBreakOrder(cueTargetDelayRelExRFRate);
            rtRelExRFSortedByCTDelay = rtRelExRF(sortCTDelayRelExRFInd);
            
            [~,sortCTDelayHoldInRFInd] = sortBreakOrder(cueTargetDelayHoldInRFRate);
            rtHoldInRFSortedByCTDelay = rtHoldInRF(sortCTDelayHoldInRFInd);
            
            [~,sortCTDelayHoldExRFInd] = sortBreakOrder(cueTargetDelayHoldExRFRate);
            rtHoldExRFSortedByCTDelay = rtHoldExRF(sortCTDelayHoldExRFInd);
            
            [~,sortTDDelayHoldInRFInd] = sortBreakOrder(targetDimDelayHoldInRFRate);
            rtHoldInRFSortedByTDDelay = rtHoldInRF(sortTDDelayHoldInRFInd);
            
            [~,sortTDDelayHoldExRFInd] = sortBreakOrder(targetDimDelayHoldExRFRate);
            rtHoldExRFSortedByTDDelay = rtHoldExRF(sortTDDelayHoldExRFInd);
            
            % make session-wise RT plots while processing the first unit
            if unitCount == 1
                assert(all(ES.UE.rt >= 0.3 & ES.UE.rt <= 0.8)); 
                checkRTStatAlpha = 0.05;
                plotFileName = sprintf('%s/%s-sessionInd%d-rtDist-v%d.png', outputDir, sessionName, sessionInd, v);
                plotRTDistribution(rtRelInRF, rtRelExRF, rtHoldInRF, rtHoldExRF, ...
                        checkRTStatAlpha, sessionName, isZeroDistractors, plotFileName);
                    
                % there should not be any difference between RTs on InRF and
                % ExRF trials. otherwise there is significant spatial bias
                % in this session.
                p = ranksum(rtRelInRF, rtRelExRF);
                if p < checkRTStatAlpha
                    warning('Release Trial median RT is significantly different between InRF and ExRF conditions (p = %0.3f)\n', p);
                end
                p = ranksum(rtHoldInRF, rtHoldExRF);
                if p < checkRTStatAlpha
                    warning('Hold Trial median RT is significantly different between InRF and ExRF conditions (p = %0.3f)\n', p);
                end
                
                assert(all(ES.UE.cueTargetDelayDur >= 450 & ES.UE.cueTargetDelayDur <= 850));
                assert(all(ES.UE.cueTargetDelayDur >= 250 & ES.UE.cueTargetDelayDur <= 1150));
                figure_tr_inch(12, 5);
                subaxis(1, 2, 1);
                histogram(ES.UE.cueTargetDelayDur, 450:25:850);
                xlabel('Cue-Target Delay Duration (ms)');
                ylabel('Number of Trials');
                
                subaxis(1, 2, 2);
                histogram(ES.UE.targetDimDelayDur, 250:25:1150);
                xlabel('Target-Dim Delay Duration (ms)');
                ylabel('Number of Trials');
                
                suptitle(sprintf('Delay Period Distributions: %s', sessionName));
                
                plotFileName = sprintf('%s/%s-sessionInd%d-delayDurDist-v%d.png', outputDir, sessionName, sessionInd, v);
                fprintf('\tSaving figure to file %s...\n', plotFileName);
                export_fig(plotFileName, '-nocrop');
            end
%             
% 
%             % rel trials, RT on trials with lowest 1/3 of firing rates in
%             % CT delay vs RT on trials with highest 1/3 of firing rates
%             % InRF and ExRF trials separately
%             
%             % rel trials CT delay
%             [~,rtRelInRFSortedByCTDelayThirdsPVal] = ttest2(rtRelInRFSortedByCTDelay(topThirdIndicesRelInRF), ...
%                     rtRelInRFSortedByCTDelay(bottomThirdIndicesRelInRF));
%             [~,rtRelExRFSortedByCTDelayThirdsPVal] = ttest2(rtRelExRFSortedByCTDelay(topThirdIndicesRelExRF), ...
%                     rtRelExRFSortedByCTDelay(bottomThirdIndicesRelExRF));
%             
%             % hold trials CT delay
%             [~,rtHoldInRFSortedByCTDelayThirdsPVal] = ttest2(rtHoldInRFSortedByCTDelay(topThirdIndicesHoldInRF), ...
%                     rtHoldInRFSortedByCTDelay(bottomThirdIndicesHoldInRF));
%             [~,rtHoldExRFSortedByCTDelayThirdsPVal] = ttest2(rtHoldExRFSortedByCTDelay(topThirdIndicesHoldExRF), ...
%                     rtHoldExRFSortedByCTDelay(bottomThirdIndicesHoldExRF));
% 
%             % hold trials TD delay
%             [~,rtHoldInRFSortedByTDDelayThirdsPVal] = ttest2(rtHoldInRFSortedByTDDelay(topThirdIndicesHoldInRF), ...
%                     rtHoldInRFSortedByTDDelay(bottomThirdIndicesHoldInRF));
%             [~,rtHoldExRFSortedByTDDelayThirdsPVal] = ttest2(rtHoldExRFSortedByTDDelay(topThirdIndicesHoldExRF), ...
%                     rtHoldExRFSortedByTDDelay(bottomThirdIndicesHoldExRF));
% 
%             figure_tr_inch(15, 10);
%             textParamsNormal = {'FontSize', 8, 'Units', 'normalized'};
%             textParamsBold = {'FontSize', 8, 'Units', 'normalized', 'FontWeight', 'bold'};
%             corrStatAlpha = 0.05;
%             
%             subaxis(2, 3, 1, 'SV', 0.1);
%             hold on;
%             histogram(rtRelInRFSortedByCTDelay(topThirdIndicesRelInRF), binEdges);
%             histogram(rtRelInRFSortedByCTDelay(bottomThirdIndicesRelInRF), binEdges);
%             title({'RT distribution of Top and Bottom Third', 'Extreme Cue-Target Delay Firing Rates Rel InRF'});
%             if rtRelInRFSortedByCTDelayThirdsPVal < corrStatAlpha
%                 textParams = textParamsBold;
%             else
%                 textParams = textParamsNormal;
%             end
%             text(0.85, 0.95, sprintf('p = %0.3f', rtRelInRFSortedByCTDelayThirdsPVal), textParams{:});
%             
%             subaxis(2, 3, 4);
%             hold on;
%             histogram(rtRelExRFSortedByCTDelay(topThirdIndicesRelExRF), binEdges);
%             histogram(rtRelExRFSortedByCTDelay(bottomThirdIndicesRelExRF), binEdges);
%             title({'RT distribution of Top and Bottom Third', 'Extreme Cue-Target Delay Firing Rates Rel ExRF'});
%             if rtRelExRFSortedByCTDelayThirdsPVal < corrStatAlpha
%                 textParams = textParamsBold;
%             else
%                 textParams = textParamsNormal;
%             end
%             text(0.85, 0.95, sprintf('p = %0.3f', rtRelExRFSortedByCTDelayThirdsPVal), textParams{:});
%             
%             subaxis(2, 3, 2);
%             hold on;
%             histogram(rtHoldInRFSortedByCTDelay(topThirdIndicesHoldInRF), binEdges);
%             histogram(rtHoldInRFSortedByCTDelay(bottomThirdIndicesHoldInRF), binEdges);
%             title({'RT distribution of Top and Bottom Third', 'Extreme Cue-Target Delay Firing Rates Hold InRF'});
%             if rtHoldInRFSortedByCTDelayThirdsPVal < corrStatAlpha
%                 textParams = textParamsBold;
%             else
%                 textParams = textParamsNormal;
%             end
%             text(0.85, 0.95, sprintf('p = %0.3f', rtHoldInRFSortedByCTDelayThirdsPVal), textParams{:});
%             
%             subaxis(2, 3, 5);
%             hold on;
%             histogram(rtHoldExRFSortedByCTDelay(topThirdIndicesHoldExRF), binEdges);
%             histogram(rtHoldExRFSortedByCTDelay(bottomThirdIndicesHoldExRF), binEdges);
%             title({'RT distribution of Top and Bottom Third', 'Extreme Cue-Target Delay Firing Rates Hold ExRF'});
%             if rtHoldExRFSortedByCTDelayThirdsPVal < corrStatAlpha
%                 textParams = textParamsBold;
%             else
%                 textParams = textParamsNormal;
%             end
%             text(0.85, 0.95, sprintf('p = %0.3f', rtHoldExRFSortedByCTDelayThirdsPVal), textParams{:});
%             
%             subaxis(2, 3, 3);
%             hold on;
%             histogram(rtHoldInRFSortedByTDDelay(topThirdIndicesHoldInRF), binEdges);
%             histogram(rtHoldInRFSortedByTDDelay(bottomThirdIndicesHoldInRF), binEdges);
%             title({'RT distribution of Top and Bottom Third', 'Extreme Target-Dim Delay Firing Rates Hold InRF'});
%             if rtHoldInRFSortedByTDDelayThirdsPVal < corrStatAlpha
%                 textParams = textParamsBold;
%             else
%                 textParams = textParamsNormal;
%             end
%             text(0.85, 0.95, sprintf('p = %0.3f', rtHoldInRFSortedByTDDelayThirdsPVal), textParams{:});
%             
%             subaxis(2, 3, 6);
%             hold on;
%             histogram(rtHoldExRFSortedByTDDelay(topThirdIndicesHoldExRF), binEdges);
%             histogram(rtHoldExRFSortedByTDDelay(bottomThirdIndicesHoldExRF), binEdges);
%             title({'RT distribution of Top and Bottom Third', 'Extreme Target-Dim Delay Firing Rates Hold ExRF'});
%             if rtHoldExRFSortedByTDDelayThirdsPVal < corrStatAlpha
%                 textParams = textParamsBold;
%             else
%                 textParams = textParamsNormal;
%             end
%             text(0.85, 0.95, sprintf('p = %0.3f', rtHoldExRFSortedByTDDelayThirdsPVal), textParams{:});
%             
%             plotFileName = sprintf('%s/%s-%s-rtVsFiringHist-v%d.png', processedDataDir, unitName, blockName, v);
%             fprintf('\tSaving figure to file %s...\n', plotFileName);
%             export_fig(plotFileName, '-nocrop');
            
            % ideally shuffle test?
            % how to deal with correlated increases in firing in a trial?

            rtFiringRateStruct(unitCount).meanRTRelInRFTopThirdFiringRateCTDelay = mean(rtRelInRFSortedByCTDelay(topThirdIndicesRelInRF));
            rtFiringRateStruct(unitCount).meanRTRelInRFBottomThirdFiringRateCTDelay = mean(rtRelInRFSortedByCTDelay(bottomThirdIndicesRelInRF));
            
            rtFiringRateStruct(unitCount).meanRTRelExRFTopThirdFiringRateCTDelay = mean(rtRelExRFSortedByCTDelay(topThirdIndicesRelExRF));
            rtFiringRateStruct(unitCount).meanRTRelExRFBottomThirdFiringRateCTDelay = mean(rtRelExRFSortedByCTDelay(bottomThirdIndicesRelExRF));
            
            rtFiringRateStruct(unitCount).meanRTHoldInRFTopThirdFiringRateCTDelay = mean(rtHoldInRFSortedByCTDelay(topThirdIndicesHoldInRF));
            rtFiringRateStruct(unitCount).meanRTHoldInRFBottomThirdFiringRateCTDelay = mean(rtHoldInRFSortedByCTDelay(bottomThirdIndicesHoldInRF));

            rtFiringRateStruct(unitCount).meanRTHoldExRFTopThirdFiringRateCTDelay = mean(rtHoldExRFSortedByCTDelay(topThirdIndicesHoldExRF));
            rtFiringRateStruct(unitCount).meanRTHoldExRFBottomThirdFiringRateCTDelay = mean(rtHoldExRFSortedByCTDelay(bottomThirdIndicesHoldExRF));
            
            rtFiringRateStruct(unitCount).meanRTHoldInRFTopThirdFiringRateTDDelay = mean(rtHoldInRFSortedByTDDelay(topThirdIndicesHoldInRF));
            rtFiringRateStruct(unitCount).meanRTHoldInRFBottomThirdFiringRateTDDelay = mean(rtHoldInRFSortedByTDDelay(bottomThirdIndicesHoldInRF));

            rtFiringRateStruct(unitCount).meanRTHoldExRFTopThirdFiringRateTDDelay = mean(rtHoldExRFSortedByTDDelay(topThirdIndicesHoldExRF));
            rtFiringRateStruct(unitCount).meanRTHoldExRFBottomThirdFiringRateTDDelay = mean(rtHoldExRFSortedByTDDelay(bottomThirdIndicesHoldExRF));
            
            [rtFiringRateStruct(unitCount).pearsonCorrCoefRelInRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValRelInRFCTDelayRT] = corr(cueTargetDelayRelInRFRate, rtRelInRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefRelExRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValRelExRFCTDelayRT] = corr(cueTargetDelayRelExRFRate, rtRelExRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldInRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldInRFCTDelayRT] = corr(cueTargetDelayHoldInRFRate, rtHoldInRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldExRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldExRFCTDelayRT] = corr(cueTargetDelayHoldExRFRate, rtHoldExRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldInRFTDDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldInRFTDDelayRT] = corr(targetDimDelayHoldInRFRate, rtHoldInRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldExRFTDDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldExRFTDDelayRT] = corr(targetDimDelayHoldExRFRate, rtHoldExRF, 'type', 'Pearson');
            
            [rtFiringRateStruct(unitCount).spearmanCorrCoefRelInRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValRelInRFCTDelayRT] = corr(cueTargetDelayRelInRFRate, rtRelInRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefRelExRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValRelExRFCTDelayRT] = corr(cueTargetDelayRelExRFRate, rtRelExRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldInRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldInRFCTDelayRT] = corr(cueTargetDelayHoldInRFRate, rtHoldInRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldExRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldExRFCTDelayRT] = corr(cueTargetDelayHoldExRFRate, rtHoldExRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldInRFTDDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldInRFTDDelayRT] = corr(targetDimDelayHoldInRFRate, rtHoldInRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldExRFTDDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldExRFTDDelayRT] = corr(targetDimDelayHoldExRFRate, rtHoldExRF, 'type', 'Spearman');
            
            %%
            plotFileName = sprintf('%s/%s-%s-rtVsFiringScatter-v%d.png', processedDataDir, unitName, blockName, v);
            plotRTFiringRateCorrelation(cueTargetDelayRelInRFRate, ...
                    cueTargetDelayRelExRFRate, ...
                    cueTargetDelayHoldInRFRate, ...
                    targetDimDelayHoldInRFRate, ...
                    cueTargetDelayHoldExRFRate, ...
                    targetDimDelayHoldExRFRate, ...
                    rtRelInRF, rtRelExRF, rtHoldInRF, rtHoldExRF, ...
                    rtFiringRateStruct(unitCount), unitName, isZeroDistractors, plotFileName);
            close;
            
            %%
            % note there may be different numbers of trials in each 
            targetDimDelayHoldInRFRateAll{unitCount} = targetDimDelayHoldInRFRate;
            rtHoldInRFAll{unitCount} = rtHoldInRF;
            
            %%            
            inRFPreCueBaseline = ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(inRFLoc);
            exRFPreCueBaseline = ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(exRFLoc);

            % baseline correct the norm factor, account for cases where
            % suppression > enhancement
            % normalize including saccadic responses
            inRFNormFactor = max([ES.maxFiringRateBySpdfInclMotor - inRFPreCueBaseline; inRFPreCueBaseline - ES.minFiringRateBySpdfInclMotor]);
            exRFNormFactor = max([ES.maxFiringRateBySpdfInclMotor - exRFPreCueBaseline; exRFPreCueBaseline - ES.minFiringRateBySpdfInclMotor]);

            spdfInfo.enterFixationSpdfInRFNorm(unitCount,:) = (ES.enterFixation.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.enterFixationSpdfExRFNorm(unitCount,:) = (ES.enterFixation.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.enterFixationSpdfInRFNormErr(unitCount,:) = ES.enterFixation.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.enterFixationSpdfExRFNormErr(unitCount,:) = ES.enterFixation.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            spdfInfo.cueOnsetSpdfInRFNorm(unitCount,:) = (ES.cueOnset.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.cueOnsetSpdfExRFNorm(unitCount,:) = (ES.cueOnset.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.cueOnsetSpdfInRFNormErr(unitCount,:) = ES.cueOnset.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.cueOnsetSpdfExRFNormErr(unitCount,:) = ES.cueOnset.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            spdfInfo.arrayOnsetHoldSpdfInRFNorm(unitCount,:) = (ES.arrayOnsetHoldBal.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.arrayOnsetHoldSpdfExRFNorm(unitCount,:) = (ES.arrayOnsetHoldBal.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.arrayOnsetHoldSpdfInRFNormErr(unitCount,:) = ES.arrayOnsetHoldBal.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.arrayOnsetHoldSpdfExRFNormErr(unitCount,:) = ES.arrayOnsetHoldBal.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            spdfInfo.targetDimSpdfInRFNorm(unitCount,:) = (ES.targetDimBal.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.targetDimSpdfExRFNorm(unitCount,:) = (ES.targetDimBal.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.targetDimSpdfInRFNormErr(unitCount,:) = ES.targetDimBal.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.targetDimSpdfExRFNormErr(unitCount,:) = ES.targetDimBal.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            spdfInfo.exitFixationSpdfInRFNorm(unitCount,:) = (ES.exitFixation.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.exitFixationSpdfExRFNorm(unitCount,:) = (ES.exitFixation.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.exitFixationSpdfInRFNormErr(unitCount,:) = ES.exitFixation.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.exitFixationSpdfExRFNormErr(unitCount,:) = ES.exitFixation.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            spdfInfo.meanNormSpdfInRFAllWindowsAll(unitCount,:) = ([...
                    ES.averageFiringRatesBySpdf.preEnterFixation.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.postEnterFixation.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueResponse.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.arrayResponseHoldBal.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimDelayBal.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimResponseBal.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.preExitFixation.byLoc(inRFLoc)] - inRFPreCueBaseline) / inRFNormFactor;

            spdfInfo.meanNormSpdfExRFAllWindowsAll(unitCount,:) = ([...
                    ES.averageFiringRatesBySpdf.preEnterFixation.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.postEnterFixation.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueResponse.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.arrayResponseHoldBal.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimDelayBal.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimResponseBal.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.preExitFixation.byLoc(exRFLoc)]  - exRFPreCueBaseline) / exRFNormFactor;

            % overwrite for each cell but they should all be the same
            enterFixationT = ES.enterFixation.t - ES.enterFixation.window(1);
            cueOnsetT = ES.cueOnset.t - ES.cueOnset.window(1);
            arrayOnsetT = ES.arrayOnset.t - ES.arrayOnset.window(1);
            targetDimT = ES.targetDimBal.t - ES.targetDimBal.window(1);
            exitFixationT = ES.exitFixation.t - ES.exitFixation.window(1);
            
            % compute time traces with fixed trial time and interpolation
            % between event-locked windows - there's actually no
            % interpolation -- the spdfs are stitched together. should be
            % smooth enough. another option is to concatenate the spike
            % times and then do the smoothing, but that requires the same
            % number of trials for each event lock, which is not the case
            concatCueOnsetWindowOffset = [-0.3 0.3];
            concatArrayOnsetWindowOffset = [-0.3 0.3];
            concatTargetDimWindowOffset = [-0.3 0.15];
            concatExitFixationWindowOffset = [-0.15 0];
            
            concatCueOnsetWindowIndices = getTimeLogicalWithTolerance(cueOnsetT, concatCueOnsetWindowOffset);
            concatArrayOnsetWindowIndices = getTimeLogicalWithTolerance(arrayOnsetT, concatArrayOnsetWindowOffset);
            concatTargetDimWindowIndices = getTimeLogicalWithTolerance(targetDimT, concatTargetDimWindowOffset);
            concatExitFixationWindowIndices = getTimeLogicalWithTolerance(exitFixationT, concatExitFixationWindowOffset);
            nConcatValues = sum(concatCueOnsetWindowIndices) + ...
                    sum(concatArrayOnsetWindowIndices) + ...
                    sum(concatTargetDimWindowIndices) + ...
                    sum(concatExitFixationWindowIndices);
            tStep = mean(diff(cueOnsetT));
            assert(all((tStep - diff(cueOnsetT)) < 1e-8));
            tStart = cueOnsetT(find(concatCueOnsetWindowIndices, 1, 'first'));
            tEnd = tStart + (nConcatValues - 1) * tStep;
            spdfInfo.concatCueOnsetT = tStart:tStep:tEnd;
            
            spdfInfo.meanSpdfInRFConcatAll(unitCount,:) = [...
                    ES.cueOnset.spdfByLoc(inRFLoc,concatCueOnsetWindowIndices) ...
                    ES.arrayOnset.spdfByLoc(inRFLoc,concatArrayOnsetWindowIndices) ...
                    ES.targetDimBal.spdfByLoc(inRFLoc,concatTargetDimWindowIndices) ...
                    ES.exitFixation.spdfByLoc(inRFLoc,concatExitFixationWindowIndices)];
                
            spdfInfo.meanSpdfExRFConcatAll(unitCount,:) = [...
                    ES.cueOnset.spdfByLoc(exRFLoc,concatCueOnsetWindowIndices) ...
                    ES.arrayOnset.spdfByLoc(exRFLoc,concatArrayOnsetWindowIndices) ...
                    ES.targetDimBal.spdfByLoc(exRFLoc,concatTargetDimWindowIndices) ...
                    ES.exitFixation.spdfByLoc(exRFLoc,concatExitFixationWindowIndices)];
            
            for k = 1:numel(ES.isLocUsed)
                if ES.isLocUsed(k)
                    if unitCount > 1
                        arrayOnsetRelBalSpikeTimesByLoc{k}(unitCount,:) = ES.arrayOnsetRelBal.spikeTimesByLoc{k};
                        arrayOnsetHoldBalSpikeTimesByLoc{k}(unitCount,:) = ES.arrayOnsetHoldBal.spikeTimesByLoc{k};
                        targetDimBalSpikeTimesByLoc{k}(unitCount,:) = ES.targetDimBal.spikeTimesByLoc{k};
                    else
                        arrayOnsetRelBalSpikeTimesByLoc{k} = ES.arrayOnsetRelBal.spikeTimesByLoc{k};
                        arrayOnsetHoldBalSpikeTimesByLoc{k} = ES.arrayOnsetHoldBal.spikeTimesByLoc{k};
                        targetDimBalSpikeTimesByLoc{k} = ES.targetDimBal.spikeTimesByLoc{k};
                    end
                end
            end
            inRFLocs(unitCount) = inRFLoc;
            exRFLocs(unitCount) = exRFLoc;
            
            % could use peak method or resampled baseline method
            arrayHoldBalLatencyInRF(unitCount) = ES.arrayOnsetHoldBal.latencyByLoc(inRFLoc);
            arrayHoldBalLatencyExRF(unitCount) = ES.arrayOnsetHoldBal.latencyByLoc(exRFLoc);
            
            targetDimBalLatencyInRF(unitCount) = ES.targetDimBal.latencyByLoc(inRFLoc);
            targetDimBalLatencyExRF(unitCount) = ES.targetDimBal.latencyByLoc(exRFLoc);
            
            fn = fieldnames(ES.averageFiringRatesBySpdf);
            for k = 1:numel(fn)
                if isfield(averageFiringRatesBySpdf, fn{k})
                    averageFiringRatesBySpdf.(fn{k})(unitCount,:) = [ES.averageFiringRatesBySpdf.(fn{k})]';
                else
                    averageFiringRatesBySpdf.(fn{k}) = [ES.averageFiringRatesBySpdf.(fn{k})]'; % no pre-allocation
                end
            end
            
            fn = fieldnames(ES.averageFiringRatesByCount);
            for k = 1:numel(fn)
                if isfield(averageFiringRatesByCount, fn{k})
                    averageFiringRatesByCount.(fn{k})(unitCount,:) = [ES.averageFiringRatesByCount.(fn{k})]';
                else
                    averageFiringRatesByCount.(fn{k}) = [ES.averageFiringRatesByCount.(fn{k})]'; % no pre-allocation
                end
            end

        end
    else
        fprintf('Skipping %s due to min firing rate requirement...\n', unitName);
    end
end

%% plot pre-saccade slopes
figure_tr_inch(15, 5);
subaxis(1, 3, 1);
plot(earlyPreExitFixationSlope, latePreExitFixationSlope, '.', 'MarkerSize', 20);
xlabel('Early pre-saccade slope');
ylabel('Late pre-saccade slope');

subaxis(1, 3, 2);
bar(muaChannelsToLoad, earlyPreExitFixationSlope);
xlabel('Early pre-saccade slope');
ylim([-300 300]);

subaxis(1, 3, 3);
bar(muaChannelsToLoad, latePreExitFixationSlope);
xlabel('Late pre-saccade slope');
ylim([-300 300]);

plotFileName = sprintf('%s/%s-sessionInd%d-preSaccadicSlopes-v%d.png', outputDir, sessionName, sessionInd, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% compute single trial population latency by combining spikes across recordings on a probe
% correlate population latency with RT
rtRelBalByLoc = cell(nLoc, 1);
rtHoldBalByLoc = cell(nLoc, 1);
for k = 1:nLoc
    rtRelBalByLoc{k} = ES.UE.rt(ES.UE.cueLoc == k & ES.UE.isRelBal);
    rtHoldBalByLoc{k} = ES.UE.rt(ES.UE.cueLoc == k & ES.UE.isHoldBal);
end
plotFileName = sprintf('%s/%s-sessionInd%d-rtVsArrayOnsetRelLatency-v%d.png', outputDir, sessionName, sessionInd, v);
plotRTLatencyCorrelation(ES.arrayOnsetRelBal, ES.isLocUsed, rtRelBalByLoc, arrayOnsetRelBalSpikeTimesByLoc, inRFLocs, plotFileName);
plotFileName = sprintf('%s/%s-sessionInd%d-rtVsTargetDimLatency-v%d.png', outputDir, sessionName, sessionInd, v);
plotRTLatencyCorrelation(ES.targetDimBal, ES.isLocUsed, rtHoldBalByLoc, targetDimBalSpikeTimesByLoc, inRFLocs, plotFileName);

%% correlate target dim delay activity (rank trial within unit) with rt (rank trial)
% not enough variance within target dim delay activity across trials
% targetDimDelayHoldInRFRankAll = [];
% rtHoldInRFRankAll = [];
% figure;
% hold on;
% for i = 1:unitCount
%     assert(numel(targetDimDelayHoldInRFRateAll{i}) == numel(rtHoldInRFAll{i}));
%     nTrialInRF = numel(rtHoldInRFAll{i});
%     if mean(targetDimDelayHoldInRFRateAll{i}) > 10
%         targetDimDelayHoldInRFRankAll = [targetDimDelayHoldInRFRankAll; tiedrank(targetDimDelayHoldInRFRateAll{i}) / nTrialInRF];
%         rtHoldInRFRankAll = [rtHoldInRFRankAll; rtHoldInRFAll{i}];
%         plot(tiedrank(targetDimDelayHoldInRFRateAll{i}) / nTrialInRF, rtHoldInRFAll{i}, '.');
%     end
% end

%% save
saveFileName = sprintf('%s/%s-sessionInd%d-muaAnalysisSummaryData-v%d.mat', outputDir, sessionName, sessionInd, v);
fprintf('\n');
fprintf('Writing file %s ...\n', saveFileName);
save(saveFileName, ...
        'unitNames', ...
        'isSignificantResponseVsBaseline', ...
        'isSignificantResponseVsPreviousPeriod', ...
        'isSignificantResponseVsBootstrapBaseline', ...
        'isSignificantResponseVsBootstrapPreviousPeriod', ...
        'isSignificantSelectivity', ...
        'cueResponseVsBaselineDirection', ...
        'infoRates', ...
        'diffRates', ...
        'attnIndices', ...
        'localization', ...
        'isInVPulvinar', ...
        'isInDPulvinar', ...
        'preCueBaselineExpFit', ...
        'earlyPreExitFixationSlope', ...
        'latePreExitFixationSlope', ...
        'spdfInfo', ...
        'enterFixationT', ...
        'cueOnsetT', ...
        'arrayOnsetT', ...
        'targetDimT', ...
        'exitFixationT', ...
        'rtFiringRateStruct', ...
        'arrayHoldBalLatencyInRF', ...
        'arrayHoldBalLatencyExRF', ...
        'targetDimBalLatencyInRF', ...
        'targetDimBalLatencyExRF', ...
        'averageFiringRatesBySpdf', ...
        'averageFiringRatesByCount', ...
        'inRFLocs', ...
        'exRFLocs');
