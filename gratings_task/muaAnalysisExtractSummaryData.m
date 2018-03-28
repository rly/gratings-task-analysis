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
cueResponseVsBootstrapBaselineDirection = zeros(nUnitsApprox, 1);
preExitFixationVsBootstrapBaselineDirection = zeros(nUnitsApprox, 1);
infoRates = nan(nUnitsApprox, 5); % 5 periods
diffRates = nan(nUnitsApprox, 3); % 2 delay periods + array response
attnIndices = nan(nUnitsApprox, 3); % 2 delay periods + array response
localization = cell(nUnitsApprox, 1);
isInVPulvinar = false(nUnitsApprox, 1);
isInDPulvinar = false(nUnitsApprox, 1);
earlyPreExitFixationSlope = nan(nUnitsApprox, 1);
latePreExitFixationSlope = nan(nUnitsApprox, 1);
earlyPreExitFixationWindowOffset = [-0.2 -0.05];
latePreExitFixationWindowOffset = [-0.05 0];

rtFiringRateStruct = struct();

nLoc = 4;
arrayOnsetRelSpikeTimesByLoc = cell(nLoc, 1);
arrayOnsetHoldSpikeTimesByLoc = cell(nLoc, 1);
targetDimSpikeTimesByLoc = cell(nLoc, 1);
inRFLocs = nan(nUnitsApprox, 1);
exRFLocs = nan(nUnitsApprox, 1);

arrayHoldResponseLatencyInRF = nan(nUnitsApprox, 1);
arrayHoldResponseLatencyExRF = nan(nUnitsApprox, 1);
arrayResponseLatencyInRF = nan(nUnitsApprox, 1);
arrayResponseLatencyExRF = nan(nUnitsApprox, 1);
targetDimResponseLatencyInRF = nan(nUnitsApprox, 1);
targetDimResponseLatencyExRF = nan(nUnitsApprox, 1);

averageFiringRatesBySpdf = struct();
averageFiringRatesByCount = struct();

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
                any(ES.averageFiringRatesBySpdf.arrayHoldResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.arrayRelResponse.byLoc >= minFiringRate) || ...
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
                    min([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.preExitFixationVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLocUsed];
                
            % significance test using rank sum test. correct for multiple
            % comparisons by comparing p < statAlpha / nLocUsed
            isSignificantResponseVsPreviousPeriod(unitCount,:) = [...
                    min([ES.cueResponseVsBaselineSignRankTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.arrayHoldResponseVsCueTargetDelaySignRankTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.targetDimResponseVsTargetDimDelaySignRankTestStatsByLoc.p]) < statAlpha / nLocUsed ...
                    min([ES.preExitFixationVsPreExitFixationEarlySignRankTestStatsByLoc.p]) < statAlpha / nLocUsed];
            
            % significance test using non-parametric test vs bootstrapped
            % baseline. correct for multiple comparisons by comparing p <
            % statAlpha / nLocUsed
            isSignificantResponseVsBootstrapBaseline(unitCount,:) = [...
                    min(ES.cueResponsePValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.cueTargetDelayPValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.targetDimDelayPValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.targetDimResponsePValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.preExitFixationPValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed];
                
            isSignificantResponseVsBootstrapPreviousPeriod(unitCount,:) = [...
                    min(ES.cueResponsePValueByBootstrapBaselineSpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.arrayHoldResponsePValueByBootstrapCueTargetDelaySpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc) < statAlpha / nLocUsed ...
                    min(ES.preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc) < statAlpha / nLocUsed];
            
            
            % look at most significant response which is not necessarily
            % InRF location
            [~,loc] = min(ES.cueResponsePValueByBootstrapBaselineSpdfByLoc);
            if ES.averageFiringRatesBySpdf.cueResponse.byLoc(loc) > ES.meanBootstrappedMeanPreCueBaselines
                cueResponseVsBootstrapBaselineDirection(unitCount) = 1;
            else
                cueResponseVsBootstrapBaselineDirection(unitCount) = -1;
            end
            [~,loc] = min(ES.preExitFixationPValueByBootstrapBaselineSpdfByLoc);
            if ES.averageFiringRatesBySpdf.preExitFixation.byLoc(loc) > ES.meanBootstrappedMeanPreCueBaselines
                preExitFixationVsBootstrapBaselineDirection(unitCount) = 1;
            else
                preExitFixationVsBootstrapBaselineDirection(unitCount) = -1;
            end
            clear loc;
            
            % determine slope of firing rate prior to exit fixation saccade
            earlyPreExitFixationSlopeStruct = computeSlopeFiringRateBySpdf(earlyPreExitFixationWindowOffset, ES.exitFixation);
            latePreExitFixationSlopeStruct = computeSlopeFiringRateBySpdf(latePreExitFixationWindowOffset, ES.exitFixation);
            earlyPreExitFixationSlope(unitCount) = earlyPreExitFixationSlopeStruct.all;
            latePreExitFixationSlope(unitCount) = latePreExitFixationSlopeStruct.all;
            
            
            isSignificantSelectivity(unitCount,:) = [...
                    ES.cueResponseInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.cueTargetDelayInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.arrayHoldResponseInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.targetDimDelayInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.targetDimResponseInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha];

            % correct for multiple comparisons: 4 loc x 5 periods
%                 baseAlpha = 0.05;
%                 isSignificantAnyTaskMod = isCell(unitCount) && ...
%                         min([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) < baseAlpha / numel([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) / 5 || ...
%                         min([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) < baseAlpha / numel([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) / 5 || ...
%                         min([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) < baseAlpha / numel([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) / 5 || ...
%                         min([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) < baseAlpha / numel([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) / 5 || ...
%                         min([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) < baseAlpha / numel([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) / 5;

            infoRates(unitCount,:) = [...
                    ES.cueResponseInfoRateStruct.infoRate ...
                    ES.cueTargetDelayInfoRateStruct.infoRate ...
                    ES.arrayHoldResponseInfoRateStruct.infoRate ...
                    ES.targetDimDelayInfoRateStruct.infoRate ...
                    ES.targetDimResponseInfoRateStruct.infoRate];

            diffRates(unitCount,:) = [...
                    ES.cueTargetDelayDiff ...
                    ES.arrayHoldResponseDiff ...
                    ES.targetDimDelayDiff];    
            
            attnIndices(unitCount,:) = [...
                    ES.cueTargetDelayAI ...
                    ES.arrayHoldResponseAI ...
                    ES.targetDimDelayAI];

            inRFLoc = ES.inRFLoc;
            exRFLoc = ES.exRFLoc;
            
            %%
            % extract firing rates (count-based) at InRF and ExRF locations
            % note: cueTargetDelayLongWindowOffset = [-0.4 0];
            % note: targetDimDelayLongWindowOffset = [-0.4 0];
            
            rtRelInRF = ES.UE.rt(ES.UE.cueLoc == inRFLoc & ~ES.UE.isHoldTrial);
            rtRelExRF = ES.UE.rt(ES.UE.cueLoc == exRFLoc & ~ES.UE.isHoldTrial);
            rtHoldInRF = ES.UE.rt(ES.UE.cueLoc == inRFLoc & ES.UE.isHoldTrial);
            rtHoldExRF = ES.UE.rt(ES.UE.cueLoc == exRFLoc & ES.UE.isHoldTrial);

            cueTargetDelayLongRelInRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongRel.trialRateByLoc{inRFLoc};
            cueTargetDelayLongRelExRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongRel.trialRateByLoc{exRFLoc};
            
            cueTargetDelayLongHoldInRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongHold.trialRateByLoc{inRFLoc};
            targetDimDelayLongHoldInRFRate = ES.averageFiringRatesByCount.targetDimDelayLong.trialRateByLoc{inRFLoc};

            cueTargetDelayLongHoldExRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongHold.trialRateByLoc{exRFLoc};
            targetDimDelayLongHoldExRFRate = ES.averageFiringRatesByCount.targetDimDelayLong.trialRateByLoc{exRFLoc};

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
            
            assert(nTrialsRelInRF == numel(cueTargetDelayLongRelInRFRate));
            assert(nTrialsRelExRF == numel(cueTargetDelayLongRelExRFRate));
            assert(all(nTrialsHoldInRF == [numel(cueTargetDelayLongHoldInRFRate) numel(targetDimDelayLongHoldInRFRate)]));
            assert(all(nTrialsHoldExRF == [numel(cueTargetDelayLongHoldExRFRate) numel(targetDimDelayLongHoldExRFRate)]));
            
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
            [~,sortCTDelayRelInRFInd] = sortBreakOrder(cueTargetDelayLongRelInRFRate);
            rtRelInRFSortedByCTDelay = rtRelInRF(sortCTDelayRelInRFInd);
            
            [~,sortCTDelayRelExRFInd] = sortBreakOrder(cueTargetDelayLongRelExRFRate);
            rtRelExRFSortedByCTDelay = rtRelExRF(sortCTDelayRelExRFInd);
            
            [~,sortCTDelayHoldInRFInd] = sortBreakOrder(cueTargetDelayLongHoldInRFRate);
            rtHoldInRFSortedByCTDelay = rtHoldInRF(sortCTDelayHoldInRFInd);
            
            [~,sortCTDelayHoldExRFInd] = sortBreakOrder(cueTargetDelayLongHoldExRFRate);
            rtHoldExRFSortedByCTDelay = rtHoldExRF(sortCTDelayHoldExRFInd);
            
            [~,sortTDDelayHoldInRFInd] = sortBreakOrder(targetDimDelayLongHoldInRFRate);
            rtHoldInRFSortedByTDDelay = rtHoldInRF(sortTDDelayHoldInRFInd);
            
            [~,sortTDDelayHoldExRFInd] = sortBreakOrder(targetDimDelayLongHoldExRFRate);
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
            
            [rtFiringRateStruct(unitCount).pearsonCorrCoefRelInRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValRelInRFCTDelayRT] = corr(cueTargetDelayLongRelInRFRate, rtRelInRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefRelExRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValRelExRFCTDelayRT] = corr(cueTargetDelayLongRelExRFRate, rtRelExRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldInRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldInRFCTDelayRT] = corr(cueTargetDelayLongHoldInRFRate, rtHoldInRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldExRFCTDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldExRFCTDelayRT] = corr(cueTargetDelayLongHoldExRFRate, rtHoldExRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldInRFTDDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldInRFTDDelayRT] = corr(targetDimDelayLongHoldInRFRate, rtHoldInRF, 'type', 'Pearson');
            [rtFiringRateStruct(unitCount).pearsonCorrCoefHoldExRFTDDelayRT,rtFiringRateStruct(unitCount).pearsonCorrCoefPValHoldExRFTDDelayRT] = corr(targetDimDelayLongHoldExRFRate, rtHoldExRF, 'type', 'Pearson');
            
            [rtFiringRateStruct(unitCount).spearmanCorrCoefRelInRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValRelInRFCTDelayRT] = corr(cueTargetDelayLongRelInRFRate, rtRelInRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefRelExRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValRelExRFCTDelayRT] = corr(cueTargetDelayLongRelExRFRate, rtRelExRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldInRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldInRFCTDelayRT] = corr(cueTargetDelayLongHoldInRFRate, rtHoldInRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldExRFCTDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldExRFCTDelayRT] = corr(cueTargetDelayLongHoldExRFRate, rtHoldExRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldInRFTDDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldInRFTDDelayRT] = corr(targetDimDelayLongHoldInRFRate, rtHoldInRF, 'type', 'Spearman');
            [rtFiringRateStruct(unitCount).spearmanCorrCoefHoldExRFTDDelayRT,rtFiringRateStruct(unitCount).spearmanCorrCoefPValHoldExRFTDDelayRT] = corr(targetDimDelayLongHoldExRFRate, rtHoldExRF, 'type', 'Spearman');
            
            %%
            plotFileName = sprintf('%s/%s-%s-rtVsFiringScatter-v%d.png', processedDataDir, unitName, blockName, v);
            plotRTFiringRateCorrelation(cueTargetDelayLongRelInRFRate, ...
                    cueTargetDelayLongRelExRFRate, ...
                    cueTargetDelayLongHoldInRFRate, ...
                    targetDimDelayLongHoldInRFRate, ...
                    cueTargetDelayLongHoldExRFRate, ...
                    targetDimDelayLongHoldExRFRate, ...
                    rtRelInRF, rtRelExRF, rtHoldInRF, rtHoldExRF, ...
                    rtFiringRateStruct(unitCount), unitName, isZeroDistractors, plotFileName);
            close;
            
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

            spdfInfo.arrayOnsetHoldSpdfInRFNorm(unitCount,:) = (ES.arrayOnsetHold.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.arrayOnsetHoldSpdfExRFNorm(unitCount,:) = (ES.arrayOnsetHold.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.arrayOnsetHoldSpdfInRFNormErr(unitCount,:) = ES.arrayOnsetHold.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.arrayOnsetHoldSpdfExRFNormErr(unitCount,:) = ES.arrayOnsetHold.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            spdfInfo.targetDimSpdfInRFNorm(unitCount,:) = (ES.targetDim.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            spdfInfo.targetDimSpdfExRFNorm(unitCount,:) = (ES.targetDim.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            spdfInfo.targetDimSpdfInRFNormErr(unitCount,:) = ES.targetDim.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            spdfInfo.targetDimSpdfExRFNormErr(unitCount,:) = ES.targetDim.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

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
                    ES.averageFiringRatesBySpdf.arrayHoldResponse.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimResponse.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.preExitFixation.byLoc(inRFLoc)] - inRFPreCueBaseline) / inRFNormFactor;

            spdfInfo.meanNormSpdfExRFAllWindowsAll(unitCount,:) = ([...
                    ES.averageFiringRatesBySpdf.preEnterFixation.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.postEnterFixation.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueResponse.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.arrayHoldResponse.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimResponse.byLoc(exRFLoc) ...
                    ES.averageFiringRatesBySpdf.preExitFixation.byLoc(exRFLoc)]  - exRFPreCueBaseline) / exRFNormFactor;

            % overwrite for each cell but they should all be the same
            enterFixationT = ES.enterFixation.t - ES.enterFixation.window(1);
            cueOnsetT = ES.cueOnset.t - ES.cueOnset.window(1);
            arrayOnsetT = ES.arrayOnset.t - ES.arrayOnset.window(1);
            targetDimT = ES.targetDim.t - ES.targetDim.window(1);
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
                    ES.targetDim.spdfByLoc(inRFLoc,concatTargetDimWindowIndices) ...
                    ES.exitFixation.spdfByLoc(inRFLoc,concatExitFixationWindowIndices)];
                
            spdfInfo.meanSpdfExRFConcatAll(unitCount,:) = [...
                    ES.cueOnset.spdfByLoc(exRFLoc,concatCueOnsetWindowIndices) ...
                    ES.arrayOnset.spdfByLoc(exRFLoc,concatArrayOnsetWindowIndices) ...
                    ES.targetDim.spdfByLoc(exRFLoc,concatTargetDimWindowIndices) ...
                    ES.exitFixation.spdfByLoc(exRFLoc,concatExitFixationWindowIndices)];
            
            for k = 1:numel(ES.isLocUsed)
                if ES.isLocUsed(k)
                    if unitCount > 1
                        arrayOnsetRelSpikeTimesByLoc{k}(unitCount,:) = ES.arrayOnsetRel.spikeTimesByLoc{k};
                        arrayOnsetHoldSpikeTimesByLoc{k}(unitCount,:) = ES.arrayOnsetHold.spikeTimesByLoc{k};
                        targetDimSpikeTimesByLoc{k}(unitCount,:) = ES.targetDim.spikeTimesByLoc{k};
                    else
                        arrayOnsetRelSpikeTimesByLoc{k} = ES.arrayOnsetRel.spikeTimesByLoc{k};
                        arrayOnsetHoldSpikeTimesByLoc{k} = ES.arrayOnsetHold.spikeTimesByLoc{k};
                        targetDimSpikeTimesByLoc{k} = ES.targetDim.spikeTimesByLoc{k};
                    end
                end
            end
            inRFLocs(unitCount) = inRFLoc;
            exRFLocs(unitCount) = exRFLoc;
            
            % could use peak method or resampled baseline method
            arrayHoldResponseLatencyInRF(unitCount) = ES.arrayOnsetHold.latencyBootByLoc(inRFLoc);
            arrayHoldResponseLatencyExRF(unitCount) = ES.arrayOnsetHold.latencyBootByLoc(exRFLoc);
            
            % could use peak method or resampled baseline method
            arrayResponseLatencyInRF(unitCount) = ES.arrayOnset.latencyBootByLoc(inRFLoc);
            arrayResponseLatencyExRF(unitCount) = ES.arrayOnset.latencyBootByLoc(exRFLoc);
            
            targetDimResponseLatencyInRF(unitCount) = ES.targetDim.latencyBootByLoc(inRFLoc);
            targetDimResponseLatencyExRF(unitCount) = ES.targetDim.latencyBootByLoc(exRFLoc);
            
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

%% single trial population latency relationship with RT

%% compute single trial population latency by combining spikes across recordings on a probe
rtRelByLoc = cell(nLoc, 1);
rtHoldByLoc = cell(nLoc, 1);
for k = 1:nLoc
    rtRelByLoc{k} = ES.UE.rt(ES.UE.cueLoc == k & ~ES.UE.isHoldTrial);
    rtHoldByLoc{k} = ES.UE.rt(ES.UE.cueLoc == k & ES.UE.isHoldTrial);
end
plotFileName = sprintf('%s/%s-sessionInd%d-rtVsArrayOnsetRelLatency-v%d.png', outputDir, sessionName, sessionInd, v);
plotRTLatencyCorrelation(ES.arrayOnset, ES.isLocUsed, rtRelByLoc, arrayOnsetRelSpikeTimesByLoc, inRFLocs, plotFileName);
plotFileName = sprintf('%s/%s-sessionInd%d-rtVsTargetDimLatency-v%d.png', outputDir, sessionName, sessionInd, v);
plotRTLatencyCorrelation(ES.targetDim, ES.isLocUsed, rtHoldByLoc, targetDimSpikeTimesByLoc, inRFLocs, plotFileName);

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
        'cueResponseVsBootstrapBaselineDirection', ...
        'preExitFixationVsBootstrapBaselineDirection', ...
        'infoRates', ...
        'diffRates', ...
        'attnIndices', ...
        'localization', ...
        'isInVPulvinar', ...
        'isInDPulvinar', ...
        'earlyPreExitFixationSlope', ...
        'latePreExitFixationSlope', ...
        'spdfInfo', ...
        'enterFixationT', ...
        'cueOnsetT', ...
        'arrayOnsetT', ...
        'targetDimT', ...
        'exitFixationT', ...
        'rtFiringRateStruct', ...
        'rtRelInRF', ...
        'rtRelExRF', ...
        'rtHoldInRF', ...
        'rtHoldExRF', ...
        'arrayHoldResponseLatencyInRF', ...
        'arrayHoldResponseLatencyExRF', ...
        'arrayResponseLatencyInRF', ...
        'arrayResponseLatencyExRF', ...
        'targetDimResponseLatencyInRF', ...
        'targetDimResponseLatencyExRF', ...
        'averageFiringRatesBySpdf', ...
        'averageFiringRatesByCount', ...
        'inRFLocs', ...
        'exRFLocs');
