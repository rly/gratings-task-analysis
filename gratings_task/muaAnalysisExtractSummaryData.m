function muaAnalysisExtractSummaryData(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)

v = 10;
tic;

% for preallocation. make sure this is an underestimate or equal to actual
% number of units saved
nUnitsApprox = 1; 

unitNames = cell(nUnitsApprox, 1);
isSignificantStats = false(nUnitsApprox, 10); % 5 periods > baseline, 5 periods info rate
infoRates = nan(nUnitsApprox, 5); % 5 periods
attnIndices = nan(nUnitsApprox, 2); % 2 delay periods
localization = cell(nUnitsApprox, 1);
isInVPulvinar = false(nUnitsApprox, 1);
meanRTHoldInRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldInRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelInRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelInRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldInRFTopThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
meanRTHoldInRFBottomThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelExRFTopThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTRelExRFBottomThirdFiringRateCTDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFTopThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
meanRTHoldExRFBottomThirdFiringRateTDDelayAll = nan(nUnitsApprox, 1);
corrCoefHoldInRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefRelInRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefHoldInRFTDDelayRTAll = nan(nUnitsApprox, 3);
corrCoefHoldExRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefRelExRFCTDelayRTAll = nan(nUnitsApprox, 3);
corrCoefHoldExRFTDDelayRTAll = nan(nUnitsApprox, 3);

nTimeEnterFixation = 1401; % hard coded temp
enterFixationSpdfInRFNorm = nan(nUnitsApprox, nTimeEnterFixation);
enterFixationSpdfExRFNorm = nan(nUnitsApprox, nTimeEnterFixation);
enterFixationSpdfInRFNormErr = nan(nUnitsApprox, nTimeEnterFixation);
enterFixationSpdfExRFNormErr = nan(nUnitsApprox, nTimeEnterFixation);
nTimeCueOnset = 1401; % hard coded temp
cueOnsetSpdfInRFNorm = nan(nUnitsApprox, nTimeCueOnset);
cueOnsetSpdfExRFNorm = nan(nUnitsApprox, nTimeCueOnset);
cueOnsetSpdfInRFNormErr = nan(nUnitsApprox, nTimeCueOnset);
cueOnsetSpdfExRFNormErr = nan(nUnitsApprox, nTimeCueOnset);
nTimeArrayOnset = 1401; % hard coded temp
arrayOnsetHoldSpdfInRFNorm = nan(nUnitsApprox, nTimeArrayOnset);
arrayOnsetHoldSpdfExRFNorm = nan(nUnitsApprox, nTimeArrayOnset);
arrayOnsetHoldSpdfInRFNormErr = nan(nUnitsApprox, nTimeArrayOnset);
arrayOnsetHoldSpdfExRFNormErr = nan(nUnitsApprox, nTimeArrayOnset);
nTimeTargetDim = 1401; % hard coded temp
targetDimSpdfInRFNorm = nan(nUnitsApprox, nTimeTargetDim);
targetDimSpdfExRFNorm = nan(nUnitsApprox, nTimeTargetDim);
targetDimSpdfInRFNormErr = nan(nUnitsApprox, nTimeTargetDim);
targetDimSpdfExRFNormErr = nan(nUnitsApprox, nTimeTargetDim);
nTimeExitFixation = 1401; % hard coded temp
exitFixationSpdfInRFNorm = nan(nUnitsApprox, nTimeExitFixation);
exitFixationSpdfExRFNorm = nan(nUnitsApprox, nTimeExitFixation);
exitFixationSpdfInRFNormErr = nan(nUnitsApprox, nTimeExitFixation);
exitFixationSpdfExRFNormErr = nan(nUnitsApprox, nTimeExitFixation);

% N time-locking periods
meanNormSpdfInRFAllWindowsAll = nan(nUnitsApprox, 9);
meanNormSpdfExRFAllWindowsAll = nan(nUnitsApprox, 9);

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
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, muaChannelsToLoad, 'Gratings', 'Gratings', 1, 0);
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

            if ismember(spikeStruct.channelID, R.pldChannels)
                localization{unitCount} = 'PLd';
            elseif ismember(spikeStruct.channelID, R.plvChannels)
                localization{unitCount} = 'PLv';
            elseif ismember(spikeStruct.channelID, R.pmChannels)
                localization{unitCount} = 'PM';
            elseif ismember(spikeStruct.channelID, R.piChannels)
                localization{unitCount} = 'PI';
            else
                localization{unitCount} = '';
            end

            if ismember(spikeStruct.channelID, R.vPulChannels)
                isInVPulvinar(unitCount) = 1;
                localization{unitCount} = 'vPul'; % TEMP
            else
                isInVPulvinar(unitCount) = 0;
            end

            assert(all(numel([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) == ...
                    [numel([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p])
                    numel([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p])
                    numel([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p])
                    numel([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p])]));
            nLoc = numel([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]);

            isSignificantStats(unitCount,:) = [...
                    min([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLoc ...
                    min([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLoc ...
                    min([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLoc ...
                    min([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLoc ...
                    min([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / nLoc ...
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

            attnIndices(unitCount,:) = [...
                    ES.cueTargetDelayAI ...
                    ES.targetDimDelayAI];

            inRFLoc = ES.inRFLoc;
            exRFLoc = ES.exRFLoc;
            
            
            

            cueTargetDelayLongHoldInRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongHold.trialRateByLoc{inRFLoc};
            targetDimDelayLongHoldInRFRate = ES.averageFiringRatesByCount.targetDimDelayLong.trialRateByLoc{inRFLoc};
            rtHoldInRF = ES.UE.rt(ES.UE.cueLoc == inRFLoc & ES.UE.isHoldTrial);

            cueTargetDelayLongHoldExRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongHold.trialRateByLoc{exRFLoc};
            targetDimDelayLongHoldExRFRate = ES.averageFiringRatesByCount.targetDimDelayLong.trialRateByLoc{exRFLoc};
            rtHoldExRF = ES.UE.rt(ES.UE.cueLoc == exRFLoc & ES.UE.isHoldTrial);

            cueTargetDelayLongRelInRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongRel.trialRateByLoc{inRFLoc};
            rtRelInRF = ES.UE.rt(ES.UE.cueLoc == inRFLoc & ~ES.UE.isHoldTrial);

            cueTargetDelayLongRelExRFRate = ES.averageFiringRatesByCount.cueTargetDelayLongRel.trialRateByLoc{exRFLoc};
            rtRelExRF = ES.UE.rt(ES.UE.cueLoc == exRFLoc & ~ES.UE.isHoldTrial);

%                 figure_tr_inch(5, 10);
%                 subaxis(2, 1, 1);
%                 hold on;
%                 plot(cueTargetDelayLongHoldInRFRate, rtHoldInRF, ...
%                         '.', 'MarkerSize', 10);
%                 plot(cueTargetDelayLongHoldExRFRate, rtHoldExRF, ...
%                         '.', 'MarkerSize', 10);
%                     
%                 subaxis(2, 1, 2);
%                 hold on;
%                 plot(targetDimDelayLongHoldInRFRate, rtHoldInRF, ...
%                         '.', 'MarkerSize', 10);
%                 plot(targetDimDelayLongHoldExRFRate, rtHoldExRF, ...
%                         '.', 'MarkerSize', 10);
%                 pause;
%                 close;

            % note: most of these values are 0 because of the short
            % time window and the sparse firing. even with using
            % sortBreakOrder, this can lead to problems
            [~,sortCTDelayHoldInRFInd] = sortBreakOrder(cueTargetDelayLongHoldInRFRate);
            [~,sortCTDelayRelInRFInd] = sortBreakOrder(cueTargetDelayLongRelInRFRate);
            [~,sortTDDelayHoldInRFInd] = sortBreakOrder(targetDimDelayLongHoldInRFRate);
            rtHoldInRFSortedByCTDelay = rtHoldInRF(sortCTDelayHoldInRFInd);
            rtRelInRFSortedByCTDelay = rtRelInRF(sortCTDelayRelInRFInd);
            rtHoldInRFSortedByTDDelay = rtHoldInRF(sortTDDelayHoldInRFInd);

            [~,sortCTDelayHoldExRFInd] = sortBreakOrder(cueTargetDelayLongHoldExRFRate);
            [~,sortCTDelayRelExRFInd] = sortBreakOrder(cueTargetDelayLongRelExRFRate);
            [~,sortTDDelayHoldExRFInd] = sortBreakOrder(targetDimDelayLongHoldExRFRate);
            rtHoldExRFSortedByCTDelay = rtHoldExRF(sortCTDelayHoldExRFInd);
            rtRelExRFSortedByCTDelay = rtRelExRF(sortCTDelayRelExRFInd);
            rtHoldExRFSortedByTDDelay = rtHoldExRF(sortTDDelayHoldExRFInd);

            nTrialsHoldInRF = numel(cueTargetDelayLongHoldInRFRate);
            topThirdIndicesHoldInRF = round(nTrialsHoldInRF * 2/3)+1:nTrialsHoldInRF;
            bottomThirdIndicesHoldInRF = 1:round(nTrialsHoldInRF * 1/3);

            nTrialsRelInRF = numel(cueTargetDelayLongRelInRFRate);
            topThirdIndicesRelInRF = round(nTrialsRelInRF * 2/3)+1:nTrialsRelInRF;
            bottomThirdIndicesRelInRF = 1:round(nTrialsRelInRF * 1/3);

            nTrialsHoldExRF = numel(cueTargetDelayLongHoldExRFRate);
            topThirdIndicesHoldExRF = round(nTrialsHoldExRF * 2/3)+1:nTrialsHoldExRF;
            bottomThirdIndicesHoldExRF = 1:round(nTrialsHoldExRF * 1/3);

            nTrialsRelExRF = numel(cueTargetDelayLongRelExRFRate);
            topThirdIndicesRelExRF = round(nTrialsRelExRF * 2/3)+1:nTrialsRelExRF;
            bottomThirdIndicesRelExRF = 1:round(nTrialsRelExRF * 1/3);

            % ct delay
%                 [h,p] = ttest2(rtInRFSortedByCTDelay(topThirdIndicesInRF), ...
%                         rtInRFSortedByCTDelay(bottomThirdIndicesInRF))
%                 [h,p] = ttest2(rtExRFSortedByCTDelay(topThirdIndicesExRF), ...
%                         rtExRFSortedByCTDelay(bottomThirdIndicesExRF))
% 
%                 % td delay
%                 [h,p] = ttest2(rtInRFSortedByTDDelay(topThirdIndicesInRF), ...
%                         rtInRFSortedByTDDelay(bottomThirdIndicesInRF))
%                 [h,p] = ttest2(rtExRFSortedByTDDelay(topThirdIndicesExRF), ...
%                         rtExRFSortedByTDDelay(bottomThirdIndicesExRF))
% 
%                 binEdges = 0.3:0.05:0.8;
% 
%                 % sorted by ct delay
%                 figure_tr_inch(10, 10);
%                 subaxis(2, 2, 1);
%                 hold on;
%                 histogram(rtInRFSortedByCTDelay(topThirdIndicesInRF), binEdges);
%                 histogram(rtInRFSortedByCTDelay(bottomThirdIndicesInRF), binEdges);
%                 title({'RT distribution of Top and Bottom Third', 'Extreme Cue-Target Delay Firing Rates InRF'});
%                 subaxis(2, 2, 3);
%                 hold on;
%                 histogram(rtExRFSortedByCTDelay(topThirdIndicesExRF), binEdges);
%                 histogram(rtExRFSortedByCTDelay(bottomThirdIndicesExRF), binEdges);
%                 title({'RT distribution of Top and Bottom Third', 'Extreme Cue-Target Delay Firing Rates ExRF'});
% 
%                 [mean(rtInRFSortedByCTDelay(topThirdIndicesInRF)) ...
%                         mean(rtInRFSortedByCTDelay(bottomThirdIndicesInRF)) ...
%                         mean(rtExRFSortedByCTDelay(topThirdIndicesExRF)) ...
%                         mean(rtExRFSortedByCTDelay(bottomThirdIndicesExRF))]
% 
%                 % sorted by td delay
%                 subaxis(2, 2, 2);
%                 hold on;
%                 histogram(rtInRFSortedByTDDelay(topThirdIndicesInRF), binEdges);
%                 histogram(rtInRFSortedByTDDelay(bottomThirdIndicesInRF), binEdges);
%                 title({'RT distribution of Top and Bottom Third', 'Extreme Target-Dim Delay Firing Rates InRF'});
%                 subaxis(2, 2, 4);
%                 hold on;
%                 histogram(rtExRFSortedByTDDelay(topThirdIndicesExRF), binEdges);
%                 histogram(rtExRFSortedByTDDelay(bottomThirdIndicesExRF), binEdges);
%                 title({'RT distribution of Top and Bottom Third', 'Extreme Target-Dim Delay Firing Rates ExRF'});
% 
%                 [mean(rtInRFSortedByTDDelay(topThirdIndicesInRF)) ...
%                         mean(rtInRFSortedByTDDelay(bottomThirdIndicesInRF)) ...
%                         mean(rtExRFSortedByTDDelay(topThirdIndicesExRF)) ...
%                         mean(rtExRFSortedByTDDelay(bottomThirdIndicesExRF))]
            % ideally shuffle test?
            % how to deal with correlated increases in firing in a trial?

            meanRTHoldInRFTopThirdFiringRateCTDelayAll(unitCount) = mean(rtHoldInRFSortedByCTDelay(topThirdIndicesHoldInRF));
            meanRTHoldInRFBottomThirdFiringRateCTDelayAll(unitCount) = mean(rtHoldInRFSortedByCTDelay(bottomThirdIndicesHoldInRF));
            meanRTRelInRFTopThirdFiringRateCTDelayAll(unitCount) = mean(rtRelInRFSortedByCTDelay(topThirdIndicesRelInRF));
            meanRTRelInRFBottomThirdFiringRateCTDelayAll(unitCount) = mean(rtRelInRFSortedByCTDelay(bottomThirdIndicesRelInRF));
            meanRTHoldInRFTopThirdFiringRateTDDelayAll(unitCount) = mean(rtHoldInRFSortedByTDDelay(topThirdIndicesHoldInRF));
            meanRTHoldInRFBottomThirdFiringRateTDDelayAll(unitCount) = mean(rtHoldInRFSortedByTDDelay(bottomThirdIndicesHoldInRF));
            meanRTHoldExRFTopThirdFiringRateCTDelayAll(unitCount) = mean(rtHoldExRFSortedByCTDelay(topThirdIndicesHoldExRF));
            meanRTHoldExRFBottomThirdFiringRateCTDelayAll(unitCount) = mean(rtHoldExRFSortedByCTDelay(bottomThirdIndicesHoldExRF));
            meanRTRelExRFTopThirdFiringRateCTDelayAll(unitCount) = mean(rtRelExRFSortedByCTDelay(topThirdIndicesRelExRF));
            meanRTRelExRFBottomThirdFiringRateCTDelayAll(unitCount) = mean(rtRelExRFSortedByCTDelay(bottomThirdIndicesRelExRF));
            meanRTHoldExRFTopThirdFiringRateTDDelayAll(unitCount) = mean(rtHoldExRFSortedByTDDelay(topThirdIndicesHoldExRF));
            meanRTHoldExRFBottomThirdFiringRateTDDelayAll(unitCount) = mean(rtHoldExRFSortedByTDDelay(bottomThirdIndicesHoldExRF));

            [corrCoefHoldInRFCTDelayRTAll(unitCount,1),corrCoefHoldInRFCTDelayRTAll(unitCount,2)] = corr(cueTargetDelayLongHoldInRFRate, rtHoldInRF);
            [corrCoefRelInRFCTDelayRTAll(unitCount,1),corrCoefRelInRFCTDelayRTAll(unitCount,2)] = corr(cueTargetDelayLongRelInRFRate, rtRelInRF);
            [corrCoefHoldInRFTDDelayRTAll(unitCount,1),corrCoefHoldInRFTDDelayRTAll(unitCount,2)] = corr(targetDimDelayLongHoldInRFRate, rtHoldInRF);
            [corrCoefHoldExRFCTDelayRTAll(unitCount,1),corrCoefHoldExRFCTDelayRTAll(unitCount,2)] = corr(cueTargetDelayLongHoldExRFRate, rtHoldExRF);
            [corrCoefRelExRFCTDelayRTAll(unitCount,1),corrCoefRelExRFCTDelayRTAll(unitCount,2)] = corr(cueTargetDelayLongRelExRFRate, rtRelExRF);
            [corrCoefHoldExRFTDDelayRTAll(unitCount,1),corrCoefHoldExRFTDDelayRTAll(unitCount,2)] = corr(targetDimDelayLongHoldExRFRate, rtHoldExRF);

            % need to use fisher xform (atanh) to test population
            % correlation coefficient
            corrCoefHoldInRFCTDelayRTAll(unitCount,3) = atanh(corrCoefHoldInRFCTDelayRTAll(unitCount,1));
            corrCoefRelInRFCTDelayRTAll(unitCount,3) = atanh(corrCoefRelInRFCTDelayRTAll(unitCount,1));
            corrCoefHoldInRFTDDelayRTAll(unitCount,3) = atanh(corrCoefHoldInRFTDDelayRTAll(unitCount,1));
            corrCoefHoldExRFCTDelayRTAll(unitCount,3) = atanh(corrCoefHoldExRFCTDelayRTAll(unitCount,1));
            corrCoefRelExRFCTDelayRTAll(unitCount,3) = atanh(corrCoefRelExRFCTDelayRTAll(unitCount,1));
            corrCoefHoldExRFTDDelayRTAll(unitCount,3) = atanh(corrCoefHoldExRFTDDelayRTAll(unitCount,1));

            
            
            
            
            inRFPreCueBaseline = ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(inRFLoc);
            exRFPreCueBaseline = ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(exRFLoc);

            % baseline correct the norm factor, account for cases where
            % suppression > enhancement
            % normalize including saccadic responses
            inRFNormFactor = max([ES.maxFiringRateBySpdfInclMotor - inRFPreCueBaseline; inRFPreCueBaseline - ES.minFiringRateBySpdfInclMotor]);
            exRFNormFactor = max([ES.maxFiringRateBySpdfInclMotor - exRFPreCueBaseline; exRFPreCueBaseline - ES.minFiringRateBySpdfInclMotor]);

            enterFixationSpdfInRFNorm(unitCount,:) = (ES.enterFixation.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            enterFixationSpdfExRFNorm(unitCount,:) = (ES.enterFixation.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            enterFixationSpdfInRFNormErr(unitCount,:) = ES.enterFixation.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            enterFixationSpdfExRFNormErr(unitCount,:) = ES.enterFixation.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            cueOnsetSpdfInRFNorm(unitCount,:) = (ES.cueOnset.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            cueOnsetSpdfExRFNorm(unitCount,:) = (ES.cueOnset.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            cueOnsetSpdfInRFNormErr(unitCount,:) = ES.cueOnset.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            cueOnsetSpdfExRFNormErr(unitCount,:) = ES.cueOnset.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            arrayOnsetHoldSpdfInRFNorm(unitCount,:) = (ES.arrayOnsetHold.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            arrayOnsetHoldSpdfExRFNorm(unitCount,:) = (ES.arrayOnsetHold.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            arrayOnsetHoldSpdfInRFNormErr(unitCount,:) = ES.arrayOnsetHold.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            arrayOnsetHoldSpdfExRFNormErr(unitCount,:) = ES.arrayOnsetHold.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            targetDimSpdfInRFNorm(unitCount,:) = (ES.targetDim.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            targetDimSpdfExRFNorm(unitCount,:) = (ES.targetDim.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            targetDimSpdfInRFNormErr(unitCount,:) = ES.targetDim.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            targetDimSpdfExRFNormErr(unitCount,:) = ES.targetDim.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            exitFixationSpdfInRFNorm(unitCount,:) = (ES.exitFixation.spdfByLoc(inRFLoc,:) - inRFPreCueBaseline) / inRFNormFactor;
            exitFixationSpdfExRFNorm(unitCount,:) = (ES.exitFixation.spdfByLoc(exRFLoc,:) - exRFPreCueBaseline) / exRFNormFactor;

            exitFixationSpdfInRFNormErr(unitCount,:) = ES.exitFixation.spdfErrByLoc(inRFLoc,:) / inRFNormFactor;
            exitFixationSpdfExRFNormErr(unitCount,:) = ES.exitFixation.spdfErrByLoc(exRFLoc,:) / exRFNormFactor;

            meanNormSpdfInRFAllWindowsAll(unitCount,:) = ([...
                    ES.averageFiringRatesBySpdf.preEnterFixation.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.postEnterFixation.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.preCueBaseline.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueResponse.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.arrayHoldResponse.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.targetDimResponse.byLoc(inRFLoc) ...
                    ES.averageFiringRatesBySpdf.preExitFixation.byLoc(inRFLoc)] - inRFPreCueBaseline) / inRFNormFactor;

            meanNormSpdfExRFAllWindowsAll(unitCount,:) = ([...
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
        end
    else
        fprintf('Skipping %s due to min firing rate requirement...\n', unitName);
    end
end

spdfInfo = var2struct(...
        enterFixationSpdfInRFNorm, ...
        enterFixationSpdfExRFNorm, ...
        enterFixationSpdfInRFNormErr, ...
        enterFixationSpdfExRFNormErr, ...
        cueOnsetSpdfInRFNorm, ...
        cueOnsetSpdfExRFNorm, ...
        cueOnsetSpdfInRFNormErr, ...
        cueOnsetSpdfExRFNormErr, ...
        arrayOnsetHoldSpdfInRFNorm, ...
        arrayOnsetHoldSpdfExRFNorm, ...
        arrayOnsetHoldSpdfInRFNormErr, ...
        arrayOnsetHoldSpdfExRFNormErr, ...
        targetDimSpdfInRFNorm, ...
        targetDimSpdfExRFNorm, ...
        targetDimSpdfInRFNormErr, ...
        targetDimSpdfExRFNormErr, ...
        exitFixationSpdfInRFNorm, ...
        exitFixationSpdfExRFNorm, ...
        exitFixationSpdfInRFNormErr, ...
        exitFixationSpdfExRFNormErr, ...
        meanNormSpdfInRFAllWindowsAll, ...
        meanNormSpdfExRFAllWindowsAll);


saveFileName = sprintf('%s/%s-muaAnalysisSummaryData-v%d.mat', outputDir, sessionName, v);
fprintf('\n');
fprintf('Writing file %s ...\n', saveFileName);
save(saveFileName, ...
        'unitNamesAll', ...
        'isSignificantStats', ...
        'infoRates', ...
        'attnIndices', ...
        'localization', ...
        'isInVPulvinar', ...
        'spdfInfo', ...
        'enterFixationT', ...
        'cueOnsetT', ...
        'arrayOnsetT', ...
        'targetDimT', ...
        'exitFixationT');