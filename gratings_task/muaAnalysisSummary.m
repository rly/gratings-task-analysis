
clear;
% processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
processedDataRootDir = 'Y:/rly/gratings-task-analysis/processed_data/';
dataDirRoot = 'Z:/ryanly/McCartney/originals';
muaDataDirRoot = 'Y:/rly/simple-mua-detection/processed_data/';
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';

v = 9;

nUnitsApprox = 900; % make sure this is an underestimate
unitNamesAll = cell(nUnitsApprox, 1);
isCell = false(nUnitsApprox, 1);
isBroadSpiking = false(nUnitsApprox, 1);
isNarrowSpiking = false(nUnitsApprox, 1);
isSignificantStats = false(nUnitsApprox, 10); % 5 periods > baseline, 5 periods info rate
infoRates = nan(nUnitsApprox, 5); % 5 periods
attnIndices = nan(nUnitsApprox, 2); % 2 delay periods
localization = cell(nUnitsApprox, 1);
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
statAlpha = 0.05/2; % reduce from 0.05 to account for multiple comparisons... but by how much?
% should also be running a lot of shuffle tests given the number of trials

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');

sessionIndAll = 1:37;
for k = 1:numel(sessionIndAll)
    sessionInd = sessionIndAll(k);

    %% load recording information
    recordingInfo = readRecordingInfo(recordingInfoFileName);
    struct2var(recordingInfo(sessionInd));
    pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

    %% load recording data
    fprintf('\n-------------------------------------------------------\n');
    fprintf('Gratings Task Analysis\n');
    fprintf('Loading %s...\n', pl2FilePath);
    fprintf('Session index: %d\n', sessionInd);

    tic;
    isLoadSpikes = 0;
    isLoadMua = 1;
    isLoadLfp = 0;
    isLoadSpkc = 0;
    isLoadDirect = 0;
    D = loadPL2(pl2FilePath, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
             spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

    processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
    if exist(processedDataDir, 'dir') == 0
        mkdir(processedDataDir);
    end
    fprintf('... done (%0.2f s).\n', toc);
    
    assert(numel(blockNames) == numel(D.blockStartTimes));
    blockName = strjoin(blockNames(gratingsTask3DIndices), '-');
    fprintf('Analyzing block names: %s.\n', blockName);    
    
    %% remove spike and event times not during task to save memory
    D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);
    
    %% summarize across population of cells    
    totalTimeOverall = sum(D.blockStopTimes(gratingsTask3DIndices) - D.blockStartTimes(gratingsTask3DIndices));
    minFiringRateOverall = 0.2;

    nUnits = numel(D.allMUAStructs);
    fprintf('-------------------------------------------------------------\n');
    fprintf('Processing %d units...\n', nUnits);

    for i = 1:nUnits
        spikeStruct = D.allMUAStructs{i};
        unitName = spikeStruct.name;
        spikeTimes = spikeStruct.ts;
        firingRateOverall = numel(spikeTimes) / totalTimeOverall;
%         fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
%                 nUnits, round(i/nUnits*100));

        if firingRateOverall >= minFiringRateOverall
            saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                    processedDataDir, unitName, blockName, v);

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
                
                unitNamesAll{unitCount} = unitName;
                
                if ismember(spikeStruct.channelID, pldChannels)
                    localization{unitCount} = 'PLd';
                elseif ismember(spikeStruct.channelID, plvChannels)
                    localization{unitCount} = 'PLv';
                elseif ismember(spikeStruct.channelID, pmChannels)
                    localization{unitCount} = 'PM';
                elseif ismember(spikeStruct.channelID, piChannels)
                    localization{unitCount} = 'PI';
                else
                    localization{unitCount} = '';
                end
        
                isCell(unitCount) = strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking');
                if isCell(unitCount)
                    isBroadSpiking(unitCount) = strcmp(spikeStruct.physClass, 'Broad-Spiking');
                    isNarrowSpiking(unitCount) = strcmp(spikeStruct.physClass, 'Narrow-Spiking');
                else
                    isBroadSpiking(unitCount) = false;
                    isNarrowSpiking(unitCount) = false;
                end

                isSignificantStats(unitCount,:) = [...
                        min([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) ...
                        min([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) ...
                        min([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) ...
                        min([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) ...
                        min([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) ...
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
            end
        end
    end
end

%% summarize
nUnitsAll = unitCount;
nCells = sum(isBroadSpiking | isNarrowSpiking);
assert(sum(isCell) == nCells);

isSignificantAnyTaskMod = isCell & any(isSignificantStats(:,1:5), 2);
isSignificantCueActivity = isCell & isSignificantStats(:,1); % note that this uses the stricter statAlpha
isSignificantAnySpatialSelectivity = isCell & any(isSignificantStats(:,6:10), 2);
isSignificantEvokedSelectivity = isCell & any(isSignificantStats(:,[6 8 10]), 2);
isSignificantDelaySelectivity = isCell & any(isSignificantStats(:,[7 9]), 2);
isSigSelectEvokedNotDelay = isCell & any(isSignificantStats(:,[6 8 10]), 2) & any(isSignificantStats(:,[7 9]), 2);
isSigSelectOnlyEvoked = isCell & any(isSignificantStats(:,[6 8 10]), 2) & ~any(isSignificantStats(:,[7 9]), 2);
isSigSelectOnlyDelay = isCell & ~any(isSignificantStats(:,[6 8 10]), 2) & any(isSignificantStats(:,[7 9]), 2);
isInPulvinar = strcmp(localization, 'PLd') | strcmp(localization, 'PLv') | strcmp(localization, 'PM') | strcmp(localization, 'PI');

fprintf('-------------------------------\n');
fprintf('%d/%d = %d%% cells were classified as broad-spiking and %d/%d = %d%% as narrow-spiking.\n', ...
        sum(isBroadSpiking & isCell), nCells, ...
        round(sum(isBroadSpiking & isCell)/nCells * 100), ...
        sum(isNarrowSpiking & isCell), nCells, ...
        round(sum(isNarrowSpiking & isCell)/nCells * 100));
fprintf('%d/%d = %d%% cells were localized to the pulvinar.\n', ...
        sum(isCell & isInPulvinar), nCells, ...
        round(sum(isCell & isInPulvinar)/nCells * 100));
fprintf('\n');

fprintf('%d/%d = %d%% cells show significant task modulation compared to baseline.\n', ...
        sum(isSignificantAnyTaskMod), nCells, ...
        round(sum(isSignificantAnyTaskMod)/nCells * 100));
fprintf('%d/%d = %d%% cells show significant spatial selectivity during some task period.\n', ...
        sum(isSignificantAnySpatialSelectivity), nCells, ...
        round(sum(isSignificantAnySpatialSelectivity)/nCells * 100));
fprintf('%d/%d = %d%% cells show significant spatial selectivity in a task period involving a visual change.\n', ...
        sum(isSignificantEvokedSelectivity), nCells, ...
        round(sum(isSignificantEvokedSelectivity)/nCells * 100));
fprintf('%d/%d = %d%% cells show significant spatial selectivity in a task period involving sustained attention.\n', ...
        sum(isSignificantDelaySelectivity), nCells, ...
        round(sum(isSignificantDelaySelectivity)/nCells * 100));
fprintf('\n');

fprintf('%d/%d = %d%% pulvinar cells show significant task modulation compared to baseline.\n', ...
        sum(isSignificantAnyTaskMod & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantAnyTaskMod & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar cells show significant spatial selectivity during some task period.\n', ...
        sum(isSignificantAnySpatialSelectivity & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantAnySpatialSelectivity & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar cells show spatial selectivity in a task period involving a visual change.\n', ...
        sum(isSignificantEvokedSelectivity & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantEvokedSelectivity & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('%d/%d = %d%% pulvinar cells show spatial selectivity in a task period involving sustained attention.\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), sum(isCell & isInPulvinar), ...
        round(sum(isSignificantDelaySelectivity & isInPulvinar)/sum(isCell & isInPulvinar) * 100));
fprintf('\n');

fprintf('Of the %d cells that show significant spatial selectivity during some task period,\n', ...
        sum(isSignificantAnySpatialSelectivity));
fprintf('\t%d (%d%%) show significant spatial selectivity after a visual change AND during attention,\n', ...
        sum(isSigSelectEvokedNotDelay), ...
        round(sum(isSigSelectEvokedNotDelay)/sum(isSignificantAnySpatialSelectivity) * 100));
fprintf('\t%d (%d%%) show significant spatial selectivity ONLY after a visual change, and\n', ...
        sum(isSigSelectOnlyEvoked), ...
        round(sum(isSigSelectOnlyEvoked)/sum(isSignificantAnySpatialSelectivity) * 100));
fprintf('\t%d (%d%%) show significant spatial selectivity ONLY during attention.\n\n', ...
        sum(isSigSelectOnlyDelay), ...
        round(sum(isSigSelectOnlyDelay)/sum(isSignificantAnySpatialSelectivity) * 100));
fprintf('Of the %d cells that show significant task modulation compared to baseline, \n\t%d (%d%%) are broad-spiking and %d (%d%%) are narrow-spiking.\n', ...
        sum(isSignificantAnyTaskMod), ...
        sum(isBroadSpiking(isSignificantAnyTaskMod)), ...
        round(sum(isBroadSpiking(isSignificantAnyTaskMod))/sum(isSignificantAnyTaskMod) * 100), ...
        sum(isNarrowSpiking(isSignificantAnyTaskMod)), ...
        round(sum(isNarrowSpiking(isSignificantAnyTaskMod))/sum(isSignificantAnyTaskMod) * 100));
fprintf('Of the %d cells that show significant spatial selectivity after a visual change, \n\t%d (%d%%) are broad-spiking and %d (%d%%) are narrow-spiking.\n', ...
        sum(isSignificantEvokedSelectivity), ...
        sum(isBroadSpiking(isSignificantEvokedSelectivity)), ...
        round(sum(isBroadSpiking(isSignificantEvokedSelectivity))/sum(isSignificantEvokedSelectivity) * 100), ...
        sum(isNarrowSpiking(isSignificantEvokedSelectivity)), ...
        round(sum(isNarrowSpiking(isSignificantEvokedSelectivity))/sum(isSignificantEvokedSelectivity) * 100));
fprintf('Of the %d cells that show significant spatial selectivity during attention, \n\t%d (%d%%) are broad-spiking and %d (%d%%) are narrow-spiking.\n', ...
        sum(isSignificantDelaySelectivity), ...
        sum(isBroadSpiking(isSignificantDelaySelectivity)), ...
        round(sum(isBroadSpiking(isSignificantDelaySelectivity))/sum(isSignificantDelaySelectivity) * 100), ...
        sum(isNarrowSpiking(isSignificantDelaySelectivity)), ...
        round(sum(isNarrowSpiking(isSignificantDelaySelectivity))/sum(isSignificantDelaySelectivity) * 100));
fprintf('\n');

fprintf('Of the %d cells in the pulvinar that show spatial selectivity during attention, \n\t%d (%d%%) are in PLd, %d (%d%%) are in PLv, %d (%d%%) are in PM, and %d (%d%%) are in PI\n', ...
        sum(isSignificantDelaySelectivity & isInPulvinar), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PM')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PM'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100), ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PI')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PI'))/sum(isSignificantDelaySelectivity & isInPulvinar) * 100));
fprintf('\n');

fprintf('\t%d/%d = %d%% PLd cells\n', ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd')), ...
        sum(isCell & strcmp(localization, 'PLd')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
fprintf('\t%d/%d = %d%% PLv cells\n', ...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv')), ...
        sum(isCell & strcmp(localization, 'PLv')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
fprintf('\t%d/%d = %d%% PM cells\n',...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PM')), ...
        sum(isCell & strcmp(localization, 'PM')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
fprintf('\t%d/%d = %d%% PI cells\n',...
        sum(isSignificantDelaySelectivity & strcmp(localization, 'PI')), ...
        sum(isCell & strcmp(localization, 'PI')), ...
        round(sum(isSignificantDelaySelectivity & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\n');

fprintf('Significant selectivity in cue-target delay:\n');
fprintf('\t%d/%d = %d%% PLd cells\n', ...
        sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PLd')), ...
        sum(isCell & strcmp(localization, 'PLd')), ...
        round(sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
fprintf('\t%d/%d = %d%% PLv cells\n', ...
        sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PLv')), ...
        sum(isCell & strcmp(localization, 'PLv')), ...
        round(sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
fprintf('\t%d/%d = %d%% PM cells\n',...
        sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PM')), ...
        sum(isCell & strcmp(localization, 'PM')), ...
        round(sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
fprintf('\t%d/%d = %d%% PI cells\n',...
        sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PI')), ...
        sum(isCell & strcmp(localization, 'PI')), ...
        round(sum(isCell & isSignificantStats(:,7) & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\n');

fprintf('Significant selectivity in target-dim delay:\n');
fprintf('\t%d/%d = %d%% PLd cells\n', ...
        sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PLd')), ...
        sum(isCell & strcmp(localization, 'PLd')), ...
        round(sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
fprintf('\t%d/%d = %d%% PLv cells\n', ...
        sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PLv')), ...
        sum(isCell & strcmp(localization, 'PLv')), ...
        round(sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
fprintf('\t%d/%d = %d%% PM cells\n',...
        sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PM')), ...
        sum(isCell & strcmp(localization, 'PM')), ...
        round(sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
fprintf('\t%d/%d = %d%% PI cells\n',...
        sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PI')), ...
        sum(isCell & strcmp(localization, 'PI')), ...
        round(sum(isCell & isSignificantStats(:,9) & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\n');

fprintf('Significant selectivity in both delay periods:\n');
fprintf('\t%d/%d = %d%% PLd cells\n', ...
        sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PLd')), ...
        sum(isCell & strcmp(localization, 'PLd')), ...
        round(sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PLd'))/sum(isCell & strcmp(localization, 'PLd')) * 100));
fprintf('\t%d/%d = %d%% PLv cells\n', ...
        sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PLv')), ...
        sum(isCell & strcmp(localization, 'PLv')), ...
        round(sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PLv'))/sum(isCell & strcmp(localization, 'PLv')) * 100));
fprintf('\t%d/%d = %d%% PM cells\n',...
        sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PM')), ...
        sum(isCell & strcmp(localization, 'PM')), ...
        round(sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PM'))/sum(isCell & strcmp(localization, 'PM')) * 100));
fprintf('\t%d/%d = %d%% PI cells\n',...
        sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PI')), ...
        sum(isCell & strcmp(localization, 'PI')), ...
        round(sum(isCell & isSignificantStats(:,7) & isSignificantStats(:,9) & strcmp(localization, 'PI'))/sum(isCell & strcmp(localization, 'PI')) * 100));
fprintf('\n');

stop
%% Pie Chart localizing attentional modulation in the pulvinar
% We found x pulvinar cells that show significant attentional modulation. 
% In which subdivision of the pulvinar are they located?
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.1, 'MB', 0.1, 'MT', 0.1, 'MR', 0.1);
set(gcf, 'Color', 'white');
hold on;

g2 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g2all = isCell & strcmp(localization, 'PLd');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PLv');
g3all = isCell & strcmp(localization, 'PLv');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g4all = isCell & strcmp(localization, 'PM');
g5 = isSignificantDelaySelectivity & strcmp(localization, 'PI');
g5all = isCell & strcmp(localization, 'PI');

x = [sum(g2)/sum(g2all) sum(g3)/sum(g3all) sum(g4)/sum(g4all) sum(g5)/sum(g5all)];
h = pie(x / sum(x));
h(1).FaceColor = [75 172 198]/255; % use colors from the ppt
h(3).FaceColor = [255 192 0]/255;
h(5).FaceColor = [255 0 0]/255;
h(7).FaceColor = [100 50 190]/255;
h(2).String = sprintf('PLd: %s', h(2).String);
h(4).String = sprintf('PLv: %s', h(4).String);
h(6).String = sprintf('PM: %s', h(6).String);
h(8).String = sprintf('PI: %s', h(8).String);
set(h(2), 'FontSize', 20);
set(h(4), 'FontSize', 20);
set(h(6), 'FontSize', 20);
set(h(8), 'FontSize', 20);
axis off;
% adjust the locations of the text by hand

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by Sig in Which 
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity;
g2 = isCell & isInPulvinar; % & any(isSignificantStats(:,7), 2) & ~any(isSignificantStats(:,9), 2);
% g3 = isCell & pulLocalization; % & ~any(isSignificantStats(:,7), 2) & any(isSignificantStats(:,9), 2);
% g4 = isCell & pulLocalization; % & any(isSignificantStats(:,7), 2) & any(isSignificantStats(:,9), 2);


plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
% plot(infoRates(g3, 2), infoRates(g3, 4), ...
%         '.', 'MarkerSize', 20, 'Color', [1 0.5 1]);
% plot(infoRates(g4, 2), infoRates(g4, 4), ...
%         '.', 'MarkerSize', 20, 'Color', [0.5 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = max([currXLim(2) currYLim(2)]);
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by PLd vs PM vs PI
figure_tr_inch(6.5, 6.5);
subaxis(1, 1, 1, 'ML', 0.22, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PLv');
g5 = isSignificantDelaySelectivity & strcmp(localization, 'PI');

cols = lines(5);
cols = [cols(1,:); cols(3,:); cols(5,:); cols(2,:)];

maxLim = 0.55;
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

% plot(infoRates(g1,2), infoRates(g1,4), ...
%         '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
% plot(infoRates(g2, 2), infoRates(g2, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(1,:));
% plot(infoRates(g3, 2), infoRates(g3, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(2,:));
% plot(infoRates(g4, 2), infoRates(g4, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(3,:));
% plot(infoRates(g5, 2), infoRates(g5, 4), ...
%         '.', 'MarkerSize', 20, 'Color', cols(4,:));
sh = scatter(infoRates(g1,2), infoRates(g1,4), 100, 0.85*ones(1, 3), 'MarkerFaceColor', 0.85*ones(3, 1));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g2,2), infoRates(g2,4), 100, cols(1,:), 'MarkerFaceColor', cols(1,:));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g3,2), infoRates(g3,4), 100, cols(2,:), 'MarkerFaceColor', cols(2,:));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g4,2), infoRates(g4,4), 100, cols(3,:), 'MarkerFaceColor', cols(3,:));
sh.MarkerFaceAlpha = 0.9;
sh = scatter(infoRates(g5,2), infoRates(g5,4), 100, cols(4,:), 'MarkerFaceColor', cols(4,:));
sh.MarkerFaceAlpha = 0.9;

currXLim = xlim();
currYLim = ylim();
%max([currXLim(2) currYLim(2)]);
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
% xlabel({'Selectivity in', 'Cue-Target Delay (bits/s)'});
% ylabel({'Selectivity in', 'Target-Dim Delay (bits/s)'});
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-PLdvsPMvsPI-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

print(sprintf('%s/posterFigs/infoRate1vsInfoRate2-byPulSubdivision-v%d.svg', processedDataRootDir, v), '-dsvg');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by PLd vs PLv vs PM vs PI
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PLv');
g5 = isSignificantDelaySelectivity & strcmp(localization, 'PI');

plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(infoRates(g3, 2), infoRates(g3, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
plot(infoRates(g4, 2), infoRates(g4, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
plot(infoRates(g5, 2), infoRates(g5, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0.7 0]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.2;
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-PLdvsPMvsPI-zoom-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by PLd vs PLv vs PM vs PI
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantEvokedSelectivity & isInPulvinar;
g2 = isSignificantEvokedSelectivity & strcmp(localization, 'PM');
g3 = isSignificantEvokedSelectivity & strcmp(localization, 'PLd');
g4 = isSignificantEvokedSelectivity & strcmp(localization, 'PLv');
g5 = isSignificantEvokedSelectivity & strcmp(localization, 'PI');

plot(meanNormSpdfInRFAllWindowsAll(g1, 4), meanNormSpdfInRFAllWindowsAll(g1, 6), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(meanNormSpdfInRFAllWindowsAll(g2, 4), meanNormSpdfInRFAllWindowsAll(g2, 6), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(meanNormSpdfInRFAllWindowsAll(g3, 4), meanNormSpdfInRFAllWindowsAll(g3, 6), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
plot(meanNormSpdfInRFAllWindowsAll(g4, 4), meanNormSpdfInRFAllWindowsAll(g4, 6), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
plot(meanNormSpdfInRFAllWindowsAll(g5, 4), meanNormSpdfInRFAllWindowsAll(g5, 6), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0.7 0]);
currXLim = xlim();
currYLim = ylim();
maxLim = max([currXLim currYLim]);
plot([-maxLim maxLim], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([-maxLim maxLim]);
ylim([-maxLim maxLim]);
axis square;
box off;
xlabel({'Normalized Firing Rate After Cue Onset'});
ylabel({'Normalized Firing Rate After Array Onset'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-cueResponseArrayResponse-PLdvsPMvsPI-zoom-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');


%% Cue-Target Delay Info vs Target-Dim Delay AI, Color by PLd vs PLv vs PM vs PI
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
% TODO for now use significance based on info rate BUT should use AI
% shuffle
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & strcmp(localization, 'PLd');
g3 = isSignificantDelaySelectivity & strcmp(localization, 'PM');
g4 = isSignificantDelaySelectivity & strcmp(localization, 'PI');

plot(attnIndices(g1, 1), attnIndices(g1, 2), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(attnIndices(g2, 1), attnIndices(g2, 2), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(attnIndices(g3, 1), attnIndices(g3, 2), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
plot(attnIndices(g4, 1), attnIndices(g4, 2), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.7; % set manually
assert(~any(any(abs(attnIndices(g1 | g2 | g3 | g4,:) > maxLim))));
plot([0 0], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([-maxLim maxLim]);
ylim([-maxLim maxLim]);
axis square;
box off;
xlabel({'Attention Index during', 'Cue-Target Delay'});
ylabel({'Attention Index during', 'Target-Dim Delay'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-AI1vsAI2-PLdvsPMvsPI-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by BS vs NS
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & isBroadSpiking & isInPulvinar;
g3 = isSignificantDelaySelectivity & isNarrowSpiking & isInPulvinar;

plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(infoRates(g3, 2), infoRates(g3, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = max([currXLim(2) currYLim(2)]);
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-BSvsNS-%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by BS vs NS - ZOOM
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & isBroadSpiking & isInPulvinar;
g3 = isSignificantDelaySelectivity & isNarrowSpiking & isInPulvinar;

plot(infoRates(g1, 2), infoRates(g1, 4), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(infoRates(g2, 2), infoRates(g2, 4), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(infoRates(g3, 2), infoRates(g3, 4), ...
        '.', 'MarkerSize', 20, 'Color', [0.2 0 1]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.2;
plot([0 maxLim], [0 maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([0 maxLim]);
ylim([0 maxLim]);
axis square;
box off;
xlabel({'Information Rate during', 'Cue-Target Delay (bits/sec)'});
ylabel({'Information Rate during', 'Target-Dim Delay (bits/sec)'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-infoRate1vsInfoRate2-BSvsNS-zoom-%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

%% Cue-Target Delay Info vs Target-Dim Delay Info, Color by BS vs NS
figure_tr_inch(4, 4);
subaxis(1, 1, 1, 'ML', 0.25, 'MB', 0.22, 'MT', 0.05);
set(gcf, 'Color', 'white');
hold on;
g1 = isCell & ~isSignificantDelaySelectivity & isInPulvinar;
g2 = isSignificantDelaySelectivity & isBroadSpiking & isInPulvinar;
g3 = isSignificantDelaySelectivity & isNarrowSpiking & isInPulvinar;

plot(attnIndices(g1, 1), attnIndices(g1, 2), ...
        '.', 'MarkerSize', 20, 'Color', 0.8*ones(3, 1));
plot(attnIndices(g2, 1), attnIndices(g2, 2), ...
        '.', 'MarkerSize', 20, 'Color', [1 0.5 0]);
plot(attnIndices(g3, 1), attnIndices(g3, 2), ...
        '.', 'MarkerSize', 20, 'Color', [0.75 0.25 0.75]);
currXLim = xlim();
currYLim = ylim();
maxLim = 0.7; % set manually
assert(~any(any(abs(attnIndices(g1 | g2 | g3 | g4,:) > maxLim))));
plot([0 0], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [0 0], '-', 'Color', 0.3 * ones(3, 1));
plot([-maxLim maxLim], [-maxLim maxLim], '-', 'Color', 0.3 * ones(3, 1));
xlim([-maxLim maxLim]);
ylim([-maxLim maxLim]);
axis square;
box off;
xlabel({'Attention Index during', 'Cue-Target Delay'});
ylabel({'Attention Index during', 'Target-Dim Delay'});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/allSessions-AI1vsAI2-BSvsNS-%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');



%% test diff RT for top third vs bottom third firing rates in delay periods
meanRTRelInRFDiffThirdFiringRateCTDelayAll = meanRTRelInRFTopThirdFiringRateCTDelayAll - meanRTRelInRFBottomThirdFiringRateCTDelayAll;
meanRTRelExRFDiffThirdFiringRateCTDelayAll = meanRTRelExRFTopThirdFiringRateCTDelayAll - meanRTRelExRFBottomThirdFiringRateCTDelayAll;
meanRTHoldInRFDiffThirdFiringRateCTDelayAll = meanRTHoldInRFTopThirdFiringRateCTDelayAll - meanRTHoldInRFBottomThirdFiringRateCTDelayAll;
meanRTHoldExRFDiffThirdFiringRateCTDelayAll = meanRTHoldExRFTopThirdFiringRateCTDelayAll - meanRTHoldExRFBottomThirdFiringRateCTDelayAll;
meanRTHoldInRFDiffThirdFiringRateTDDelayAll = meanRTHoldInRFTopThirdFiringRateTDDelayAll - meanRTHoldInRFBottomThirdFiringRateTDDelayAll;
meanRTHoldExRFDiffThirdFiringRateTDDelayAll = meanRTHoldExRFTopThirdFiringRateTDDelayAll - meanRTHoldExRFBottomThirdFiringRateTDDelayAll;

meanRTRelInRFDiffThirdFiringRateCTDelayPul = meanRTRelInRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTRelExRFDiffThirdFiringRateCTDelayPul = meanRTRelExRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTHoldInRFDiffThirdFiringRateCTDelayPul = meanRTHoldInRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTHoldExRFDiffThirdFiringRateCTDelayPul = meanRTHoldExRFDiffThirdFiringRateCTDelayAll(isCell & isInPulvinar);
meanRTHoldInRFDiffThirdFiringRateTDDelayPul = meanRTHoldInRFDiffThirdFiringRateTDDelayAll(isCell & isInPulvinar);
meanRTHoldExRFDiffThirdFiringRateTDDelayPul = meanRTHoldExRFDiffThirdFiringRateTDDelayAll(isCell & isInPulvinar);

[~,p] = ttest(meanRTRelInRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTRelExRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTHoldInRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTHoldExRFDiffThirdFiringRateCTDelayPul)
[~,p] = ttest(meanRTHoldInRFDiffThirdFiringRateTDDelayPul)
[~,p] = ttest(meanRTHoldExRFDiffThirdFiringRateTDDelayPul)


cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);
binEdges = -0.1:0.01:0.1;

figure_tr_inch(9, 6);
set(gcf, 'Color', 'w');
plotHs = nan(6, 1);
plotHs(1) = subaxis(2, 3, 1);
hold on;
hist1 = histogram(meanRTRelInRFDiffThirdFiringRateCTDelayPul, binEdges);
hist1.FaceColor = inRFCol;
plotHs(2) = subaxis(2, 3, 2);
hold on;
hist2 = histogram(meanRTHoldInRFDiffThirdFiringRateCTDelayPul, binEdges);
hist2.FaceColor = inRFCol;
plotHs(3) = subaxis(2, 3, 3);
hold on;
hist3 = histogram(meanRTHoldInRFDiffThirdFiringRateTDDelayPul, binEdges);
hist3.FaceColor = inRFCol;
plotHs(4) = subaxis(2, 3, 4);
hold on;
hist4 = histogram(meanRTRelExRFDiffThirdFiringRateCTDelayPul, binEdges);
hist4.FaceColor = exRFCol;
plotHs(5) = subaxis(2, 3, 5);
hold on;
hist5 = histogram(meanRTHoldExRFDiffThirdFiringRateCTDelayPul, binEdges);
hist5.FaceColor = exRFCol;
plotHs(6) = subaxis(2, 3, 6);
hold on;
hist6 = histogram(meanRTHoldExRFDiffThirdFiringRateTDDelayPul, binEdges);
hist6.FaceColor = exRFCol;

% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);

% note: interestingly, sessions 1-4 showed more significant
% difference shifted right than sessions 5-6,8-9 which showed more
% significant difference shifted left on hold trials
% both groups showed more sig diff shifted right on release trials

% TODO create model, compare how much variability in RT is explained by CT
% delay and TD delay, using partial determination analysis, which removes
% the correlation between the two. 

%% 
% median/third split may be more appropriate than a regression since 
% firing rates and the correlation of firing rate and RT may be more
% stepwise? there may also be a natural split, e.g. zero vs nonzero firing
% average correlation coefficient is not the way to do this. better to have
% different sessions in the model:
% https://stats.stackexchange.com/questions/8019/averaging-correlation-values
% or transform the r values using Fisher transform

corrCoefHoldInRFCTDelayRTPul = corrCoefHoldInRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefRelInRFCTDelayRTPul = corrCoefRelInRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefHoldInRFTDDelayRTPul = corrCoefHoldInRFTDDelayRTAll(isCell & isInPulvinar,3);
corrCoefHoldExRFCTDelayRTPul = corrCoefHoldExRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefRelExRFCTDelayRTPul = corrCoefRelExRFCTDelayRTAll(isCell & isInPulvinar,3);
corrCoefHoldExRFTDDelayRTPul = corrCoefHoldExRFTDDelayRTAll(isCell & isInPulvinar,3);

[~,p] = ttest(corrCoefRelInRFCTDelayRTPul)
[~,p] = ttest(corrCoefRelExRFCTDelayRTPul)
[~,p] = ttest(corrCoefHoldInRFCTDelayRTPul)
[~,p] = ttest(corrCoefHoldExRFCTDelayRTPul)
[~,p] = ttest(corrCoefHoldInRFTDDelayRTPul)
[~,p] = ttest(corrCoefHoldExRFTDDelayRTPul)


cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);
binEdges = -0.3:0.05:0.3;

figure_tr_inch(9, 6);
set(gcf, 'Color', 'w');
plotHs = nan(6, 1);
plotHs(1) = subaxis(2, 3, 1);
hold on;
hist1 = histogram(corrCoefRelInRFCTDelayRTPul, binEdges);
hist1.FaceColor = inRFCol;
plotHs(2) = subaxis(2, 3, 2);
hold on;
hist2 = histogram(corrCoefHoldInRFCTDelayRTPul, binEdges);
hist2.FaceColor = inRFCol;
plotHs(3) = subaxis(2, 3, 3);
hold on;
hist3 = histogram(corrCoefHoldInRFTDDelayRTPul, binEdges);
hist3.FaceColor = inRFCol;
plotHs(4) = subaxis(2, 3, 4);
hold on;
hist4 = histogram(corrCoefRelExRFCTDelayRTPul, binEdges);
hist4.FaceColor = exRFCol;
plotHs(5) = subaxis(2, 3, 5);
hold on;
hist5 = histogram(corrCoefHoldExRFCTDelayRTPul, binEdges);
hist5.FaceColor = exRFCol;
plotHs(6) = subaxis(2, 3, 6);
hold on;
hist6 = histogram(corrCoefHoldExRFTDDelayRTPul, binEdges);
hist6.FaceColor = exRFCol;

% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);

%% per-condition baseline-corrected normalized mean
fprintf('\n');
subdivisions = {'PM', 'PLd', 'PLv', 'PI', 'PUL', 'all'};
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(nUnitsAll, 1);
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = isInPulvinar;
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    condition = isCell & isInSubdivision;% & all(targetDimSpdfInRFNorm > -1,2);
    
    enterFixationSpdfInRFNormSub = trimNanRows(enterFixationSpdfInRFNorm(condition,:));
    enterFixationSpdfExRFNormSub = trimNanRows(enterFixationSpdfExRFNorm(condition,:));
    cueOnsetSpdfInRFNormSub = trimNanRows(cueOnsetSpdfInRFNorm(condition,:));
    cueOnsetSpdfExRFNormSub = trimNanRows(cueOnsetSpdfExRFNorm(condition,:));
    arrayOnsetHoldSpdfInRFNormSub = trimNanRows(arrayOnsetHoldSpdfInRFNorm(condition,:));
    arrayOnsetHoldSpdfExRFNormSub = trimNanRows(arrayOnsetHoldSpdfExRFNorm(condition,:));
    targetDimSpdfInRFNormSub = trimNanRows(targetDimSpdfInRFNorm(condition,:));
    targetDimSpdfExRFNormSub = trimNanRows(targetDimSpdfExRFNorm(condition,:));
    exitFixationSpdfInRFNormSub = trimNanRows(exitFixationSpdfInRFNorm(condition,:));
    exitFixationSpdfExRFNormSub = trimNanRows(exitFixationSpdfExRFNorm(condition,:));

    fprintf('\t%s: %d cells\n', subdivision, sum(condition));

    enterFixationT = ES.enterFixationT - ES.enterFixationWindow(1);
    cueOnsetT = ES.cueOnsetT - ES.cueOnsetWindow(1);
    arrayOnsetT = ES.arrayOnsetT - ES.arrayOnsetWindow(1);
    targetDimT = ES.targetDimT - ES.targetDimWindow(1);
    exitFixationT = ES.exitFixationT - ES.exitFixationWindow(1);

    plotFileName = sprintf('%s/allSessions-%s-meanSpdfs3-v%d.png', processedDataRootDir, subdivision, v);

    quickSpdfAllEvents3InARowPopMean(cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
            arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
            targetDimSpdfExRFNormSub, cueOnsetT, arrayOnsetT, targetDimT, plotFileName);

    plotFileName = sprintf('%s/allSessions-%s-meanSpdfs5-v%d.png', processedDataRootDir, subdivision, v);
    
    quickSpdfAllEvents5InARowPopMean(enterFixationSpdfInRFNormSub, enterFixationSpdfExRFNormSub, ...
            cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
            arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
            targetDimSpdfExRFNormSub, exitFixationSpdfInRFNormSub, exitFixationSpdfExRFNormSub, ...
            enterFixationT, cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, plotFileName);
    
    print(sprintf('%s/allSessions-%s-meanSpdfs5-v%d.svg', processedDataRootDir, subdivision, v), '-dsvg');
end

%% mega figure of tiny bc normalized plots per unit by subdivision
fprintf('\n');
subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(nUnitsAll, 1);
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = isInPulvinar;
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    condition = isCell & isInSubdivision;
    unitNamesSub = unitNamesAll(condition);
    
    enterFixationSpdfInRFNormSub = (enterFixationSpdfInRFNorm(condition,:));
    enterFixationSpdfExRFNormSub = (enterFixationSpdfExRFNorm(condition,:));
    cueOnsetSpdfInRFNormSub = (cueOnsetSpdfInRFNorm(condition,:));
    cueOnsetSpdfExRFNormSub = (cueOnsetSpdfExRFNorm(condition,:));
    arrayOnsetHoldSpdfInRFNormSub = (arrayOnsetHoldSpdfInRFNorm(condition,:));
    arrayOnsetHoldSpdfExRFNormSub = (arrayOnsetHoldSpdfExRFNorm(condition,:));
    targetDimSpdfInRFNormSub = (targetDimSpdfInRFNorm(condition,:));
    targetDimSpdfExRFNormSub = (targetDimSpdfExRFNorm(condition,:));
    exitFixationSpdfInRFNormSub = (exitFixationSpdfInRFNorm(condition,:));
    exitFixationSpdfExRFNormSub = (exitFixationSpdfExRFNorm(condition,:));
    
    enterFixationSpdfInRFNormErrSub = (enterFixationSpdfInRFNormErr(condition,:));
    enterFixationSpdfExRFNormErrSub = (enterFixationSpdfExRFNormErr(condition,:));
    cueOnsetSpdfInRFNormErrSub = (cueOnsetSpdfInRFNormErr(condition,:));
    cueOnsetSpdfExRFNormErrSub = (cueOnsetSpdfExRFNormErr(condition,:));
    arrayOnsetHoldSpdfInRFNormErrSub = (arrayOnsetHoldSpdfInRFNormErr(condition,:));
    arrayOnsetHoldSpdfExRFNormErrSub = (arrayOnsetHoldSpdfExRFNormErr(condition,:));
    targetDimSpdfInRFNormErrSub = (targetDimSpdfInRFNormErr(condition,:));
    targetDimSpdfExRFNormErrSub = (targetDimSpdfExRFNormErr(condition,:));
    exitFixationSpdfInRFNormErrSub = (exitFixationSpdfInRFNormErr(condition,:));
    exitFixationSpdfExRFNormErrSub = (exitFixationSpdfExRFNormErr(condition,:));

    fprintf('\t%s: %d cells\n', subdivision, sum(condition));
    
    enterFixationT = ES.enterFixationT - ES.enterFixationWindow(1);
    cueOnsetT = ES.cueOnsetT - ES.cueOnsetWindow(1);
    arrayOnsetT = ES.arrayOnsetT - ES.arrayOnsetWindow(1);
    targetDimT = ES.targetDimT - ES.targetDimWindow(1);
    exitFixationT = ES.exitFixationT - ES.exitFixationWindow(1);
    
    titleBase = sprintf('%s Cells: Enter Fixation', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf1-enterFixation-v%d', processedDataRootDir, subdivision, v);
    makeTinyPlotsOfPopulation(enterFixationSpdfInRFNormSub, enterFixationSpdfInRFNormErrSub, ...
            enterFixationSpdfExRFNormSub, enterFixationSpdfExRFNormErrSub, enterFixationT, unitNamesSub, titleBase, plotFileBaseName);
    
    titleBase = sprintf('%s Cells: Cue Onset', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf2-cueOnset-v%d', processedDataRootDir, subdivision, v);
    makeTinyPlotsOfPopulation(cueOnsetSpdfInRFNormSub, cueOnsetSpdfInRFNormErrSub, ...
            cueOnsetSpdfExRFNormSub, cueOnsetSpdfExRFNormErrSub, cueOnsetT, unitNamesSub, titleBase, plotFileBaseName);
    
    titleBase = sprintf('%s Cells: Array Onset Hold', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf3-arrayOnsetHold-v%d', processedDataRootDir, subdivision, v);
    makeTinyPlotsOfPopulation(arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfInRFNormErrSub, ...
            arrayOnsetHoldSpdfExRFNormSub, arrayOnsetHoldSpdfExRFNormErrSub, arrayOnsetT, unitNamesSub, titleBase, plotFileBaseName);
    
    titleBase = sprintf('%s Cells: Target Dim', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf4-targetDim-v%d', processedDataRootDir, subdivision, v);
    makeTinyPlotsOfPopulation(targetDimSpdfInRFNormSub, targetDimSpdfInRFNormErrSub, ...
            targetDimSpdfExRFNormSub, targetDimSpdfExRFNormErrSub, targetDimT, unitNamesSub, titleBase, plotFileBaseName);
        
    titleBase = sprintf('%s Cells: Exit Fixation', subdivision);
    plotFileBaseName = sprintf('%s/allSessions-%s-tinyPop-meanSpdf5-exitFixation-v%d', processedDataRootDir, subdivision, v);
    makeTinyPlotsOfPopulation(exitFixationSpdfInRFNormSub, exitFixationSpdfInRFNormErrSub, ...
            exitFixationSpdfExRFNormSub, exitFixationSpdfExRFNormErrSub, exitFixationT, unitNamesSub, titleBase, plotFileBaseName);
    
end

%% PCA on activity space by cell
superdivision = isCell & isInPulvinar;
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca([meanNormSpdfInRFAllWindowsAll(superdivision,[1:2 4:end])]);
fprintf('\n');
fprintf('PCA: %d variables, %d observations\n', size(pcaScore, 2), size(pcaScore, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', pcaPctExplained(1) + pcaPctExplained(2));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;

cols = lines(5);
cols = [cols(1,:); cols(3,:); cols(5,:); cols(2,:)];

subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
localizationSuper = localization(superdivision);
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    isInSubdivision = strcmp(localizationSuper, subdivision);
    fprintf('\t%s: %d cells\n', subdivision, sum(isInSubdivision));

    sh = scatter(pcaScore(isInSubdivision,1), pcaScore(isInSubdivision,2), 100, cols(i,:), 'MarkerFaceColor', cols(i,:));
    sh.MarkerFaceAlpha = 0.9;
end

% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% xlim([-25 120]);
% ylim([-30 60]);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');
% set(gca, 'FontWeight', 'bold');


plotFileName = sprintf('%s/allSessions-pcaByActivity-splitBySubdivision-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

print(sprintf('%s/posterFigs/allSessions-pcaByActivity-splitBySubdivision-v%d.svg', processedDataRootDir, v), '-dsvg');

%% PCA on activity space by cell
superdivision = isCell & isInPulvinar;
[pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca([meanNormSpdfInRFAllWindowsAll(superdivision,[1:2 4:end])]);
fprintf('\n');
fprintf('PCA: %d variables, %d observations\n', size(pcaScore, 2), size(pcaScore, 1));
fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', pcaPctExplained(1) + pcaPctExplained(2));

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;

cols = lines(5);
cols = [cols(1,:); cols(3,:); cols(5,:); cols(2,:)];

subdivisions = {'PM', 'PLd', 'PLv', 'PI'};
localizationSuper = localization(superdivision);
isSignificantDelaySelectivitySuper = isSignificantDelaySelectivity(superdivision);
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    isInSubdivision = strcmp(localizationSuper, subdivision);
    fprintf('\t%s: %d cells\n', subdivision, sum(isInSubdivision));

    sh = scatter(pcaScore(isInSubdivision,1), pcaScore(isInSubdivision,2), 100, cols(i,:), 'MarkerFaceColor', cols(i,:));
    sh.MarkerFaceAlpha = 0.9;
end

% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% xlim([-25 120]);
% ylim([-30 60]);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');
% set(gca, 'FontWeight', 'bold');


plotFileName = sprintf('%s/allSessions-pcaByActivity-splitBySubdivision-delaySigHighlight-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

print(sprintf('%s/posterFigs/allSessions-pcaByActivity-splitBySubdivision-delaySigHighlight-v%d.svg', processedDataRootDir, v), '-dsvg');

%% use t-sne (random position each run)

tsneVals = tsne([meanNormSpdfInRFAllWindowsAll(superdivision,[1:2 4:end])]);

figure_tr_inch(7.5, 7.5);
subaxis(1, 1, 1, 'MB', 0.14, 'MT', 0.03, 'ML', 0.16)
hold on;
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    isInSubdivision = strcmp(localizationSuper, subdivision);
    fprintf('\t%s: %d cells\n', subdivision, sum(isInSubdivision));

    sh = scatter(tsneVals(isInSubdivision,1), tsneVals(isInSubdivision,2), 100, cols(i,:), 'MarkerFaceColor', cols(i,:));
    sh.MarkerFaceAlpha = 0.9;
end


% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% xlim([-25 120]);
% ylim([-30 60]);
set(gca, 'box', 'off');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 26);
set(gca, 'FontName', 'Calibri');
% set(gca, 'FontWeight', 'bold');
