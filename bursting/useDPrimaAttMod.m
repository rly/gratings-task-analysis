clc
clear all
close all
cd('E:\grating-data-diffbranch')
nUnitsDir = dir('*.mat');

nUnits = length(nUnitsDir);
cd('E:\grating-analysis-output')

for uniti = 1:nUnits
    tic
    
    cd('E:\grating-data-diffbranch')
    load(nUnitsDir(uniti).name)
    if ismember(unitStruct.channelID,R.dPulChannels)
        pulLocV(uniti) = 0;
    elseif ismember(unitStruct.channelID,R.vPulChannels)
        pulLocV(uniti) = 1;
    end
    
    chanNum(uniti) = unitStruct.channelID;
    startTime = unitStruct.ts(1); endTime = unitStruct.ts(end);
    cueOnset = UE.cueOnset;
    baselineOnset = UE.cueOnset - 0.300;
    tsInRange = baselineOnset >= startTime & cueOnset <= endTime;
    cueOnset = cueOnset(tsInRange); baselineOnset = baselineOnset(tsInRange);
    
    spikeTimes2use = unitStruct.ts; %spikeStruct.tsTask;
    firingRate(uniti) = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));
    isi_allts = diff(spikeTimes2use);
    perc_shortIsi_allts(uniti) = (sum(isi_allts<.005)*100) / numel(spikeTimes2use);
    
    % compute FR count
    for cuei = 1:numel(baselineOnset)
        nwinr = cueOnset(cuei);
        nwinl = baselineOnset(cuei);
        baseline.spikeCount(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        baseline.frCount(cuei) = baseline.spikeCount(cuei) / round((nwinr - nwinl),3);
    end
    
    allBaseline(uniti).meanSpikeCount = mean(baseline.spikeCount);
    allBaseline(uniti).medianSpikeCount = median(baseline.spikeCount);
    allBaseline(uniti).meanFRCount = mean(baseline.frCount);
    allBaseline(uniti).medianFRCount = median(baseline.frCount);
    
    cueOnsetTrialInRF = UE.cueOnset(UE.cueLoc == 3);%UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3);
    cueOnsetTrialOutRF = UE.cueOnset(UE.cueLoc == 4);%UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 4);
    if ~isempty(cueOnsetTrialInRF)
        % compute visually responsiveness
        for cuei = 1:numel(cueOnsetTrialInRF)
            nwinr = cueOnsetTrialInRF(cuei);
            nwinl = cueOnsetTrialInRF(cuei) - 0.175;
            spikeCount.cueInRFBaseline(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            frCount.cueInRFBaseline(cuei) = spikeCount.cueInRFBaseline(cuei) / round((nwinr - nwinl),3);
            
            nwinr = cueOnsetTrialInRF(cuei)+1.000;
            nwinl = cueOnsetTrialInRF(cuei) - 0.300;
            spikeCountWhole.cueInRFBaseline(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            frCountWhole.cueInRFBaseline(cuei) = spikeCountWhole.cueInRFBaseline(cuei) / round((nwinr - nwinl),3);
        end
        for cuei = 1:numel(cueOnsetTrialInRF)
            nwinr = cueOnsetTrialInRF(cuei) + 0.200;
            nwinl = cueOnsetTrialInRF(cuei) + 0.025;
            spikeCount.cueInRF(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            frCount.cueInRF(cuei) = spikeCount.cueInRF(cuei) / round((nwinr - nwinl),3);
        end
        if ~isempty(cueOnsetTrialOutRF)
            for cuei = 1:numel(cueOnsetTrialOutRF)
                nwinr = cueOnsetTrialOutRF(cuei);
                nwinl = cueOnsetTrialOutRF(cuei) - 0.175;
                spikeCount.cueOutRFBaseline(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                frCount.cueOutRFBaseline(cuei) = spikeCount.cueOutRFBaseline(cuei) / round((nwinr - nwinl),3);
                
                nwinr = cueOnsetTrialInRF(cuei)+1.000;
                nwinl = cueOnsetTrialInRF(cuei) - 0.300;
                spikeCountWhole.cueOutRFBaseline(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                frCountWhole.cueOutRFBaseline(cuei) = spikeCountWhole.cueOutRFBaseline(cuei) / round((nwinr - nwinl),3);
            end
            for cuei = 1:numel(cueOnsetTrialOutRF)
                nwinr = cueOnsetTrialOutRF(cuei) + 0.200;
                nwinl = cueOnsetTrialOutRF(cuei) + 0.025;
                spikeCount.cueOutRF(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                frCount.cueOutRF(cuei) = spikeCount.cueOutRF(cuei) / round((nwinr - nwinl),3);
            end
            meanSCOutRFCount(uniti) = mean(spikeCount.cueOutRF);
            meanSCOutRFBlCount(uniti) = mean(spikeCount.cueOutRFBaseline);
            [pSCCOsum(uniti)] = signrank(frCount.cueOutRFBaseline,frCount.cueOutRF);
        end
        meanSCInRFCount(uniti) = mean(spikeCount.cueInRF);
        meanSCInRFBlCount(uniti) = mean(spikeCount.cueInRFBaseline);
        [pSCCIsum(uniti)] = signrank(frCount.cueInRFBaseline,frCount.cueInRF);
        
        meanFROutRFCountWhole(uniti) = mean(frCountWhole.cueOutRFBaseline);
        meanFRInRFCountWhole(uniti) = mean(frCountWhole.cueInRFBaseline);
        
        
        targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
        targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
        arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1);
        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3);
        
        tsInRangeAttOut = cueOnsetTrialAttOut >= startTime & targetDimTrialAttOut <= endTime;
        targetDimTrialAttOut = targetDimTrialAttOut(tsInRangeAttOut);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(tsInRangeAttOut);
        cueOnsetTrialAttOut = cueOnsetTrialAttOut(tsInRangeAttOut);
        tsInRangeAttIn = cueOnsetTrialAttIn >= startTime & targetDimTrialAttIn <= endTime;
        targetDimTrialAttIn = targetDimTrialAttIn(tsInRangeAttIn);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(tsInRangeAttIn);
        cueOnsetTrialAttIn = cueOnsetTrialAttIn(tsInRangeAttIn);
        
        % compute attentional modulation
        isiAttInCA = []; isiAttOutCA = [];
        binTrainAttIn = zeros(numel(cueOnsetTrialAttIn),ceil(max(arrayOnsetTrialAttIn - (cueOnsetTrialAttIn+0.200))*1000));
        baselineTrainAttIn = zeros(numel(cueOnsetTrialAttIn),300);
        for cuei = 1:numel(cueOnsetTrialAttIn)
            nwinr = arrayOnsetTrialAttIn(cuei);
            nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
            spikeCount.attIn(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            isiAttInCA = [isiAttInCA; diff(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))];
            frCount.attIn(cuei) = spikeCount.attIn(cuei) / round((nwinr - nwinl),3);
            spikeTimesInRangeTrial = [];
            spikeTimesInRangeTrial = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr) - nwinl;
            if ~isempty(spikeTimesInRangeTrial)
                spikeTimesInRangeTrial = spikeTimesInRangeTrial(round(spikeTimesInRangeTrial*1000)>0);
                binTrainAttIn(cuei,round(spikeTimesInRangeTrial*1000)) = 1;
            end
            spikeTimesInRangeBaseline = [];
            nwinl = cueOnsetTrialAttIn(cuei)-0.3;
            nwinr = cueOnsetTrialAttIn(cuei);
            spikeTimesInRangeBaseline = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr) - nwinl;
            if ~isempty(spikeTimesInRangeBaseline)
                spikeTimesInRangeBaseline = spikeTimesInRangeBaseline(round(spikeTimesInRangeBaseline*1000)>0);
                baselineTrainAttIn(cuei,round(spikeTimesInRangeBaseline*1000)) = 1;
            end
        end
        binTrainAttOut = zeros(numel(cueOnsetTrialAttOut),ceil(max(arrayOnsetTrialAttOut - (cueOnsetTrialAttOut+0.200))*1000));
        baselineTrainAttOut = zeros(numel(cueOnsetTrialAttOut),300);
        for cuei = 1:numel(cueOnsetTrialAttOut)
            nwinr = arrayOnsetTrialAttOut(cuei);
            nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
            spikeCount.attOut(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            isiAttOutCA = [isiAttOutCA; diff(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))];
            frCount.attOut(cuei) = spikeCount.attOut(cuei) / round((nwinr - nwinl),3);
            spikeTimesInRangeTrial = [];
            spikeTimesInRangeTrial = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr) - nwinl;
            if ~isempty(spikeTimesInRangeTrial)
                spikeTimesInRangeTrial = spikeTimesInRangeTrial(round(spikeTimesInRangeTrial*1000)>0);
                binTrainAttOut(cuei,round(spikeTimesInRangeTrial*1000)) = 1;
            end
            spikeTimesInRangeBaseline = [];
            nwinl = cueOnsetTrialAttOut(cuei)-0.3;
            nwinr = cueOnsetTrialAttOut(cuei);
            spikeTimesInRangeBaseline = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr) - nwinl;
            if ~isempty(spikeTimesInRangeBaseline)
                spikeTimesInRangeBaseline = spikeTimesInRangeBaseline(round(spikeTimesInRangeBaseline*1000)>0);
                baselineTrainAttOut(cuei,round(spikeTimesInRangeBaseline*1000)) = 1;
            end
        end
        meanFRAttInCount(uniti) = mean(frCount.attIn);
        meanFRAttOutCount(uniti) = mean(frCount.attOut);
        [pSCATTsum(uniti)] = ranksum(frCount.attIn,frCount.attOut);
        
        cd('E:\To Tara\')
%         BaseOn=repmat(squeeze(mean(baselineTrainAttIn,2)),1,size(binTrainAttIn,2));
%         binTrainAttIn=binTrainAttIn-BaseOn;
%         BaseOff=repmat(squeeze(mean(baselineTrainAttOut,2)),1,size(binTrainAttOut,2));
%         binTrainAttOut=binTrainAttOut-BaseOff;
        shortestCA = min([round(min(arrayOnsetTrialAttIn - (cueOnsetTrialAttIn+0.200))*1000) round(min(arrayOnsetTrialAttOut - (cueOnsetTrialAttOut+0.200))*1000)]);
        [DprimeCueDelayNoBl(uniti),pValueCueDelayNoBl(uniti)]=getDprime(mean(binTrainAttIn(:,1:shortestCA),2),mean(binTrainAttOut(:,1:shortestCA),2));
        
        BaseOn=repmat(squeeze(mean(baselineTrainAttIn,2)),1,size(binTrainAttIn,2));
        binTrainAttIn=binTrainAttIn-BaseOn;
        BaseOff=repmat(squeeze(mean(baselineTrainAttOut,2)),1,size(binTrainAttOut,2));
        binTrainAttOut=binTrainAttOut-BaseOff;
        shortestCA = min([round(min(arrayOnsetTrialAttIn - (cueOnsetTrialAttIn+0.200))*1000) round(min(arrayOnsetTrialAttOut - (cueOnsetTrialAttOut+0.200))*1000)]);
        [DprimeCueDelay(uniti),pValueCueDelay(uniti)]=getDprime(mean(binTrainAttIn(:,1:shortestCA),2),mean(binTrainAttOut(:,1:shortestCA),2));

%         % compute attentional modulation array target delay
%         isiAttInCA = []; isiAttOutCA = [];
%         binTrainAttInTdim = zeros(numel(arrayOnsetTrialAttIn),ceil(max(targetDimTrialAttIn - (arrayOnsetTrialAttIn))*1000));
%         for cuei = 1:numel(cueOnsetTrialAttIn)
%             nwinr = targetDimTrialAttIn(cuei);
%             nwinl = arrayOnsetTrialAttIn(cuei) + 0.200;
%             spikeCount.attIn(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
%             isiAttInCA = [isiAttInCA; diff(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))];
%             frCount.attIn(cuei) = spikeCount.attIn(cuei) / round((nwinr - nwinl),3);
%             spikeTimesInRangeTrial = [];
%             spikeTimesInRangeTrial = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr) - nwinl;
%             if ~isempty(spikeTimesInRangeTrial)
%                 spikeTimesInRangeTrial = spikeTimesInRangeTrial(round(spikeTimesInRangeTrial*1000)>0);
%                 binTrainAttInTdim(cuei,round(spikeTimesInRangeTrial*1000)) = 1;
%             end
%         end
%         binTrainAttOutTdim = zeros(numel(arrayOnsetTrialAttOut),ceil(max(targetDimTrialAttOut - (arrayOnsetTrialAttOut))*1000));
%         for cuei = 1:numel(arrayOnsetTrialAttOut)
%             nwinr = targetDimTrialAttOut(cuei);
%             nwinl = arrayOnsetTrialAttOut(cuei) + 0.200;
%             spikeCount.attOut(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
%             isiAttOutCA = [isiAttOutCA; diff(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))];
%             frCount.attOut(cuei) = spikeCount.attOut(cuei) / round((nwinr - nwinl),3);
%             spikeTimesInRangeTrial = [];
%             spikeTimesInRangeTrial = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr) - nwinl;
%             if ~isempty(spikeTimesInRangeTrial)
%                 spikeTimesInRangeTrial = spikeTimesInRangeTrial(round(spikeTimesInRangeTrial*1000)>0);
%                 binTrainAttOutTdim(cuei,round(spikeTimesInRangeTrial*1000)) = 1;
%             end
%         end
%         meanFRAttInCountTdim(uniti) = mean(frCount.attIn);
%         meanFRAttOutCountTdim(uniti) = mean(frCount.attOut);
%         [pSCATTsumTdim(uniti)] = ranksum(frCount.attIn,frCount.attOut);
%         
%         shortestCA = min([round(min(targetDimTrialAttIn - (arrayOnsetTrialAttIn))*1000) round(min(targetDimTrialAttOut - (arrayOnsetTrialAttOut))*1000)]);
%         [DprimeArrayDelayNoBl(uniti),pValueArrayDelayNoBl(uniti)]=getDprime(mean(binTrainAttInTdim(:,1:shortestCA),2),mean(binTrainAttOutTdim(:,1:shortestCA),2));
%         
%         BaseOn=repmat(squeeze(mean(baselineTrainAttIn,2)),1,size(binTrainAttInTdim,2));
%         binTrainAttInTdim=binTrainAttInTdim-BaseOn;
%         BaseOff=repmat(squeeze(mean(baselineTrainAttOut,2)),1,size(binTrainAttOutTdim,2));
%         binTrainAttOutTdim=binTrainAttOutTdim-BaseOff;
%         shortestCA = min([round(min(targetDimTrialAttIn - (arrayOnsetTrialAttIn))*1000) round(min(targetDimTrialAttOut - (arrayOnsetTrialAttOut))*1000)]);
%         [DprimeArrayDelay(uniti),pValueArrayDelay(uniti)]=getDprime(mean(binTrainAttInTdim(:,1:shortestCA),2),mean(binTrainAttOutTdim(:,1:shortestCA),2));

        if pValueCueDelay(uniti)<.05
            kernelSigma = .01;
            cueOnsetfun.window = [1.5 1.5]; % seconds before, after
            cueOnsetfun.spdfWindowOffset = [-1.4 1.4]; % tighter window for spdf to avoid edge effects
            cueOnsetfun = createTimeLockedSpdf(spikeTimes2use, UE.cueOnset, UE.cueOnsetByLoc, cueOnsetfun, kernelSigma, startTime, endTime);
            spdfAttIn(uniti,:) = cueOnsetfun.spdfByLoc(3,:); spdfAttOut(uniti,:) = cueOnsetfun.spdfByLoc(1,:);
        end
        
        
        cd('C:\Users\DELL\Documents\GitHub\gratings-task-analysis\bursting')
        spikesInTrialsAttIn = zeros(numel(cueOnsetTrialAttIn),300); spikesInTrialsAttOut = zeros(numel(cueOnsetTrialAttOut),300);
        isiAttOut = []; isiAttIn = [];
        % calculate SPNA and BRI
        if length(pSCCOsum) == length(pSCCIsum)
            if (pSCCIsum(uniti) <= .05 && meanSCInRFCount(uniti) > meanSCInRFBlCount(uniti)) || (pSCCOsum(uniti) <= .05 && meanSCOutRFCount(uniti) > meanSCOutRFBlCount(uniti))
                % find spikes in cue-array period
                idx = 1;
                for cuei = 1:numel(cueOnsetTrialAttIn)
                    %                 nwinr = arrayOnsetTrialAttIn(cuei);
                    nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
                    nwinr = nwinl + 0.300;
                    if ~isempty(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))
                        tmpSpikeTimes = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                        tmpSpikeTimes = (tmpSpikeTimes - nwinl) * 1000;
                        spikeInd = round(tmpSpikeTimes); spikeInd = spikeInd(spikeInd>0);
                        isiAttIn = [isiAttIn; diff(round(tmpSpikeTimes))];
                        spikesInTrialsAttIn(cuei,spikeInd) = 1;
                        spikeCountIn(idx) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                        idx = idx + 1;
                        clear tmpSpikeTimes
                        %             else
                        %                 spikeTimesInTrial(cuei).AttIn = 0;
                    end
                end
                idx = 1;
                for cuei = 1:numel(cueOnsetTrialAttOut)
                    %                 nwinr = arrayOnsetTrialAttOut(cuei);
                    nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
                    nwinr = nwinl + 0.300;
                    if ~isempty(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))
                        tmpSpikeTimes = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                        tmpSpikeTimes = (tmpSpikeTimes - nwinl) * 1000;
                        spikeInd = round(tmpSpikeTimes); spikeInd = spikeInd(spikeInd>0);
                        isiAttOut = [isiAttOut; diff(round(tmpSpikeTimes))];
                        spikesInTrialsAttOut(cuei,spikeInd) = 1;
                        spikeCountOut(idx) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                        idx = idx + 1;
                        clear tmpSpikeTimes
                        %             else
                        %                 spikeTimesInTrial(cuei).AttOut = 0;
                    end
                end
                condiInfo(uniti,:) = [numel(cueOnsetTrialAttIn) numel(cueOnsetTrialAttOut) sum(spikeCountIn) sum(spikeCountOut)];
                spikeFs = 1000;
                maxLag = spikeFs / 4; % 20ms
                if size(spikesInTrialsAttIn,1)>1 && size(spikesInTrialsAttOut,1)>1
                [SPNA(uniti).AttIn, meanAutoCorr(uniti).AttIn, meanCrossCorr(uniti).AttIn, sdCrossCorr(uniti).AttIn] = computeSPNA(spikesInTrialsAttIn, maxLag);
                [SPNA(uniti).AttOut, meanAutoCorr(uniti).AttOut, meanCrossCorr(uniti).AttOut, sdCrossCorr(uniti).AttOut] = computeSPNA(spikesInTrialsAttOut, maxLag);
                [SPNA(uniti).All, meanAutoCorr(uniti).All, meanCrossCorr(uniti).All, sdCrossCorr(uniti).All] = computeSPNA([spikesInTrialsAttIn; spikesInTrialsAttOut], maxLag);
                BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
                BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
                BRI(uniti).AttIn = mean(SPNA(uniti).AttIn(BRILagStartInd:BRILagEndInd));
                BRI(uniti).AttOut = mean(SPNA(uniti).AttOut(BRILagStartInd:BRILagEndInd));
                BRI(uniti).All = mean(SPNA(uniti).All(BRILagStartInd:BRILagEndInd));
                allISIAttIn{uniti} = isiAttIn;
                allISIAttOut{uniti} = isiAttOut;
                end
                clear spikeTimesInTrial
            end
        elseif (pSCCIsum(uniti) <= .05 && meanSCInRFCount(uniti) > meanSCInRFBlCount(uniti))
            % find spikes in cue-array period
            idx = 1;
            for cuei = 1:numel(cueOnsetTrialAttIn)
                %                 nwinr = arrayOnsetTrialAttIn(cuei);
                nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
                nwinr = nwinl + 0.300;
                if ~isempty(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))
                    tmpSpikeTimes = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                    tmpSpikeTimes = (tmpSpikeTimes - nwinl) * 1000;
                    spikeInd = round(tmpSpikeTimes); spikeInd = spikeInd(spikeInd>0);
                    isiAttIn = [isiAttIn; diff(round(tmpSpikeTimes))];
                    spikesInTrialsAttIn(cuei,spikeInd) = 1;
                    spikeCountIn(idx) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                    idx = idx + 1;
                    clear tmpSpikeTimes
                    %             else
                    %                 spikeTimesInTrial(cuei).AttIn = 0;
                end
            end
            idx = 1;
            for cuei = 1:numel(cueOnsetTrialAttOut)
                %                 nwinr = arrayOnsetTrialAttOut(cuei);
                nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
                nwinr = nwinl + 0.300;
                if ~isempty(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))
                    tmpSpikeTimes = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                    tmpSpikeTimes = (tmpSpikeTimes - nwinl) * 1000;
                    spikeInd = round(tmpSpikeTimes); spikeInd = spikeInd(spikeInd>0);
                    isiAttOut = [isiAttOut; diff(round(tmpSpikeTimes))];
                    spikesInTrialsAttOut(cuei,spikeInd) = 1;
                    spikeCountOut(idx) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                    idx = idx + 1;
                    clear tmpSpikeTimes
                    %             else
                    %                 spikeTimesInTrial(cuei).AttOut = 0;
                end
            end
            condiInfo(uniti,:) = [numel(cueOnsetTrialAttIn) numel(cueOnsetTrialAttOut) sum(spikeCountIn) sum(spikeCountOut)];
            spikeFs = 1000;
            maxLag = spikeFs / 4; % 20ms
            %             if length(spikeTimesInTrialAttIn) > 1 && length(spikeTimesInTrialAttOut) > 1
            [SPNA(uniti).AttIn, meanAutoCorr(uniti).AttIn, meanCrossCorr(uniti).AttIn, sdCrossCorr(uniti).AttIn] = computeSPNA(spikesInTrialsAttIn, maxLag);
            [SPNA(uniti).AttOut, meanAutoCorr(uniti).AttOut, meanCrossCorr(uniti).AttOut, sdCrossCorr(uniti).AttOut] = computeSPNA(spikesInTrialsAttOut, maxLag);
            [SPNA(uniti).All, meanAutoCorr(uniti).All, meanCrossCorr(uniti).All, sdCrossCorr(uniti).All] = computeSPNA([spikesInTrialsAttIn; spikesInTrialsAttOut], maxLag);
            BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
            BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
            BRI(uniti).AttIn = mean(SPNA(uniti).AttIn(BRILagStartInd:BRILagEndInd));
            BRI(uniti).AttOut = mean(SPNA(uniti).AttOut(BRILagStartInd:BRILagEndInd));
            BRI(uniti).All = mean(SPNA(uniti).All(BRILagStartInd:BRILagEndInd));
            allISIAttIn{uniti} = isiAttIn;
            allISIAttOut{uniti} = isiAttOut;
            %             end
            clear spikeTimesInTrial
        end
    end

fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
toc
end
% cd('E:\grating-analysis-output')
visResp = (pSCCIsum <= .05 & meanSCInRFCount > meanSCInRFBlCount) | ((pSCCOsum <= .05 & pSCCOsum > 0) & meanSCOutRFCount > meanSCOutRFBlCount);
attMod = pSCATTsum < .05;
pulLocV = boolean(pulLocV);
% save('infoUnits','visResp','attMod','p*','mean*')
sum(pValueCueDelay<.05)
meanFRAll = (meanFROutRFCountWhole+meanFRInRFCountWhole)/2;%(meanFRAttInCount+meanFRAttOutCount)/2;

%%
figure
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttOut(visResp&pValueCueDelay<.05,:)))
hold on
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttIn(visResp&pValueCueDelay<.05,:)))
legend({'AttOut' 'AttIn'})


figure
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttOut((visResp&pValueCueDelayNoBl<.05&meanFRAttOutCount>5),:)))
hold on
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttIn((visResp&pValueCueDelayNoBl<.05&meanFRAttOutCount>5),:)))
legend({'AttOut' 'AttIn'})


figure
subplot(311)
histogram(meanFRAttOutCount(visResp&pValueCueDelay<.05),0:1:20)
hold on
histogram(meanFRAttInCount(visResp&pValueCueDelay<.05),0:1:20)
subplot(312)
histogram(meanFRAttOutCount(pulLocV&visResp&pValueCueDelay<.05),0:1:20)
hold on
histogram(meanFRAttInCount(pulLocV&visResp&pValueCueDelay<.05),0:1:20)
subplot(313)
histogram(meanFRAttOutCount(~pulLocV&visResp&pValueCueDelay<.05),0:1:20)
hold on
histogram(meanFRAttInCount(~pulLocV&visResp&pValueCueDelay<.05),0:1:20)
[pdPulRyanFR,bdPulRyanFR] = signrank(meanFRAttOutCount(~pulLocV&visResp&pValueCueDelay<.05),meanFRAttInCount(~pulLocV&visResp&pValueCueDelay<.05))
[pvPulRyanFR,bvPulRyanFR] = signrank(meanFRAttOutCount(pulLocV&visResp&pValueCueDelay<.05),meanFRAttInCount(pulLocV&visResp&pValueCueDelay<.05))
[pPulRyanFR,bPulRyanFR] = signrank(meanFRAttOutCount(visResp&pValueCueDelay<.05),meanFRAttInCount(visResp&pValueCueDelay<.05))

[pdPulRyanFRm1,bdPulRyanFRm1] = signrank(meanFRAttOutCount(monkey&~pulLocV&visResp&pValueCueDelay<.05),meanFRAttInCount(monkey&~pulLocV&visResp&pValueCueDelay<.05))
[pvPulRyanFRm1,bvPulRyanFRm1] = signrank(meanFRAttOutCount(monkey&pulLocV&visResp&pValueCueDelay<.05),meanFRAttInCount(monkey&pulLocV&visResp&pValueCueDelay<.05))
[pPulRyanFRm1,bPulRyanFRm1] = signrank(meanFRAttOutCount(monkey&visResp&pValueCueDelay<.05),meanFRAttInCount(monkey&visResp&pValueCueDelay<.05))

[pdPulRyanFRm2,bdPulRyanFRm2] = signrank(meanFRAttOutCount(~monkey&~pulLocV&visResp&pValueCueDelay<.05),meanFRAttInCount(~monkey&~pulLocV&visResp&pValueCueDelay<.05))
[pvPulRyanFRm2,bvPulRyanFRm2] = signrank(meanFRAttOutCount(~monkey&pulLocV&visResp&pValueCueDelay<.05),meanFRAttInCount(~monkey&pulLocV&visResp&pValueCueDelay<.05))
[pPulRyanFRm2,bPulRyanFRm2] = signrank(meanFRAttOutCount(~monkey&visResp&pValueCueDelay<.05),meanFRAttInCount(~monkey&visResp&pValueCueDelay<.05))


%%
figure
subplot(211)
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttOut(visResp&pValueCueDelay<.05,:)))
hold on
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttIn(visResp&pValueCueDelay<.05,:)))
subplot(212)
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttOut(visResp&pValueCueDelayNoBl<.05,:)))
hold on
plot(linspace(cueOnsetfun.spdfWindowOffset(1),cueOnsetfun.spdfWindowOffset(2),length(cueOnsetfun.t)),mean(spdfAttIn(visResp&pValueCueDelayNoBl<.05,:)))
legend({'AttOut' 'AttIn'})




%%
figure;
histogram([BRI(visResp&firingRate>1).All],round(min([BRI(visResp&firingRate>1).All]),1)-0.1:0.1:round(max([BRI(visResp&firingRate>1).All]),1)+0.1)
title('BRI over CA 200-500ms')

% calculate BRI index
BRI_index = (([BRI(visResp&firingRate>1).AttIn] + 1) - ([BRI(visResp&firingRate>1).AttOut] + 1)) ./ (([BRI(visResp&firingRate>1).AttIn] + 1) + ([BRI(visResp&firingRate>1).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index(pulLocV(visResp&firingRate>1)), round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index(pulLocV(visResp&firingRate>1))))])
subplot(313)
hist_briD = histogram(BRI_index(~pulLocV(visResp&firingRate>1)), round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index(~pulLocV(visResp&firingRate>1))))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);


% for dprime att mod
clear BRI_index*
BRI_index = (([BRI(visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(~pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(~pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(~pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(~pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyan,bPulRyan] = signrank(BRI_index);
[pvPulRyan,bvPulRyan] = signrank(BRI_index_vPul);
[pdPulRyan,bdPulRyan] = signrank(BRI_index_dPul);

%% take dprime direction into account
clear BRI_index*
BRI_index = (([BRI(DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) - ([BRI(DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1)) ./ (([BRI(DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) + ([BRI(DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) - ([BRI(pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1)) ./ (([BRI(pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) + ([BRI(pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(~pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) - ([BRI(~pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1)) ./ (([BRI(~pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) + ([BRI(~pulLocV&DprimeCueDelayNoBl>0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyan,bPulRyan] = signrank(BRI_index);
[pvPulRyan,bvPulRyan] = signrank(BRI_index_vPul);
[pdPulRyan,bdPulRyan] = signrank(BRI_index_dPul);


clear BRI_index*
BRI_index = (([BRI(DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) - ([BRI(DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1)) ./ (([BRI(DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) + ([BRI(DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) - ([BRI(pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1)) ./ (([BRI(pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) + ([BRI(pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(~pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) - ([BRI(~pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1)) ./ (([BRI(~pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttIn] + 1) + ([BRI(~pulLocV&DprimeCueDelayNoBl<0&visResp&pValueCueDelayNoBl<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyan,bPulRyan] = signrank(BRI_index);
[pvPulRyan,bvPulRyan] = signrank(BRI_index_vPul);
[pdPulRyan,bdPulRyan] = signrank(BRI_index_dPul);














clear BRI_index*
BRI_index = (([BRI(DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(~pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(~pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(~pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(~pulLocV&DprimeCueDelay>0&visResp&pValueCueDelay<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyan,bPulRyan] = signrank(BRI_index);
[pvPulRyan,bvPulRyan] = signrank(BRI_index_vPul);
[pdPulRyan,bdPulRyan] = signrank(BRI_index_dPul);


clear BRI_index*
BRI_index = (([BRI(DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(~pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(~pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(~pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(~pulLocV&DprimeCueDelay<0&visResp&pValueCueDelay<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyan,bPulRyan] = signrank(BRI_index);
[pvPulRyan,bvPulRyan] = signrank(BRI_index_vPul);
[pdPulRyan,bdPulRyan] = signrank(BRI_index_dPul);

%%
figure
subplot(131)
histogram([BRI(visResp&pValueCueDelay<.05).AttOut],-0.3:.01:0.4)
hold on
histogram([BRI(visResp&pValueCueDelay<.05).AttIn],-0.3:.01:0.4)
subplot(132)
histogram([BRI(pulLocV&visResp&pValueCueDelay<.05).AttOut],-0.3:.01:0.4)
hold on
histogram([BRI(pulLocV&visResp&pValueCueDelay<.05).AttIn],-0.3:.01:0.4)
subplot(133)
histogram([BRI(~pulLocV&visResp&pValueCueDelay<.05).AttOut],-0.3:.01:0.4)
hold on
histogram([BRI(~pulLocV&visResp&pValueCueDelay<.05).AttIn],-0.3:.01:0.4)

% check the monkeys 
for uniti = 1:numel(nUnitsDir)
    if nUnitsDir(uniti).name(1) == 'F'
        monkey(uniti) = 0;
    else
        monkey(uniti) = 1;
    end
end


clear BRI_index*
BRI_index = (([BRI(monkey&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(monkey&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(monkey&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(monkey&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(monkey&pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(monkey&pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(monkey&pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(monkey&pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(monkey&~pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(monkey&~pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(monkey&~pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(monkey&~pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyanm1,bPulRyanm1] = signrank(BRI_index);
[pvPulRyanm1,bvPulRyanm1] = signrank(BRI_index_vPul);
[pdPulRyanm1,bdPulRyanm1] = signrank(BRI_index_dPul);

clear BRI_index*
BRI_index = (([BRI(~monkey&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(~monkey&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(~monkey&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(~monkey&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_vPul = (([BRI(~monkey&pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(~monkey&pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(~monkey&pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(~monkey&pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1));
BRI_index_dPul = (([BRI(~monkey&~pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) - ([BRI(~monkey&~pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1)) ./ (([BRI(~monkey&~pulLocV&visResp&pValueCueDelay<.05).AttIn] + 1) + ([BRI(~monkey&~pulLocV&visResp&pValueCueDelay<.05).AttOut] + 1));
figure;
subplot(311)
hist_bri = histogram(BRI_index, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in pul n=' num2str(length(BRI_index))])
subplot(312)
hist_briV = histogram(BRI_index_vPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in vPul n=' num2str(length(BRI_index_vPul))])
subplot(313)
hist_briD = histogram(BRI_index_dPul, round(min(BRI_index),2)-.01:.01:round(max(BRI_index),2)+.01,'FaceColor',[1 0.8 0]);
title(['BRI Index in dPul n=' num2str(length(BRI_index_dPul))])
subplot(311); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(312); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);
subplot(313); ylim([0 max([hist_bri.Values hist_briV.Values hist_briD.Values])]);

[pPulRyanm2,bPulRyanm2] = signrank(BRI_index);
[pvPulRyanm2,bvPulRyanm2] = signrank(BRI_index_vPul);
[pdPulRyanm2,bdPulRyanm2] = signrank(BRI_index_dPul);



%%




figure
histogram(meanFRAll(visResp),0:1:25)

figure
histogram(meanFRAll(visResp&pulLocV),0:1:25,'FaceColor',[0 1 1])
hold on
histogram(meanFRAll(visResp&~pulLocV),0:1:25,'FaceColor',[1 0 1])


figure
histogram(meanFROutRFCountWhole(visResp),0:1:25)
hold on
histogram(meanFRInRFCountWhole(visResp),0:1:25)

figure
subplot(121)
histogram(meanFROutRFCountWhole(visResp&~pulLocV),0:1:25)
hold on
histogram(meanFRInRFCountWhole(visResp&~pulLocV),0:1:25)
subplot(122)
histogram(meanFROutRFCountWhole(visResp&pulLocV),0:1:25)
hold on
histogram(meanFRInRFCountWhole(visResp&pulLocV),0:1:25)


figure
histogram(meanFROutRFCountWhole(pValueCueDelay<.05),0:1:25)
hold on
histogram(meanFRInRFCountWhole(pValueCueDelay<.05),0:1:25)

figure
subplot(121)
histogram(meanFROutRFCountWhole(pValueCueDelay<.05&~pulLocV),0:1:25)
hold on
histogram(meanFRInRFCountWhole(pValueCueDelay<.05&~pulLocV),0:1:25)
subplot(122)
histogram(meanFROutRFCountWhole(pValueCueDelay<.05&pulLocV),0:1:25)
hold on
histogram(meanFRInRFCountWhole(pValueCueDelay<.05&pulLocV),0:1:25)




figure
histogram(meanFROutRFCountWhole(visResp&pValueCueDelay<.05),0:1:25)
hold on
histogram(meanFRInRFCountWhole(visResp&pValueCueDelay<.05),0:1:25)

figure
subplot(121)
histogram(meanFROutRFCountWhole(visResp&pValueCueDelay<.05&~pulLocV),0:1:25)
hold on
histogram(meanFRInRFCountWhole(visResp&pValueCueDelay<.05&~pulLocV),0:1:25)
subplot(122)
histogram(meanFROutRFCountWhole(visResp&pValueCueDelay<.05&pulLocV),0:1:25)
hold on
histogram(meanFRInRFCountWhole(visResp&pValueCueDelay<.05&pulLocV),0:1:25)

%%
figure
subplot(311)
plot(DprimeCueDelay(visResp),DprimeArrayDelay(visResp),'o')
hold on
plot([0 0],[-1.5 2.5],'--k')
plot([-1.5 2.5],[0 0],'--k')
subplot(312)
plot(DprimeCueDelay(visResp&~pulLocV),DprimeArrayDelay(visResp&~pulLocV),'o')
hold on
plot([0 0],[-1.5 2.5],'--k')
plot([-1.5 2.5],[0 0],'--k')
subplot(313)
plot(DprimeCueDelay(visResp&pulLocV),DprimeArrayDelay(visResp&pulLocV),'o')
hold on
plot([0 0],[-1.5 2.5],'--k')
plot([-1.5 2.5],[0 0],'--k')





% eof