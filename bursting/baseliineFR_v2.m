clc
clear all
close all
cd('E:\grating-data-diffbranch')
nUnitsDir = dir('*.mat');

nUnits = length(nUnitsDir);
cd('E:\grating-analysis-output')

for uniti = 1:nUnits
    tic
    clearvars -except uniti nUnitsDir chanNum pulLocV firingRate allBaseline...
        meanSCInRFCount meanSCInRFBlCount meanSCOutRFCount meanSCOutRFBlCount ...
        pSCCIsum pSCCOsum meanSCAttInCount meanSCAttOutCount pSCATTsum ...
        condiInfo SPNA BRI nUnits allISIAttIn allISIAttOut perc_shortIsi_allts ...
        perc_isiAttInCA  perc_isiAttOutCA comb_isiCA perc_isiAllCA       
    
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
        
        targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
        targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
        arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
        
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
        for cuei = 1:numel(cueOnsetTrialAttIn)
            nwinr = arrayOnsetTrialAttIn(cuei);
            nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
            spikeCount.attIn(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            isiAttInCA = [isiAttInCA; diff(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))];
            frCount.attIn(cuei) = spikeCount.attIn(cuei) / round((nwinr - nwinl),3);
        end
        for cuei = 1:numel(cueOnsetTrialAttOut)
            nwinr = arrayOnsetTrialAttOut(cuei);
            nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
            spikeCount.attOut(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
            isiAttOutCA = [isiAttOutCA; diff(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))];
            frCount.attOut(cuei) = spikeCount.attOut(cuei) / round((nwinr - nwinl),3);
        end
        meanFRAttInCount(uniti) = mean(frCount.attIn);
        meanFRAttOutCount(uniti) = mean(frCount.attOut);
        [pSCATTsum(uniti)] = ranksum(frCount.attIn,frCount.attOut);
        
        perc_isiAttInCA(uniti) = (sum(isiAttInCA<.005)*100)/numel(isiAttInCA);
        perc_isiAttOutCA(uniti) = (sum(isiAttOutCA<.005)*100)/numel(isiAttOutCA);
        comb_isiCA = [isiAttInCA; isiAttOutCA];
        perc_isiAllCA(uniti) = (sum(comb_isiCA<.005)*100)/numel(comb_isiCA);
        
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
    
    
%     tempResolve = 1;
%     withinBurstTime = 5; preQuietPeriod = 100;
%     [firingRate(uniti), percClassicTwoSpikeBurst(uniti), percClassicThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLongClassic(uniti),...
%         numBurstShortClassic(uniti), frAttInPerTrial,frAttOutPerTrial,...
%         percClassicTwoSpikeBurstAttIn(uniti),percClassicTwoSpikeBurstAttOut(uniti),...
%         baselineBurst(uniti),baselineBurstAttIn(uniti),baselineBurstAttOut(uniti),...
%         postCueBurst(uniti),postCueBurstAttIn(uniti),postCueBurstAttOut(uniti),...
%         postArrayBurst(uniti), postArrayBurstAttIn(uniti),postArrayBurstAttOut(uniti),...
%         cueOnsetClassic(uniti),arrayOnsetClassic(uniti),targetDimOnsetClassic(uniti)] = histClassicAttCondi_homepc(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
%                     arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, targetDimTrialAttIn, targetDimTrialAttOut, tempResolve, withinBurstTime, preQuietPeriod, UE, startTime, endTime);

    fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
    toc
end
cd('E:\grating-analysis-output')
visResp = (pSCCIsum <= .05 & meanSCInRFCount > meanSCInRFBlCount) | ((pSCCOsum <= .05 & pSCCOsum > 0) & meanSCOutRFCount > meanSCOutRFBlCount);
attMod = pSCATTsum <= .05;
pulLocV = boolean(pulLocV);
% save('infoUnits','visResp','attMod','p*','mean*')

%%
chanNumAdj = chanNum; chanNumAdj(chanNumAdj>32) = chanNumAdj(chanNumAdj>32)-32;
figure
subplot(211)
% plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
plot([allBaseline(~pulLocV &firingRate>1).meanFRCount],chanNumAdj(~pulLocV &firingRate>1)+32,'k*')
hold on
plot([allBaseline(pulLocV &firingRate>1).meanFRCount],chanNumAdj(pulLocV &firingRate>1),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
% abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
abHistD = histogram([allBaseline(~pulLocV &firingRate>1).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram([allBaseline(pulLocV &firingRate>1).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
% plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
%     [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(~pulLocV &firingRate>1).meanFRCount]) median([allBaseline(~pulLocV &firingRate>1).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median([allBaseline(pulLocV &firingRate>1).meanFRCount]) median([allBaseline(pulLocV &firingRate>1).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

figure
subplot(211)
% plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
plot([allBaseline(~pulLocV & visResp &firingRate>1).meanFRCount],chanNumAdj(~pulLocV & visResp&firingRate>1)+32,'k*')
hold on
plot([allBaseline(pulLocV & visResp&firingRate>1).meanFRCount],chanNumAdj(pulLocV & visResp&firingRate>1),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
% abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
abHistD = histogram([allBaseline(~pulLocV & visResp&firingRate>1).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram([allBaseline(pulLocV & visResp&firingRate>1).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
% plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
%     [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(~pulLocV & visResp&firingRate>1).meanFRCount]) median([allBaseline(~pulLocV & visResp&firingRate>1).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median([allBaseline(pulLocV & visResp&firingRate>1).meanFRCount]) median([allBaseline(pulLocV & visResp&firingRate>1).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

figure
subplot(211)
% plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
plot([allBaseline(~pulLocV & attMod&firingRate>1).meanFRCount],chanNumAdj(~pulLocV & attMod&firingRate>1)+32,'k*')
hold on
plot([allBaseline(pulLocV & attMod&firingRate>1).meanFRCount],chanNumAdj(pulLocV & attMod&firingRate>1),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
% abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
abHistD = histogram([allBaseline(~pulLocV & attMod&firingRate>1).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram([allBaseline(pulLocV & attMod&firingRate>1).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
% plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
%     [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(~pulLocV & attMod&firingRate>1).meanFRCount]) median([allBaseline(~pulLocV & attMod&firingRate>1).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median([allBaseline(pulLocV & attMod&firingRate>1).meanFRCount]) median([allBaseline(pulLocV & attMod&firingRate>1).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')
%%
% figure
% histogram(perc_shortIsi_allts(firingRate>1),linspace(1,60,60))
% ylabel('Number of neurons')
% xlabel('Percentage of ISIs <5ms (%)')
% title('Examining bimodality in short-ISIs')

figure
subplot(211)
histogram(perc_isiAllCA(firingRate>1),linspace(0,60,61),'FaceColor',[0.5 0.5 0.5])
ylabel('Number of neurons')
title('Examining bimodality in short-ISIs in CA')
subplot(212)
histogram(perc_isiAttOutCA(firingRate>1),linspace(0,60,61))
hold on
histogram(perc_isiAttInCA(firingRate>1),linspace(0,60,61))
ylabel('Number of neurons')
legend({'AttOut' 'AttIn'})
xlabel('Percentage of ISIs <5ms (%)')

figure
histogram(perc_isiAllCA(pulLocV & firingRate>1),linspace(0,60,61),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
histogram(perc_isiAllCA(~pulLocV & firingRate>1),linspace(0,60,61),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
ylabel('Probability')
title('Examining bimodality in short-ISIs in CA')
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
xlabel('Percentage of ISIs <5ms (%)')

%% plot BRI index
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

% calculate the attention index
visUnitsV = find(visResp & pulLocV); visUnitsD = find(visResp & ~pulLocV);
visUnitsAll = [visUnitsV visUnitsD];
idx = 1; counter = 0;
for uni = 1:length(visUnitsAll)
    attIn_bFract(uni) = sum(allISIAttIn{visUnitsAll(uni)}<=4) / length(allISIAttIn{visUnitsAll(uni)});
    attOut_bFract(uni) = sum(allISIAttOut{visUnitsAll(uni)}<=4) / length(allISIAttOut{visUnitsAll(uni)});
    if attIn_bFract(uni)>0 && attOut_bFract(uni)>0
        frac_index(idx) = (attIn_bFract(uni)-attOut_bFract(uni)) / (attIn_bFract(uni)+attOut_bFract(uni));
        idx = idx +1;
    else
        counter = counter + 1;
    end
end
figure; histogram(frac_index,round(min(frac_index),1)-.1:.1:round(max(frac_index),1)+.1,'FaceColor',[1 0.8 0]);
title(['Burst Fraction Index n=' num2str(length(frac_index))])

idx = 1; counter = 0;
for uni = 1:length(visUnitsV)
    attIn_bFractV(uni) = sum(allISIAttIn{visUnitsV(uni)}<=4) / length(allISIAttIn{visUnitsV(uni)});
    attOut_bFractV(uni) = sum(allISIAttOut{visUnitsV(uni)}<=4) / length(allISIAttOut{visUnitsV(uni)});
    if attIn_bFractV(uni)>0 && attOut_bFractV(uni)>0
        frac_indexV(idx) = (attIn_bFractV(uni)-attOut_bFractV(uni)) / (attIn_bFractV(uni)+attOut_bFractV(uni));
        idx = idx +1;
    else
        counter = counter + 1;
    end
end
figure; histogram(frac_indexV,round(min(frac_index),1)-.1:.1:round(max(frac_index),1)+.1,'FaceColor',[1 0.8 0]);
title(['Burst Fraction Index vPul n=' num2str(length(frac_indexV))])

idx = 1; counter = 0;
for uni = 1:length(visUnitsD)
    attIn_bFractD(uni) = sum(allISIAttIn{visUnitsD(uni)}<=4) / length(allISIAttIn{visUnitsD(uni)});
    attOut_bFractD(uni) = sum(allISIAttOut{visUnitsD(uni)}<=4) / length(allISIAttOut{visUnitsD(uni)});
    if attIn_bFractD(uni)>0 && attOut_bFractD(uni)>0
        frac_indexD(idx) = (attIn_bFractD(uni)-attOut_bFractD(uni)) / (attIn_bFractD(uni)+attOut_bFractD(uni));
        idx = idx +1;
    else
        counter = counter + 1;
    end
end
figure; histogram(frac_indexD,round(min(frac_index),1)-.1:.1:round(max(frac_index),1)+.1,'FaceColor',[1 0.8 0]);
title(['Burst Fraction Index dPul n=' num2str(length(frac_indexD))])


for uni = 1:length(visUnitsV)
    figure(uni)
    histogram(allISIAttOut{visUnitsV(uni)},linspace(0,50,51),'normalization','Probability')
    hold on
    histogram(allISIAttIn{visUnitsV(uni)},linspace(0,50,51),'normalization','Probability')
end
for unii = 1:length(visUnitsD)
    figure(uni + unii)
    histogram(allISIAttOut{visUnitsD(unii)},linspace(0,50,51),'normalization','Probability')
    hold on
    histogram(allISIAttIn{visUnitsD(unii)},linspace(0,50,51),'normalization','Probability')
    title('dPul')
end





% eof