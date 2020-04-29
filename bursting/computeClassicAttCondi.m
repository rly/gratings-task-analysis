clc
clear all 
close all
cd('/Users/labmanager/Documents/MATLAB/Data locally SQ3')
nUnitsDir = dir('spikeTimes2use_UE_v15*');

nUnits = length(nUnitsDir);
%cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
cd('/Volumes/HDD1/BurstSep4allHDD/')
meanFrAttInPerTrial = nan(1,nUnits);
meanFrAttOutPerTrial = nan(1,nUnits);
firingRate = nan(1,nUnits);
percClassicTwoSpikeBurst = nan(1,nUnits);
percClassicThreeSpikeBurst = nan(1,nUnits);
numInterSpikeTimes = nan(1,nUnits);
numBurstLong = nan(1,nUnits);
numBurstShort = nan(1,nUnits);
percClassicTwoSpikeBurstAttIn = nan(1,nUnits);
percClassicThreeSpikeBurstAttIn = nan(1,nUnits);
numInterSpikeTimesAttIn = nan(1,nUnits);
numBurstLongAttIn = nan(1,nUnits);
numBurstShortAttIn = nan(1,nUnits);
percClassicTwoSpikeBurstAttOut = nan(1,nUnits);
percClassicThreeSpikeBurstAttOut = nan(1,nUnits);
numInterSpikeTimesAttOut = nan(1,nUnits);
numBurstLongAttOut = nan(1,nUnits);
numBurstShortAttOut = nan(1,nUnits);
%discrep = zeros(1,nUnits);

medianFRAttInCount = nan(1,nUnits);
medianFRAttInBlCount = nan(1,nUnits);
medianFRAttOutCount = nan(1,nUnits);
medianFRAttOutBlCount = nan(1,nUnits);
pCI = nan(1,nUnits);
pCO = nan(1,nUnits);

for uniti = 1:nUnits
    tic
%         if uniti == 21
%             continue
%         end
    cd('/Users/labmanager/Documents/MATLAB/Data locally SQ3')
    load(nUnitsDir(uniti).name)
    load(['spikeTimes2use_UE_' nUnitsDir(uniti).name(23:end)],'pulLoc')
%     if pulLoc ~= 'u'

    if pulLoc == 'u'
        pulLocAll(uniti) = 0;
    else
        pulLocAll(uniti) = 1;
        if pulLoc == 'v'
            pulLocV(uniti) = 1;
        else
            pulLocV(uniti) = 0;
        end
    end
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
        startTime = spikeStruct.unitStartTime; endTime = spikeStruct.unitEndTime;
        firstSpikeTimes = spikeStruct.tsTask(1);

        targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
        targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut >= startTime & targetDimTrialAttOut <= endTime);
        targetDimTrialAttOut = targetDimTrialAttOut - firstSpikeTimes;
        targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
        targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn >= startTime & targetDimTrialAttIn <= endTime);
        targetDimTrialAttIn = targetDimTrialAttIn - firstSpikeTimes;
        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut >= startTime & arrayOnsetTrialAttOut <= endTime);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut - firstSpikeTimes;
        arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn >= startTime & arrayOnsetTrialAttIn <= endTime);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn - firstSpikeTimes;
        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
        cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut >= startTime & cueOnsetTrialAttOut <= endTime);
        cueOnsetTrialAttOut = cueOnsetTrialAttOut - firstSpikeTimes;
        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
        cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn >= startTime & cueOnsetTrialAttIn <= endTime);
        cueOnsetTrialAttIn = cueOnsetTrialAttIn - firstSpikeTimes;
        
        cueOnsetTrialInRF = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3);
        cueOnsetTrialInRF = cueOnsetTrialInRF(cueOnsetTrialInRF >= startTime & cueOnsetTrialInRF <= endTime);
        cueOnsetTrialInRF = cueOnsetTrialInRF - firstSpikeTimes;
        cueOnsetTrialOutRF = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 4);
        cueOnsetTrialOutRF = cueOnsetTrialOutRF(cueOnsetTrialOutRF >= startTime & cueOnsetTrialOutRF <= endTime);
        cueOnsetTrialOutRF = cueOnsetTrialOutRF - firstSpikeTimes;

        targetDimTrialAttOut = targetDimTrialAttOut(cueOnsetTrialAttOut>0);
        targetDimTrialAttIn = targetDimTrialAttIn(cueOnsetTrialAttIn>0);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(cueOnsetTrialAttOut>0);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(cueOnsetTrialAttIn>0);
        cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut>0);
        cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn>0);
        cueOnsetTrialInRF = cueOnsetTrialInRF(cueOnsetTrialInRF>0);
        cueOnsetTrialOutRF = cueOnsetTrialOutRF(cueOnsetTrialOutRF>0);

        spikeTimes2use = spikeStruct.tsTask;
        
        spikeTimes2use = spikeTimes2use - firstSpikeTimes;
        % compute FR count
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
        
        medianDiffSCAttInCount(uniti) = median(spikeCount.cueInRF - spikeCount.cueInRFBaseline);
        medianDiffSCAttOutCount(uniti) = median(spikeCount.cueOutRF - spikeCount.cueOutRFBaseline);

        meanSCAttInCount(uniti) = mean(spikeCount.cueInRF);
        meanSCAttInBlCount(uniti) = mean(spikeCount.cueInRFBaseline);
        meanSCAttOutCount(uniti) = mean(spikeCount.cueOutRF);
        meanSCAttOutBlCount(uniti) = mean(spikeCount.cueOutRFBaseline);
        
        [pFRCI(uniti)] = signrank(frCount.cueInRFBaseline,frCount.cueInRF);
        [pFRCO(uniti)] = signrank(frCount.cueOutRFBaseline,frCount.cueOutRF);
        [pSCCI(uniti)] = signrank(spikeCount.cueInRFBaseline,spikeCount.cueInRF);
        [pSCCO(uniti)] = signrank(spikeCount.cueOutRFBaseline,spikeCount.cueOutRF);
        
        [pFRCIsum(uniti)] = ranksum(frCount.cueInRFBaseline,frCount.cueInRF);
        [pFRCOsum(uniti)] = ranksum(frCount.cueOutRFBaseline,frCount.cueOutRF);
        [pSCCIsum(uniti)] = ranksum(spikeCount.cueInRFBaseline,spikeCount.cueInRF);
        [pSCCOsum(uniti)] = ranksum(spikeCount.cueOutRFBaseline,spikeCount.cueOutRF);
        
        % calculate SPNA and BRI
        if (pSCCI(uniti) <= .05 && meanSCAttInCount(uniti) > meanSCAttInBlCount(uniti)) || (pSCCO(uniti) <= .05 && meanSCAttOutCount(uniti) > meanSCAttOutBlCount(uniti))
            % find spikes in cue-array period
            for cuei = 1:numel(cueOnsetTrialAttIn)
                nwinr = arrayOnsetTrialAttIn(cuei);
                nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
                if ~isempty(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))
                    spikeTimesInTrial(cuei).AttIn = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                else
                    spikeTimesInTrial(cuei).AttIn = 0;
                end
            end
            for cuei = 1:numel(cueOnsetTrialAttOut)
                nwinr = arrayOnsetTrialAttOut(cuei);
                nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
                if ~isempty(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr))
                    spikeTimesInTrial(cuei).AttOut = spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
                else
                    spikeTimesInTrial(cuei).AttOut = 0;
                end
            end
            spikeFs = 1000;
            maxLag = spikeFs / 100; % 10ms
            [SPNA(uniti).AttIn, meanAutoCorr(uniti).AttIn, meanCrossCorr(uniti).AttIn, sdCrossCorr(uniti).AttIn] = computeSPNA(cat(1,spikeTimesInTrial.AttIn), maxLag);
            [SPNA(uniti).AttOut, meanAutoCorr(uniti).AttOut, meanCrossCorr(uniti).AttOut, sdCrossCorr(uniti).AttOut] = computeSPNA(cat(1,spikeTimesInTrial.AttOut), maxLag);
            BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
            BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
            BRI(uniti).AttIn = mean(SPNA(uniti).AttIn(BRILagStartInd:BRILagEndInd));
            BRI(uniti).AttOut = mean(SPNA(uniti).AttOut(BRILagStartInd:BRILagEndInd));
            clear spikeTimesInTrial
        end
        tempResolve = 1;
        
        withinBurstTime = 5; preQuietPeriod = 100;
        spikeTimes2use = spikeStruct.tsTask;
%         [firingRate(uniti), percClassicTwoSpikeBurst(uniti), percClassicThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLong(uniti),...
%             numBurstShort(uniti), frAttInPerTrial,frAttOutPerTrial,...
%             percClassicTwoSpikeBurstAttIn(uniti), percClassicThreeSpikeBurstAttIn(uniti),...
%             numInterSpikeTimesAttIn(uniti), numBurstLongAttIn(uniti), numBurstShortAttIn(uniti),...
%             percClassicTwoSpikeBurstAttOut(uniti), percClassicThreeSpikeBurstAttOut(uniti),...
%             numInterSpikeTimesAttOut(uniti), numBurstLongAttOut(uniti), numBurstShortAttOut(uniti)] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
%                         arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, tempResolve);
        [firingRate(uniti), percClassicTwoSpikeBurst(uniti), percClassicThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLongClassic(uniti),...
            numBurstShortClassic(uniti), frAttInPerTrial,frAttOutPerTrial,...
            percClassicTwoSpikeBurstAttIn(uniti),percClassicTwoSpikeBurstAttOut(uniti),...
            baselineBurst(uniti),baselineBurstAttIn(uniti),baselineBurstAttOut(uniti),...
            postCueBurst(uniti),postCueBurstAttIn(uniti),postCueBurstAttOut(uniti),...
            postArrayBurst(uniti), postArrayBurstAttIn(uniti),postArrayBurstAttOut(uniti),...
            cueOnsetClassic(uniti),arrayOnsetClassic(uniti),targetDimOnsetClassic(uniti)] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, targetDimTrialAttIn, targetDimTrialAttOut, tempResolve, withinBurstTime, preQuietPeriod, UE, startTime, endTime);
        [firingRate(uniti), percWDTwoSpikeBurst(uniti), percWDThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLong(uniti),...
            numBurstShort(uniti), frAttInPerTrial,frAttOutPerTrial,...
            percWDTwoSpikeBurstAttIn(uniti),percWDTwoSpikeBurstAttOut(uniti),...
            baselineBurstWD(uniti),baselineBurstAttInWD(uniti),baselineBurstAttOutWD(uniti),...
            postCueBurstWD(uniti),postCueBurstAttInWD(uniti),postCueBurstAttOutWD(uniti),...
            postArrayBurstWD(uniti), postArrayBurstAttInWD(uniti),postArrayBurstAttOutWD(uniti),...
            cueOnsetWD(uniti),arrayOnsetWD(uniti),targetDimOnsetWD(uniti),...
            percRyanAttIn(uniti),percRyanAttOut(uniti)] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, targetDimTrialAttIn, targetDimTrialAttOut, tempResolve, withinBurstTime, UE, startTime, endTime);
      
        withinBurstTime = 20; preQuietPeriod = 100;
        [~, percClassicTwoSpikeBurstLong(uniti), percClassicThreeSpikeBurstLong(uniti), ~, ~,~, ~,~,...
            percClassicTwoSpikeBurstAttInLong(uniti),percClassicTwoSpikeBurstAttOutLong(uniti),...
            baselineBurstLong(uniti),baselineBurstAttInLong(uniti),baselineBurstAttOutLong(uniti),...
            postCueBurstLong(uniti),postCueBurstAttInLong(uniti),postCueBurstAttOutLong(uniti),...
            postArrayBurstLong(uniti), postArrayBurstAttInLong(uniti),postArrayBurstAttOutLong(uniti),...
            cueOnsetClassicLong(uniti),arrayOnsetClassicLong(uniti),targetDimOnsetClassicLong(uniti)] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, targetDimTrialAttIn, targetDimTrialAttOut, tempResolve, withinBurstTime, preQuietPeriod, UE, startTime, endTime);
        [~, percWDTwoSpikeBurstLong(uniti), percWDThreeSpikeBurstLong(uniti), ~, ~,~, ~,~,...
            percWDTwoSpikeBurstAttInLong(uniti),percWDTwoSpikeBurstAttOutLong(uniti),...
            baselineBurstWDLong(uniti),baselineBurstAttInWDLong(uniti),baselineBurstAttOutWDLong(uniti),...
            postCueBurstWDLong(uniti),postCueBurstAttInWDLong(uniti),postCueBurstAttOutWDLong(uniti),...
            postArrayBurstWDLong(uniti), postArrayBurstAttInWDLong(uniti),postArrayBurstAttOutWDLong(uniti),...
            cueOnsetWDLong(uniti),arrayOnsetWDLong(uniti),targetDimOnsetWDLong(uniti),~,~] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, targetDimTrialAttIn, targetDimTrialAttOut, tempResolve, withinBurstTime, UE, startTime, endTime);
        meanFrAttInPerTrial(uniti) = mean(frAttInPerTrial);
        meanFrAttOutPerTrial(uniti) = mean(frAttOutPerTrial);
%     else
%         continue
%     end
    clear nSpikesAttInPerTrial nSpikesAttOutPerTrial
    fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
    toc
end

pulLocAll = boolean(pulLocAll);
visResp = (pSCCI <= .05 & meanSCAttInCount > meanSCAttInBlCount) | (pSCCO <= .05 & meanSCAttOutCount > meanSCAttOutBlCount);

% calculate Ramcharan et al
nanmean(baselineBurst(pulLocAll & visResp)*2)
nanstd(baselineBurst(pulLocAll & visResp)*2)/sqrt(sum(pulLocAll & visResp))
pulLocV
pulLocV = [pulLocV 0 0 0 0]
pulLocV = boolean(pulLocV);
nanmean(baselineBurst(pulLocV & visResp)*2)
nanstd(baselineBurst(pulLocV & visResp)*2)/sqrt(sum(pulLocV & visResp))

figure
histogram(percRyanAttOut(pulLocAll&visResp),linspace(0,30,100))
hold on
histogram(percRyanAttIn(pulLocAll&visResp),linspace(0,30,100))
legend({'AttOut' 'AttIn'})

figure
histogram(percClassicTwoSpikeBurstAttOut(pulLocAll & visResp & (percClassicTwoSpikeBurstAttOut>0 | percClassicTwoSpikeBurstAttIn>0))*2,linspace(0,10,100))
hold on
histogram(percClassicTwoSpikeBurstAttIn(pulLocAll & visResp & (percClassicTwoSpikeBurstAttIn>0 | percClassicTwoSpikeBurstAttOut>0))*2,linspace(0,10,100))
title('Conservative approach cue-array')
xlabel('Burst percentage')
ylabel('Count')
legend({'AttOut' 'AttIn'})

cueOnsetWD = struct2table(cueOnsetWD); cueOnsetClassic = struct2table(cueOnsetClassic);
cueOnsetWDLong = struct2table(cueOnsetWDLong); cueOnsetClassicLong = struct2table(cueOnsetClassicLong);
arrayOnsetWD = struct2table(arrayOnsetWD); arrayOnsetClassic = struct2table(arrayOnsetClassic);
arrayOnsetWDLong = struct2table(arrayOnsetWDLong); arrayOnsetClassicLong = struct2table(arrayOnsetClassicLong);
targetDimOnsetWD = struct2table(targetDimOnsetWD); targetDimOnsetClassic = struct2table(targetDimOnsetClassic);
targetDimOnsetWDLong = struct2table(targetDimOnsetWDLong); targetDimOnsetClassicLong = struct2table(targetDimOnsetClassicLong);

use2norm = baselineBurst(pulLocAll & visResp) + postCueBurst(pulLocAll & visResp) + postArrayBurst(pulLocAll & visResp);
use2normLong = baselineBurstLong(pulLocAll & visResp) + postCueBurstLong(pulLocAll & visResp) + postArrayBurstLong(pulLocAll & visResp);
use2normWD = baselineBurstWD(pulLocAll & visResp) + postCueBurstWD(pulLocAll & visResp) + postArrayBurstWD(pulLocAll & visResp);
use2normWDLong = baselineBurstWDLong(pulLocAll & visResp) + postCueBurstWDLong(pulLocAll & visResp) + postArrayBurstWDLong(pulLocAll & visResp);

baselineBurstNorm = baselineBurst(pulLocAll & visResp) ./ use2norm;
postCueBurstNorm = postCueBurst(pulLocAll & visResp) ./ use2norm;
postArrayBurstNorm = postArrayBurst(pulLocAll & visResp) ./ use2norm;
baselineBurstLongNorm = baselineBurstLong(pulLocAll & visResp) ./ use2normLong;
postCueBurstLongNorm = postCueBurstLong(pulLocAll & visResp) ./ use2normLong;
postArrayBurstLongNorm = postArrayBurstLong(pulLocAll & visResp) ./ use2normLong;
baselineBurstWDNorm = baselineBurstWD(pulLocAll & visResp) ./ use2normWD;
postCueBurstWDNorm = postCueBurstWD(pulLocAll & visResp) ./ use2normWD;
postArrayBurstWDNorm = postArrayBurstWD(pulLocAll & visResp) ./ use2normWD;
baselineBurstWDLongNorm = baselineBurstWDLong(pulLocAll & visResp) ./ use2normWDLong;
postCueBurstWDLongNorm = postCueBurstWDLong(pulLocAll & visResp) ./ use2normWDLong;
postArrayBurstWDLongNorm = postArrayBurstWDLong(pulLocAll & visResp) ./ use2normWDLong;

summedBurstWD = baselineBurstWD(pulLocAll & visResp) + postCueBurstWD(pulLocAll & visResp) + postArrayBurstWD(pulLocAll & visResp);
sum((summedBurstWD/3)>1) * 100 / sum(pulLocAll & visResp)

figure
errorbar(1:3,[nanmean(baselineBurstWDNorm),nanmean(postCueBurstWDNorm),nanmean(postArrayBurstWDNorm)],...
    [nanstd(baselineBurstWDNorm)/sqrt(sum(~isnan(baselineBurstWDNorm))),...
    nanstd(postCueBurstWDNorm)/sqrt(sum(~isnan(postCueBurstWDNorm))),...
    nanstd(postArrayBurstWDNorm)/sqrt(sum(~isnan(postArrayBurstWDNorm)))],'k')
hold on
errorbar(1.02:3.02,[nanmean(baselineBurstNorm),nanmean(postCueBurstNorm),nanmean(postArrayBurstNorm)],...
    [nanstd(baselineBurstNorm)/sqrt(sum(~isnan(baselineBurstNorm))),...
    nanstd(postCueBurstNorm)/sqrt(sum(~isnan(postCueBurstNorm))),...
    nanstd(postArrayBurstNorm)/sqrt(sum(~isnan(postArrayBurstNorm)))],'m')
errorbar(1.04:3.04,[nanmean(baselineBurstWDLongNorm),nanmean(postCueBurstWDLongNorm),nanmean(postArrayBurstWDLongNorm)],...
    [nanstd(baselineBurstWDLongNorm)/sqrt(sum(~isnan(baselineBurstWDLongNorm))),...
    nanstd(postCueBurstWDLongNorm)/sqrt(sum(~isnan(postCueBurstWDLongNorm))),...
    nanstd(postArrayBurstWDLongNorm)/sqrt(sum(~isnan(postArrayBurstWDLongNorm)))],'--k')
errorbar(1.06:3.06,[nanmean(baselineBurstLongNorm),nanmean(postCueBurstLongNorm),nanmean(postArrayBurstLongNorm)],...
    [nanstd(baselineBurstLongNorm)/sqrt(sum(~isnan(baselineBurstLongNorm))),...
    nanstd(postCueBurstLongNorm)/sqrt(sum(~isnan(postCueBurstLongNorm))),...
    nanstd(postArrayBurstLongNorm)/sqrt(sum(~isnan(postArrayBurstLongNorm)))],'--m')
title('replication of WD plot VisResp Pul')
xticks(1:3)
xticklabels({'Baseline' 'PreArray' 'PreTD'})
xlim([0.5 3.5])
legend({'<=5' '<=5 + 100ms' '<=20ms' '<=20 + 100ms'})

kernelSigma = .01;
psthWindow = [0.7 0.7];
nTime = fix(10*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
psthT = linspace(0, sum(psthWindow), nTime); % starts at 0
t = psthT - psthWindow(1);

figure
boundedline(t,mean(cueOnsetWD.spdf(pulLocAll & visResp,:)),std(cueOnsetWD.spdf(pulLocAll & visResp,:))/sqrt(sum(pulLocAll & visResp)),'alpha','c')
hold on
boundedline(t,mean(arrayOnsetWD.spdf(pulLocAll & visResp,:)),std(arrayOnsetWD.spdf(pulLocAll & visResp,:))/sqrt(sum(pulLocAll & visResp)),'alpha','y')
boundedline(t,mean(targetDimOnsetWD.spdf(pulLocAll & visResp,:)),std(targetDimOnsetWD.spdf(pulLocAll & visResp,:))/sqrt(sum(pulLocAll & visResp)),'alpha','g')
xlim([-0.3 0])
xlabel('Time in s')
% legend({'' 'PreCue' '' '' 'PreArray' 'PreTargetDim'})
title('Womelsdorf criterion')

figure
boundedline(t,mean(cueOnsetClassic.spdf(pulLocAll & visResp,:)),std(cueOnsetClassic.spdf(pulLocAll & visResp,:))/sqrt(sum(pulLocAll & visResp)),'alpha','c')
hold on
boundedline(t,mean(arrayOnsetClassic.spdf(pulLocAll & visResp,:)),std(arrayOnsetClassic.spdf(pulLocAll & visResp,:))/sqrt(sum(pulLocAll & visResp)),'alpha','y')
boundedline(t,mean(targetDimOnsetClassic.spdf(pulLocAll & visResp,:)),std(targetDimOnsetClassic.spdf(pulLocAll & visResp,:))/sqrt(sum(pulLocAll & visResp)),'alpha','g')
xlim([-0.3 0])
xlabel('Time in s')
% legend({'PreCue' 'PreArray' 'PreTargetDim'})
title('Conservative criterion')

% report number of burst per unit


cueOnsetSpdf = cueOnsetClassic{(pulLocAll & visResp),{'spdfByLoc'}};
cueOnsetSpdf = permute(cat(3,cueOnsetSpdf{:}),[3 1 2]);
arrayOnsetSpdf = arrayOnsetClassic{(pulLocAll & visResp),{'spdfByLoc'}};
arrayOnsetSpdf = permute(cat(3,arrayOnsetSpdf{:}),[3 1 2]);
targetDimOnsetSpdf = targetDimOnsetClassic{(pulLocAll & visResp),{'spdfByLoc'}};
targetDimOnsetSpdf = permute(cat(3,targetDimOnsetSpdf{:}),[3 1 2]);

figure
boundedline(t,squeeze(mean(cueOnsetSpdf(:,1,:),1))',squeeze(std(cueOnsetSpdf(:,1,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha')
hold on
boundedline(t,squeeze(mean(cueOnsetSpdf(:,3,:),1))',squeeze(std(cueOnsetSpdf(:,3,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha','r')
xlim([-0.3 0])
title('PreCue')

figure
boundedline(t,squeeze(mean(arrayOnsetSpdf(:,1,:),1))',squeeze(std(arrayOnsetSpdf(:,1,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha')
hold on
boundedline(t,squeeze(mean(arrayOnsetSpdf(:,3,:),1))',squeeze(std(arrayOnsetSpdf(:,3,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha','r')
xlim([-0.3 0])
title('PreArray')

figure
boundedline(t,squeeze(mean(targetDimOnsetSpdf(:,1,:),1))',squeeze(std(targetDimOnsetSpdf(:,1,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha')
hold on
boundedline(t,squeeze(mean(targetDimOnsetSpdf(:,3,:),1))',squeeze(std(targetDimOnsetSpdf(:,3,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha','r')
xlim([-0.3 0])
title('PreTargetDim')

cueOnsetSpdf = cueOnsetWD{(pulLocAll & visResp),{'spdfByLoc'}};
cueOnsetSpdf = permute(cat(3,cueOnsetSpdf{:}),[3 1 2]);
arrayOnsetSpdf = arrayOnsetWD{(pulLocAll & visResp),{'spdfByLoc'}};
arrayOnsetSpdf = permute(cat(3,arrayOnsetSpdf{:}),[3 1 2]);
targetDimOnsetSpdf = targetDimOnsetWD{(pulLocAll & visResp),{'spdfByLoc'}};
targetDimOnsetSpdf = permute(cat(3,targetDimOnsetSpdf{:}),[3 1 2]);

figure
boundedline(t,squeeze(mean(cueOnsetSpdf(:,1,:),1))',squeeze(std(cueOnsetSpdf(:,1,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha')
hold on
boundedline(t,squeeze(mean(cueOnsetSpdf(:,3,:),1))',squeeze(std(cueOnsetSpdf(:,3,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha','r')
xlim([-0.3 0])
title('PreCue')

figure
boundedline(t,squeeze(mean(arrayOnsetSpdf(:,1,:),1))',squeeze(std(arrayOnsetSpdf(:,1,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha')
hold on
boundedline(t,squeeze(mean(arrayOnsetSpdf(:,3,:),1))',squeeze(std(arrayOnsetSpdf(:,3,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha','r')
xlim([-0.3 0])
title('PreArray')

figure
boundedline(t,squeeze(mean(targetDimOnsetSpdf(:,1,:),1))',squeeze(std(targetDimOnsetSpdf(:,1,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha')
hold on
boundedline(t,squeeze(mean(targetDimOnsetSpdf(:,3,:),1))',squeeze(std(targetDimOnsetSpdf(:,3,:),[],1)/sqrt(sum(pulLocAll & visResp)))','alpha','r')
xlim([-0.3 0])
title('PreTargetDim')


figure
plot(t, mean(cueOnsetWD.spdf(pulLocAll & visResp,:)),'m')
hold on
plot(t,mean(cueOnsetClassic.spdf(pulLocAll & visResp,:)),'k')
plot(t,mean(cueOnsetWDLong.spdf(pulLocAll & visResp,:)),'--m')
plot(t,mean(cueOnsetClassicLong.spdf(pulLocAll & visResp,:)),'--k')
legend({'<=5' '<=5 + 100ms' '<=20ms' '<=20 + 100ms'})
title('Cue-locked PSTH of Bursts')

figure
plot(t, mean(arrayOnsetWD.spdf(pulLocAll & visResp,:)),'m')
hold on
plot(t,mean(arrayOnsetClassic.spdf(pulLocAll & visResp,:)),'k')
plot(t,mean(arrayOnsetWDLong.spdf(pulLocAll & visResp,:)),'--m')
plot(t,mean(arrayOnsetClassicLong.spdf(pulLocAll & visResp,:)),'--k')
legend({'<=5' '<=5 + 100ms' '<=20ms' '<=20 + 100ms'})
title('Array-locked PSTH of Bursts')

figure
plot(t, mean(targetDimOnsetWD.spdf(pulLocAll & visResp,:)),'m')
hold on
plot(t,mean(targetDimOnsetClassic.spdf(pulLocAll & visResp,:)),'k')
plot(t,mean(targetDimOnsetWDLong.spdf(pulLocAll & visResp,:)),'--m')
plot(t,mean(targetDimOnsetClassicLong.spdf(pulLocAll & visResp,:)),'--k')
legend({'<=5' '<=5 + 100ms' '<=20ms' '<=20 + 100ms'})
title('TargetDim-locked PSTH of Bursts')

% plot BRI 
figure
histogram(BRI.AttIn(pulLocAll&VisResp)

% histogram of firing rate during task and for attention conditions cue-array
% period
figure
subplot(211)
histogram(firingRate(pulLocAll & visResp),linspace(0,30,60))
ylabel('Count')
title('FR based on all spikes')
subplot(212)
histogram(meanFrAttOutPerTrial(pulLocAll & visResp),linspace(0,30,60))
hold on
histogram(meanFrAttInPerTrial(pulLocAll & visResp),linspace(0,30,60))
legend({'AttOut' 'AttIn'})
title('FR in cue-array')
xlabel('Frequency in Hz')
ylabel('Count')

firingRateVR = firingRate(pulLocAll & visResp);
meanFrAttOutPerTrialVR = meanFrAttOutPerTrial(pulLocAll & visResp);
meanFrAttInPerTrialVR = meanFrAttInPerTrial(pulLocAll & visResp);
figure
subplot(211)
histogram(firingRateVR((summedBurstWD/3)>1),linspace(0,30,60))
ylabel('Count')
title('FR based on all spikes')
subplot(212)
histogram(meanFrAttOutPerTrialVR((summedBurstWD/3)>1),linspace(0,30,60))
hold on
histogram(meanFrAttInPerTrialVR((summedBurstWD/3)>1),linspace(0,30,60))
legend({'AttOut' 'AttIn'})
title('FR in cue-array')
xlabel('Frequency in Hz')
ylabel('Count')

figure
errorbar(1:3,[nanmean(baselineBurst(pulLocAll & visResp)),nanmean(postCueBurst(pulLocAll & visResp)),nanmean(postArrayBurst(pulLocAll & visResp))],...
    [nanstd(baselineBurst(pulLocAll & visResp))/sqrt(sum(~isnan(baselineBurst(pulLocAll & visResp)))),...
    nanstd(postCueBurst(pulLocAll & visResp))/sqrt(sum(~isnan(postCueBurst(pulLocAll & visResp)))),...
    nanstd(postArrayBurst(pulLocAll & visResp))/sqrt(sum(~isnan(postArrayBurst(pulLocAll & visResp))))],'b')
hold on
errorbar(1:3,[nanmean(baselineBurstLong(pulLocAll & visResp)),nanmean(postCueBurstLong(pulLocAll & visResp)),nanmean(postArrayBurstLong(pulLocAll & visResp))],...
    [nanstd(baselineBurstLong(pulLocAll & visResp))/sqrt(sum(~isnan(baselineBurstLong(pulLocAll & visResp)))),...
    nanstd(postCueBurstLong(pulLocAll & visResp))/sqrt(sum(~isnan(postCueBurstLong(pulLocAll & visResp)))),...
    nanstd(postArrayBurstLong(pulLocAll & visResp))/sqrt(sum(~isnan(postArrayBurstLong(pulLocAll & visResp))))],'--b')
errorbar(1:3,[nanmean(baselineBurstWD(pulLocAll & visResp)),nanmean(postCueBurstWD(pulLocAll & visResp)),nanmean(postArrayBurstWD(pulLocAll & visResp))],...
    [nanstd(baselineBurstWD(pulLocAll & visResp))/sqrt(sum(~isnan(baselineBurstWD(pulLocAll & visResp)))),...
    nanstd(postCueBurstWD(pulLocAll & visResp))/sqrt(sum(~isnan(postCueBurstWD(pulLocAll & visResp)))),...
    nanstd(postArrayBurstWD(pulLocAll & visResp))/sqrt(sum(~isnan(postArrayBurstWD(pulLocAll & visResp))))],'r')
errorbar(1:3,[nanmean(baselineBurstWDLong(pulLocAll & visResp)),nanmean(postCueBurstWDLong(pulLocAll & visResp)),nanmean(postArrayBurstWDLong(pulLocAll & visResp))],...
    [nanstd(baselineBurstWDLong(pulLocAll & visResp))/sqrt(sum(~isnan(baselineBurstWDLong(pulLocAll & visResp)))),...
    nanstd(postCueBurstWDLong(pulLocAll & visResp))/sqrt(sum(~isnan(postCueBurstWDLong(pulLocAll & visResp)))),...
    nanstd(postArrayBurstWDLong(pulLocAll & visResp))/sqrt(sum(~isnan(postArrayBurstWDLong(pulLocAll & visResp))))],'--r')


errorbar(1:3,[mean(percClassicTwoSpikeBurst) mean(percClassicTwoSpikeBurstAttIn),...
    mean(percClassicTwoSpikeBurstAttOut)],[std(percClassicTwoSpikeBurst)/sqrt(length(percClassicTwoSpikeBurst)),...
    std(percClassicTwoSpikeBurstAttIn)/sqrt(length(percClassicTwoSpikeBurst)),std(percClassicTwoSpikeBurstAttOut)/sqrt(length(percClassicTwoSpikeBurst))],'b')
hold on
errorbar(1:3,[mean(percWDTwoSpikeBurst) mean(percWDTwoSpikeBurstAttIn),...
    mean(percWDTwoSpikeBurstAttOut)],[std(percWDTwoSpikeBurst)/sqrt(length(percWDTwoSpikeBurst)),...
    std(percWDTwoSpikeBurstAttIn)/sqrt(length(percWDTwoSpikeBurst)),std(percWDTwoSpikeBurstAttOut)/sqrt(length(percWDTwoSpikeBurst))],'r')
errorbar(1:3,[mean(percClassicTwoSpikeBurstLong) mean(percClassicTwoSpikeBurstAttInLong),...
    mean(percClassicTwoSpikeBurstAttOutLong)],[std(percClassicTwoSpikeBurstLong)/sqrt(length(percClassicTwoSpikeBurstLong)),...
    std(percClassicTwoSpikeBurstAttInLong)/sqrt(length(percClassicTwoSpikeBurstLong)),std(percClassicTwoSpikeBurstAttOutLong)/sqrt(length(percClassicTwoSpikeBurstLong))],'--b')
errorbar(1:3,[mean(percWDTwoSpikeBurstLong) mean(percWDTwoSpikeBurstAttInLong),...
    mean(percWDTwoSpikeBurstAttOutLong)],[std(percWDTwoSpikeBurstLong)/sqrt(length(percWDTwoSpikeBurstLong)),...
    std(percWDTwoSpikeBurstAttInLong)/sqrt(length(percWDTwoSpikeBurstLong)),std(percWDTwoSpikeBurstAttOutLong)/sqrt(length(percWDTwoSpikeBurstLong))],'--r')
title('replication of WD plot ALL UNITS')
xticks(1:3)
xticklabels({'All' 'AttIn' 'AttOut'})
xlim([0.5 3.5])
legend({'<=5 + 100ms' '<=5ms' '<=20 + 100ms' '<=20ms'})
    
figure
errorbar(1:3,[mean(percClassicTwoSpikeBurst(boolean(pulLocAll))) mean(percClassicTwoSpikeBurstAttIn(boolean(pulLocAll))),...
    mean(percClassicTwoSpikeBurstAttOut(boolean(pulLocAll)))],[std(percClassicTwoSpikeBurst(boolean(pulLocAll)))/sqrt(length(percClassicTwoSpikeBurst(boolean(pulLocAll)))),...
    std(percClassicTwoSpikeBurstAttIn(boolean(pulLocAll)))/sqrt(length(percClassicTwoSpikeBurst(boolean(pulLocAll)))),std(percClassicTwoSpikeBurstAttOut(boolean(pulLocAll)))/sqrt(length(percClassicTwoSpikeBurst(boolean(pulLocAll))))],'b')
hold on
errorbar(1:3,[mean(percWDTwoSpikeBurst(boolean(pulLocAll))) mean(percWDTwoSpikeBurstAttIn(boolean(pulLocAll))),...
    mean(percWDTwoSpikeBurstAttOut(boolean(pulLocAll)))],[std(percWDTwoSpikeBurst(boolean(pulLocAll)))/sqrt(length(percWDTwoSpikeBurst(boolean(pulLocAll)))),...
    std(percWDTwoSpikeBurstAttIn(boolean(pulLocAll)))/sqrt(length(percWDTwoSpikeBurst(boolean(pulLocAll)))),std(percWDTwoSpikeBurstAttOut(boolean(pulLocAll)))/sqrt(length(percWDTwoSpikeBurst(boolean(pulLocAll))))],'r')
errorbar(1:3,[mean(percClassicTwoSpikeBurstLong(boolean(pulLocAll))) mean(percClassicTwoSpikeBurstAttInLong(boolean(pulLocAll))),...
    mean(percClassicTwoSpikeBurstAttOutLong(boolean(pulLocAll)))],[std(percClassicTwoSpikeBurstLong(boolean(pulLocAll)))/sqrt(length(percClassicTwoSpikeBurstLong(boolean(pulLocAll)))),...
    std(percClassicTwoSpikeBurstAttInLong(boolean(pulLocAll)))/sqrt(length(percClassicTwoSpikeBurstLong(boolean(pulLocAll)))),std(percClassicTwoSpikeBurstAttOutLong(boolean(pulLocAll)))/sqrt(length(percClassicTwoSpikeBurstLong(boolean(pulLocAll))))],'--b')
errorbar(1:3,[mean(percWDTwoSpikeBurstLong(boolean(pulLocAll))) mean(percWDTwoSpikeBurstAttInLong(boolean(pulLocAll))),...
    mean(percWDTwoSpikeBurstAttOutLong(boolean(pulLocAll)))],[std(percWDTwoSpikeBurstLong(boolean(pulLocAll)))/sqrt(length(percWDTwoSpikeBurstLong(boolean(pulLocAll)))),...
    std(percWDTwoSpikeBurstAttInLong(boolean(pulLocAll)))/sqrt(length(percWDTwoSpikeBurstLong(boolean(pulLocAll)))),std(percWDTwoSpikeBurstAttOutLong(boolean(pulLocAll)))/sqrt(length(percWDTwoSpikeBurstLong(boolean(pulLocAll))))],'--r')
title('replication of WD plot PUL UNITS')
xticks(1:3)
xticklabels({'All' 'AttIn' 'AttOut'})
xlim([0.5 3.5])
legend({'<=5 + 100ms' '<=5ms' '<=20 + 100ms' '<=20ms'})

figure
errorbar(1:3,[mean(percClassicTwoSpikeBurst(boolean(pulLocAll)&visResp)) mean(percClassicTwoSpikeBurstAttIn(boolean(pulLocAll)&visResp)),...
    mean(percClassicTwoSpikeBurstAttOut(boolean(pulLocAll)&visResp))],[std(percClassicTwoSpikeBurst(boolean(pulLocAll)&visResp))/sqrt(length(percClassicTwoSpikeBurst(boolean(pulLocAll)&visResp))),...
    std(percClassicTwoSpikeBurstAttIn(boolean(pulLocAll)&visResp))/sqrt(length(percClassicTwoSpikeBurst(boolean(pulLocAll)&visResp))),std(percClassicTwoSpikeBurstAttOut(boolean(pulLocAll)&visResp))/sqrt(length(percClassicTwoSpikeBurst(boolean(pulLocAll)&visResp)))],'b')
hold on
errorbar(1:3,[mean(percWDTwoSpikeBurst(boolean(pulLocAll)&visResp)) mean(percWDTwoSpikeBurstAttIn(boolean(pulLocAll)&visResp)),...
    mean(percWDTwoSpikeBurstAttOut(boolean(pulLocAll)&visResp))],[std(percWDTwoSpikeBurst(boolean(pulLocAll)&visResp))/sqrt(length(percWDTwoSpikeBurst(boolean(pulLocAll)&visResp))),...
    std(percWDTwoSpikeBurstAttIn(boolean(pulLocAll)&visResp))/sqrt(length(percWDTwoSpikeBurst(boolean(pulLocAll)&visResp))),std(percWDTwoSpikeBurstAttOut(boolean(pulLocAll)&visResp))/sqrt(length(percWDTwoSpikeBurst(boolean(pulLocAll)&visResp)))],'r')
errorbar(1:3,[mean(percClassicTwoSpikeBurstLong(boolean(pulLocAll)&visResp)) mean(percClassicTwoSpikeBurstAttInLong(boolean(pulLocAll)&visResp)),...
    mean(percClassicTwoSpikeBurstAttOutLong(boolean(pulLocAll)&visResp))],[std(percClassicTwoSpikeBurstLong(boolean(pulLocAll)&visResp))/sqrt(length(percClassicTwoSpikeBurstLong(boolean(pulLocAll)&visResp))),...
    std(percClassicTwoSpikeBurstAttInLong(boolean(pulLocAll)&visResp))/sqrt(length(percClassicTwoSpikeBurstLong(boolean(pulLocAll)&visResp))),std(percClassicTwoSpikeBurstAttOutLong(boolean(pulLocAll)&visResp))/sqrt(length(percClassicTwoSpikeBurstLong(boolean(pulLocAll)&visResp)))],'--b')
errorbar(1:3,[mean(percWDTwoSpikeBurstLong(boolean(pulLocAll)&visResp)) mean(percWDTwoSpikeBurstAttInLong(boolean(pulLocAll)&visResp)),...
    mean(percWDTwoSpikeBurstAttOutLong(boolean(pulLocAll)&visResp))],[std(percWDTwoSpikeBurstLong(boolean(pulLocAll)&visResp))/sqrt(length(percWDTwoSpikeBurstLong(boolean(pulLocAll)&visResp))),...
    std(percWDTwoSpikeBurstAttInLong(boolean(pulLocAll)&visResp))/sqrt(length(percWDTwoSpikeBurstLong(boolean(pulLocAll)&visResp))),std(percWDTwoSpikeBurstAttOutLong(boolean(pulLocAll)&visResp))/sqrt(length(percWDTwoSpikeBurstLong(boolean(pulLocAll)&visResp)))],'--r')
title('replication of WD plot PUL VisResp UNITS')
xticks(1:3)
xticklabels({'All' 'AttIn' 'AttOut'})
xlim([0.5 3.5])
legend({'<=5 + 100ms' '<=5ms' '<=20 + 100ms' '<=20ms'})


figure
subplot(211)
histogram(percClassicThreeSpikeBurst,linspace(0,1,110))
title('Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicThreeSpikeBurstAttOut,linspace(0,1,110))
hold on
histogram(percClassicThreeSpikeBurstAttIn,linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percClassicThreeSpikeBurst(percClassicThreeSpikeBurst>0),linspace(0,1,110))
title('Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicThreeSpikeBurstAttOut(percClassicThreeSpikeBurstAttOut>0),linspace(0,1,110))
hold on
histogram(percClassicThreeSpikeBurstAttIn(percClassicThreeSpikeBurstAttIn>0),linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percClassicTwoSpikeBurst,linspace(0,5,110))
title('Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicTwoSpikeBurstAttOut,linspace(0,5,110))
hold on
histogram(percClassicTwoSpikeBurstAttIn,linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percClassicTwoSpikeBurst(percClassicTwoSpikeBurst>0),linspace(0,5,110))
title('Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicTwoSpikeBurstAttOut(percClassicTwoSpikeBurstAttOut>0),linspace(0,5,110))
hold on
histogram(percClassicTwoSpikeBurstAttIn(percClassicTwoSpikeBurstAttIn>0),linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})


figure
subplot(211)
histogram(firingRate,linspace(0,30,220))
title('Firing Rate')
ylabel('Count')
subplot(212)
histogram(meanFrAttOutPerTrial,linspace(0,30,220))
hold on
histogram(meanFrAttInPerTrial,linspace(0,30,220))
ylabel('Count'); xlabel('Frequency in Hz')
legend({'AttOut' 'AttIn'})

% visually responsive for cue location with higher FR post cue compared to
% baseline

figure
subplot(211)
histogram(percClassicThreeSpikeBurst(visResp),linspace(0,1,110))
title('VisResp Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicThreeSpikeBurstAttOut(visResp),linspace(0,1,110))
hold on
histogram(percClassicThreeSpikeBurstAttIn(visResp),linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percClassicThreeSpikeBurst(percClassicThreeSpikeBurst>0 & visResp),linspace(0,1,110))
title('VisResp Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicThreeSpikeBurstAttOut(percClassicThreeSpikeBurstAttOut>0 & visResp),linspace(0,1,110))
hold on
histogram(percClassicThreeSpikeBurstAttIn(percClassicThreeSpikeBurstAttIn>0 & visResp),linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percClassicTwoSpikeBurst(visResp),linspace(0,5,110))
title('VisResp Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicTwoSpikeBurstAttOut(visResp),linspace(0,5,110))
hold on
histogram(percClassicTwoSpikeBurstAttIn(visResp),linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percClassicTwoSpikeBurst(percClassicTwoSpikeBurst>0 & visResp),linspace(0,5,110))
title('VisResp Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percClassicTwoSpikeBurstAttOut(percClassicTwoSpikeBurstAttOut>0 & visResp),linspace(0,5,110))
hold on
histogram(percClassicTwoSpikeBurstAttIn(percClassicTwoSpikeBurstAttIn>0 & visResp),linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})


figure
subplot(211)
histogram(firingRate(visResp),linspace(0,30,220))
title('VisResp Firing Rate')
ylabel('Count')
subplot(212)
histogram(meanFrAttOutPerTrial(visResp),linspace(0,30,220))
hold on
histogram(meanFrAttInPerTrial(visResp),linspace(0,30,220))
ylabel('Count'); xlabel('Frequency in Hz')
legend({'AttOut' 'AttIn'})



%% Copy Womelsdorf method
clc
clear all 
close all
cd('/Users/labmanager/Documents/MATLAB/Data locally SQ3')
nUnitsDir = dir('spikeTimes2use_UE_v14*');

nUnits = length(nUnitsDir);
%cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
cd('/Volumes/HDD1/BurstSep4allHDD/')
meanFrAttInPerTrial = nan(1,nUnits);
meanFrAttOutPerTrial = nan(1,nUnits);
firingRate = nan(1,nUnits);
percWDTwoSpikeBurst = nan(1,nUnits);
percWDThreeSpikeBurst = nan(1,nUnits);
numInterSpikeTimes = nan(1,nUnits);
numBurstLong = nan(1,nUnits);
numBurstShort = nan(1,nUnits);
percWDTwoSpikeBurstAttIn = nan(1,nUnits);
percWDThreeSpikeBurstAttIn = nan(1,nUnits);
numInterSpikeTimesAttIn = nan(1,nUnits);
numBurstLongAttIn = nan(1,nUnits);
numBurstShortAttIn = nan(1,nUnits);
percWDTwoSpikeBurstAttOut = nan(1,nUnits);
percWDThreeSpikeBurstAttOut = nan(1,nUnits);
numInterSpikeTimesAttOut = nan(1,nUnits);
numBurstLongAttOut = nan(1,nUnits);
numBurstShortAttOut = nan(1,nUnits);
%discrep = zeros(1,nUnits);

medianFRAttInCount = nan(1,nUnits);
medianFRAttInBlCount = nan(1,nUnits);
medianFRAttOutCount = nan(1,nUnits);
medianFRAttOutBlCount = nan(1,nUnits);
pCI = nan(1,nUnits);
pCO = nan(1,nUnits);

for uniti = 1:nUnits
%         if uniti == 21
%             continue
%         end
    cd('/Users/labmanager/Documents/MATLAB/Data locally SQ3')
    load(nUnitsDir(uniti).name)
    load(['spikeTimes2use_UE_' nUnitsDir(uniti).name(23:end)],'pulLoc')
    if pulLoc ~= 'u'
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
        startTime = spikeStruct.unitStartTime; endTime = spikeStruct.unitEndTime;
        firstSpikeTimes = spikeStruct.tsTask(1);

        targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
        targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut >= startTime & targetDimTrialAttOut <= endTime);
        targetDimTrialAttOut = targetDimTrialAttOut - firstSpikeTimes;
        targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
        targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn >= startTime & targetDimTrialAttIn <= endTime);
        targetDimTrialAttIn = targetDimTrialAttIn - firstSpikeTimes;
        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut >= startTime & arrayOnsetTrialAttOut <= endTime);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut - firstSpikeTimes;
        arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn >= startTime & arrayOnsetTrialAttIn <= endTime);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn - firstSpikeTimes;
        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
        cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut >= startTime & cueOnsetTrialAttOut <= endTime);
        cueOnsetTrialAttOut = cueOnsetTrialAttOut - firstSpikeTimes;
        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
        cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn >= startTime & cueOnsetTrialAttIn <= endTime);
        cueOnsetTrialAttIn = cueOnsetTrialAttIn - firstSpikeTimes;

        targetDimTrialAttOut = targetDimTrialAttOut(cueOnsetTrialAttOut>0);
        targetDimTrialAttIn = targetDimTrialAttIn(cueOnsetTrialAttIn>0);
        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(cueOnsetTrialAttOut>0);
        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(cueOnsetTrialAttIn>0);
        cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut>0);
        cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn>0);

        spikeTimes2use = spikeStruct.tsTask;
        
        % compute FR count
        for cuei = 1:numel(UE.cueOnsetByLoc{3})
            nwinr = UE.cueOnsetByLoc{3}(cuei);
            nwinl = UE.cueOnsetByLoc{3}(cuei) - 0.175;
            frCount.cueInRFBaseline(cuei) = length(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr));
        end
        for cuei = 1:numel(UE.cueOnsetByLoc{3})
            nwinr = UE.cueOnsetByLoc{3}(cuei) + 0.200;
            nwinl = UE.cueOnsetByLoc{3}(cuei) + 0.025;
            frCount.cueInRF(cuei) = length(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr));
        end
        
        for cuei = 1:numel(UE.cueOnsetByLoc{4})
            nwinr = UE.cueOnsetByLoc{4}(cuei);
            nwinl = UE.cueOnsetByLoc{4}(cuei) - 0.175;
            frCount.cueOutRFBaseline(cuei) = length(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr));
        end
        for cuei = 1:numel(UE.cueOnsetByLoc{4})
            nwinr = UE.cueOnsetByLoc{4}(cuei) + 0.200;
            nwinl = UE.cueOnsetByLoc{4}(cuei) + 0.025;
            frCount.cueOutRF(cuei) = length(spikeTimes2use(spikeTimes2use>nwinl & spikeTimes2use<nwinr));
        end
        
        medianFRAttInCount(uniti) = mean(frCount.cueInRF);
        medianFRAttInBlCount(uniti) = mean(frCount.cueInRFBaseline);
        medianFRAttOutCount(uniti) = mean(frCount.cueOutRF);
        medianFRAttOutBlCount(uniti) = mean(frCount.cueOutRFBaseline);
        
        [pCI(uniti)] = signrank(frCount.cueInRFBaseline,frCount.cueInRF);
        [pCO(uniti)] = signrank(frCount.cueOutRFBaseline,frCount.cueOutRF);
                
        [firingRate(uniti), percWDTwoSpikeBurst(uniti), percWDThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLong(uniti),...
            numBurstShort(uniti), frAttInPerTrial,frAttOutPerTrial,...
            percWDTwoSpikeBurstAttIn(uniti), percWDThreeSpikeBurstAttIn(uniti),...
            numInterSpikeTimesAttIn(uniti), numBurstLongAttIn(uniti), numBurstShortAttIn(uniti),...
            percWDTwoSpikeBurstAttOut(uniti), percWDThreeSpikeBurstAttOut(uniti),...
            numInterSpikeTimesAttOut(uniti), numBurstLongAttOut(uniti), numBurstShortAttOut(uniti)] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut);

        meanFrAttInPerTrial(uniti) = mean(frAttInPerTrial);
        meanFrAttOutPerTrial(uniti) = mean(frAttOutPerTrial);
    else
        continue
    end
    clear nSpikesAttInPerTrial nSpikesAttOutPerTrial
    fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
end

figure
subplot(211)
histogram(percWDThreeSpikeBurst,linspace(0,1,110))
title('Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percWDThreeSpikeBurstAttOut,linspace(0,1,110))
hold on
histogram(percWDThreeSpikeBurstAttIn,linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percWDThreeSpikeBurst(percWDThreeSpikeBurst>0),linspace(0,1,110))
title('Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percWDThreeSpikeBurstAttOut(percWDThreeSpikeBurstAttOut>0),linspace(0,1,110))
hold on
histogram(percWDThreeSpikeBurstAttIn(percWDThreeSpikeBurstAttIn>0),linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percWDTwoSpikeBurst,linspace(0,5,110))
title('Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percWDTwoSpikeBurstAttOut,linspace(0,5,110))
hold on
histogram(percWDTwoSpikeBurstAttIn,linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percWDTwoSpikeBurst(percWDTwoSpikeBurst>0),linspace(0,5,110))
title('Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percWDTwoSpikeBurstAttOut(percWDTwoSpikeBurstAttOut>0),linspace(0,5,110))
hold on
histogram(percWDTwoSpikeBurstAttIn(percWDTwoSpikeBurstAttIn>0),linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})


figure
subplot(211)
histogram(firingRate,linspace(0,30,220))
title('Firing Rate')
ylabel('Count')
subplot(212)
histogram(meanFrAttOutPerTrial,linspace(0,30,220))
hold on
histogram(meanFrAttInPerTrial,linspace(0,30,220))
ylabel('Count'); xlabel('Frequency in Hz')
legend({'AttOut' 'AttIn'})

% visually responsive for cue location with higher FR post cue compared to
% baseline
visResp = (pCI <= .05 & medianFRAttInCount > medianFRAttInBlCount) | (pCO <= .05 & medianFRAttOutCount > medianFRAttOutBlCount);

figure
subplot(211)
histogram(percWDThreeSpikeBurst(visResp),linspace(0,1,110))
title('VisResp Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percWDThreeSpikeBurstAttOut(visResp),linspace(0,1,110))
hold on
histogram(percWDThreeSpikeBurstAttIn(visResp),linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percWDThreeSpikeBurst(percWDThreeSpikeBurst>0 & visResp),linspace(0,1,110))
title('VisResp Classic Bursts that contain at least 3 spikes')
ylabel('Count')
subplot(212)
histogram(percWDThreeSpikeBurstAttOut(percWDThreeSpikeBurstAttOut>0 & visResp),linspace(0,1,110))
hold on
histogram(percWDThreeSpikeBurstAttIn(percWDThreeSpikeBurstAttIn>0 & visResp),linspace(0,1,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percWDTwoSpikeBurst(visResp),linspace(0,5,110))
title('VisResp Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percWDTwoSpikeBurstAttOut(visResp),linspace(0,5,110))
hold on
histogram(percWDTwoSpikeBurstAttIn(visResp),linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})

figure
subplot(211)
histogram(percWDTwoSpikeBurst(percWDTwoSpikeBurst>0 & visResp),linspace(0,5,110))
title('VisResp Classic Bursts that contain at least 2 spikes')
ylabel('Count')
subplot(212)
histogram(percWDTwoSpikeBurstAttOut(percWDTwoSpikeBurstAttOut>0 & visResp),linspace(0,5,110))
hold on
histogram(percWDTwoSpikeBurstAttIn(percWDTwoSpikeBurstAttIn>0 & visResp),linspace(0,5,110))
ylabel('Count')
xlabel('Percentage in %')
legend({'AttOut' 'AttIn'})


figure
subplot(211)
histogram(firingRate(visResp),linspace(0,30,220))
title('VisResp Firing Rate')
ylabel('Count')
subplot(212)
histogram(meanFrAttOutPerTrial(visResp),linspace(0,30,220))
hold on
histogram(meanFrAttInPerTrial(visResp),linspace(0,30,220))
ylabel('Count'); xlabel('Frequency in Hz')
legend({'AttOut' 'AttIn'})


figure
histogram(firingRate(percWDTwoSpikeBurst>=1),linspace(1,50,100))
hold on
histogram(firingRate(percWDTwoSpikeBurst<1),linspace(1,50,100))
legend({'Burst >= 1%' 'Burst<1%'})
title('Firing rate all units')
xlabel('Firing rate in Hz')
ylabel('Count')

figure
histogram(firingRate(visResp & percWDTwoSpikeBurst>=1),linspace(1,50,100))
hold on
histogram(firingRate(visResp & percWDTwoSpikeBurst<1),linspace(1,50,100))
legend({'Burst >= 1%' 'Burst<1%'})
title('Firing rate visResp units')
xlabel('Firing rate in Hz')
ylabel('Count')




% eof