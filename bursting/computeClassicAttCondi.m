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

meanFRAttInCount = nan(1,nUnits);
meanFRAttInBlCount = nan(1,nUnits);
meanFRAttOutCount = nan(1,nUnits);
meanFRAttOutBlCount = nan(1,nUnits);
pCI = nan(1,nUnits);
pCO = nan(1,nUnits);

for uniti = 1:nUnits
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
        
        meanFRAttInCount(uniti) = mean(frCount.cueInRF);
        meanFRAttInBlCount(uniti) = mean(frCount.cueInRFBaseline);
        meanFRAttOutCount(uniti) = mean(frCount.cueOutRF);
        meanFRAttOutBlCount(uniti) = mean(frCount.cueOutRFBaseline);
        
        [pCI(uniti)] = signrank(frCount.cueInRFBaseline,frCount.cueInRF);
        [pCO(uniti)] = signrank(frCount.cueOutRFBaseline,frCount.cueOutRF);
        
        tempResolve = 1;
        
%         [firingRate(uniti), percClassicTwoSpikeBurst(uniti), percClassicThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLong(uniti),...
%             numBurstShort(uniti), frAttInPerTrial,frAttOutPerTrial,...
%             percClassicTwoSpikeBurstAttIn(uniti), percClassicThreeSpikeBurstAttIn(uniti),...
%             numInterSpikeTimesAttIn(uniti), numBurstLongAttIn(uniti), numBurstShortAttIn(uniti),...
%             percClassicTwoSpikeBurstAttOut(uniti), percClassicThreeSpikeBurstAttOut(uniti),...
%             numInterSpikeTimesAttOut(uniti), numBurstLongAttOut(uniti), numBurstShortAttOut(uniti)] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
%                         arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, tempResolve);
        [firingRate(uniti), percClassicTwoSpikeBurst(uniti), percClassicThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLong(uniti),...
            numBurstShort(uniti), frAttInPerTrial,frAttOutPerTrial,...
            percClassicTwoSpikeBurstAttIn(uniti),percClassicTwoSpikeBurstAttOut(uniti)] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, tempResolve);

        meanFrAttInPerTrial(uniti) = mean(frAttInPerTrial);
        meanFrAttOutPerTrial(uniti) = mean(frAttOutPerTrial);
%     else
%         continue
%     end
    clear nSpikesAttInPerTrial nSpikesAttOutPerTrial
    fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
end

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
visResp = (pCI <= .05 & meanFRAttInCount > meanFRAttInBlCount) | (pCO <= .05 & meanFRAttOutCount > meanFRAttOutBlCount);

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

meanFRAttInCount = nan(1,nUnits);
meanFRAttInBlCount = nan(1,nUnits);
meanFRAttOutCount = nan(1,nUnits);
meanFRAttOutBlCount = nan(1,nUnits);
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
        
        meanFRAttInCount(uniti) = mean(frCount.cueInRF);
        meanFRAttInBlCount(uniti) = mean(frCount.cueInRFBaseline);
        meanFRAttOutCount(uniti) = mean(frCount.cueOutRF);
        meanFRAttOutBlCount(uniti) = mean(frCount.cueOutRFBaseline);
        
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
visResp = (pCI <= .05 & meanFRAttInCount > meanFRAttInBlCount) | (pCO <= .05 & meanFRAttOutCount > meanFRAttOutBlCount);

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