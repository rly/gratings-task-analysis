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
discrep = zeros(1,nUnits);
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
        [firingRate(uniti), percClassicTwoSpikeBurst(uniti), percClassicThreeSpikeBurst(uniti), numInterSpikeTimes(uniti), numBurstLong(uniti),...
            numBurstShort(uniti), discrep(uniti), frAttInPerTrial,frAttOutPerTrial,...
            percClassicTwoSpikeBurstAttIn(uniti), percClassicThreeSpikeBurstAttIn(uniti),...
            numInterSpikeTimesAttIn(uniti), numBurstLongAttIn(uniti), numBurstShortAttIn(uniti),...
            percClassicTwoSpikeBurstAttOut(uniti), percClassicThreeSpikeBurstAttOut(uniti),...
            numInterSpikeTimesAttOut(uniti), numBurstLongAttOut(uniti), numBurstShortAttOut(uniti)] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
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



% eof