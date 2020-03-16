function [epochedSpikeTimesAttIn, epochedSpikeTimesAttOut,...
    frAttInPerTrial, frAttOutPerTrial,frBaselinePerTrial, poissonBinaryAttIn,...
    poissonBinaryAttOut,datatmpShuffledAttIn,datatmpShuffledAttOut,...
    spikesInTrialsAttIn,spikesInTrialsAttOut] = createEpochedSpikeTimes(allSpikeTimes,...
    endTimeAttIn, endTimeAttOut, startTimeAttOut,startTimeAttIn, startTime, endTime,...
    firstSpikeTimes,SCM,Poisson,Shuffled,nTrialsAttIn)

% Cuts the spikeTimes in trial epochs based on start and end times of trial
% events (e.g. cue onset and array onset or target dim) and spits out both 
% spike times and binary data cut into epochs.
% 
% If data is already epoched it can also cut out the trial times of
% interest, e.g. 200ms post cue onset.
%
% With SCM input set to 1 it will also perform spike count matching across
% trials (i.e. ditching a random trial until the spike counts between
% different attention conditions are matched).



% For debugging:
% endTimeAttIn = arrayOnsetTrialAttIn;
% endTimeAttOut = arrayOnsetTrialAttOut;
% startTimeAttOut = cueOnsetTrialAttOut;
% startTimeAttIn = cueOnsetTrialAttIn;
% allSpikeTimes = spikeTimes2use;
% allSpikeTimes = [SpikeLFPInfo_all.CueOnSpike{1,1}{1}; SpikeLFPInfo_all.CueOffSpike{1,1}{1}];
% nTrialsAttIn = size(SpikeLFPInfo_all.CueOnSpike{1,1}{1},1);
% startTimeAttIn = 200;
% endTimeAttIn = 600;
% startTimeAttOut = 200;
% endTimeAttOut = 600;
% [epochedSpikeTimesAttIn, epochedSpikeTimesAttOut,...
%     nSpikesAttInPerTrial, nSpikesAttOutPerTrial,poissonBinaryAttIn,...
%     poissonBinaryAttOut,datatmpShuffledAttIn,datatmpShuffledAttOut] = createEpochedSpikeTimes([SpikeLFPInfo_all.CueOnSpike{loci,sessioni}{uniti}; SpikeLFPInfo_all.CueOffSpike{loci,sessioni}{uniti}],...
%     200, 200, 200,200, 0, 0,...
%     0,SCM,Poisson,Shuffled,nTrialsAttIn)

Fs = 1000;

if size(allSpikeTimes,1) == 1
    startTimeAttOut = startTimeAttOut(startTimeAttOut >= startTime & startTimeAttOut <= endTime);
    startTimeAttOut = startTimeAttOut - firstSpikeTimes; startTimeAttOuttmp = startTimeAttOut;
    startTimeAttOut = startTimeAttOut(startTimeAttOuttmp>0);

    startTimeAttIn = startTimeAttIn(startTimeAttIn >= startTime & startTimeAttIn <= endTime);
    startTimeAttIn = startTimeAttIn - firstSpikeTimes; startTimeAttIntmp = startTimeAttIn;
    startTimeAttIn = startTimeAttIn(startTimeAttIntmp>0);

    endTimeAttOut = endTimeAttOut(endTimeAttOut >= startTime & endTimeAttOut <= endTime);
    endTimeAttOut = endTimeAttOut - firstSpikeTimes;
    endTimeAttOut = endTimeAttOut(startTimeAttOuttmp>0);

    endTimeAttIn = endTimeAttIn(endTimeAttIn >= startTime & endTimeAttIn <= endTime);
    endTimeAttIn = endTimeAttIn - firstSpikeTimes;
    endTimeAttIn = endTimeAttIn(startTimeAttIntmp>0);
    
    data_pts = round(allSpikeTimes(1),3):1/1000:round(allSpikeTimes(end),3);
    binarySpikeTrain = zeros(1,length(data_pts));
    binarySpikeTrain(ismembertol(data_pts,round(allSpikeTimes,3),.00000001)) = 1;
else
    if size(allSpikeTimes,2) == 1601
        timeOnTrial = -600:1:1000;
        [~,timeidx(1)] = min(abs(timeOnTrial - startTimeAttIn));
        [~,timeidx(2)] = min(abs(timeOnTrial - endTimeAttIn));
    end
    binarySpikeTrain = [];
    for n = 1:size(allSpikeTimes,1)
        binarySpikeTrain = [binarySpikeTrain allSpikeTimes(n,:)];
    end
    allSpikeTimesNonBin = find(binarySpikeTrain);
    for n = 1:nTrialsAttIn
        startTimeAttIntmp(n) = timeidx(1) + ((n-1)*size(allSpikeTimes,2)); 
        endTimeAttIntmp(n) = timeidx(2) + ((n-1)*size(allSpikeTimes,2));
    end
    for n = nTrialsAttIn:size(allSpikeTimes,1)
        startTimeAttOuttmp(n) = timeidx(1) + ((n-1)*size(allSpikeTimes,2)); 
        endTimeAttOuttmp(n) = timeidx(2) + ((n-1)*size(allSpikeTimes,2));
    end
    startTimeAttIn = startTimeAttIntmp / Fs; endTimeAttIn = endTimeAttIntmp / Fs;
    startTimeAttOut = startTimeAttOuttmp(startTimeAttOuttmp>0)  / Fs; 
    endTimeAttOut = endTimeAttOuttmp(endTimeAttOuttmp>0) / Fs;
    allSpikeTimesBin = allSpikeTimes;
    allSpikeTimes = []; allSpikeTimes = allSpikeTimesNonBin / Fs;
end
    
% calculate trial based spiking
data = binarySpikeTrain;
NEAttIn=length(startTimeAttIn); NEAttOut=length(startTimeAttOut);
nEAttIn=floor(startTimeAttIn*Fs)+1; nEAttOut=floor(startTimeAttOut*Fs)+1;
for n=1:NEAttIn
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttIn(n)*Fs);
    indx=nEAttIn(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        nSpikesAttInPerTrial(n) = sum(data(indx));
        frAttInPerTrial(n) = sum(data(indx)) / (length(data(indx)) / Fs);
        epochedSpikeTimesAttIn(n).spikeTimes=find(data(indx));
        epochedSpikeTimesAttIn(n).binData=data(indx);
        clear indx nwinr
    end
    if nEAttIn(n) > 500
        indx=(nEAttIn(n)-500):(nEAttIn(n)-201); % baseline time window
        if length(indx) >1 && indx(end) < size(data,2)
            frAttInBaselinePerTrial(n) = sum(data(indx)) / (length(data(indx)) / Fs);
            clear indx 
        end
    end
end
startTimeAttIn = startTimeAttIn(1:length(nSpikesAttInPerTrial)); NEAttIn=length(startTimeAttIn);
endTimeAttIn = endTimeAttIn(1:length(nSpikesAttInPerTrial));
for n=1:NEAttOut;
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        frAttOutPerTrial(n) = sum(data(indx)) / (length(data(indx)) / Fs);
        epochedSpikeTimesAttOut(n).spikeTimes=find(data(indx));
        epochedSpikeTimesAttOut(n).binData=data(indx);
        clear indx
    end
    if nEAttOut(n) > 300
        indx=(nEAttOut(n)-300):(nEAttOut(n)-1); % baseline time window
        if length(indx) >1 && indx(end) < size(data,2)
            frAttOutBaselinePerTrial(n) = sum(data(indx)) / (length(data(indx)) / Fs);
            clear indx 
        end
    end
end
startTimeAttOut = startTimeAttOut(1:length(nSpikesAttOutPerTrial)); NEAttOut=length(startTimeAttOut);
endTimeAttOut = endTimeAttOut(1:length(nSpikesAttOutPerTrial));
frBaselinePerTrial = [frAttOutBaselinePerTrial frAttInBaselinePerTrial];

if SCM
    % Perform spike count matching
    nSpikesAttIn = sum(nSpikesAttInPerTrial);
    nSpikesAttOut = sum(nSpikesAttOutPerTrial);
    counter = 0;
    if nSpikesAttIn > nSpikesAttOut
        trialSel = randperm(NEAttIn);
        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
            clear nSpikesAttInPerTrialNw epochedSpikeTimesAttInNw
            counter = counter + 1;
            NEAttInNw = NEAttIn - counter;
            startTimeAttInNw = startTimeAttIn(trialSel(1:NEAttInNw)); endTimeAttInNw = endTimeAttIn(trialSel(1:NEAttInNw));
            for n=1:NEAttInNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(endTimeAttInNw(n)*Fs);
                indx=round(startTimeAttInNw(n)*Fs):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    nSpikesAttInPerTrialNw(n) = sum(data(indx));
                    epochedSpikeTimesAttInNw(n).spikeTimes=find(data(indx));
                    epochedSpikeTimesAttInNw(n).binData=data(indx);
                end
                clear indx
            end
            if exist('nSpikesAttInPerTrialNw','var')
                nSpikesAttInPerTrialNw2use = nSpikesAttInPerTrialNw;
                nSpikesAttIn = sum(nSpikesAttInPerTrialNw2use);
                clear nSpikesAttInPerTrialNw2use 
            else
                break
            end
        end
        epochedSpikeTimesAttIn = epochedSpikeTimesAttInNw;
        nSpikesAttInPerTrial = nSpikesAttInPerTrialNw;
        NEAttOutNw = NEAttOut; endTimeAttOutNw = endTimeAttOut; startTimeAttOutNw = startTimeAttOut;
    elseif nSpikesAttOut > nSpikesAttIn
        trialSel = randperm(NEAttOut);
        while nSpikesAttOut + mean(nSpikesAttOutPerTrial) / 2 > nSpikesAttIn
            clear nSpikesAttOutPerTrialNw epochedSpikeTimesAttOutNw
            counter = counter + 1;
            NEAttOutNw = NEAttOut - counter;
            startTimeAttOutNw = startTimeAttOut(trialSel(1:NEAttOutNw)); endTimeAttOutNw = endTimeAttOut(trialSel(1:NEAttOutNw));
            for n=1:NEAttOutNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(endTimeAttOutNw(n)*Fs);
                indx=round(startTimeAttOutNw(n)*Fs):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    nSpikesAttOutPerTrialNw(n) = sum(data(indx));
                    epochedSpikeTimesAttOutNw(n).spikeTimes=find(data(indx));
                    epochedSpikeTimesAttOutNw(n).binData=data(indx);
                end
                clear indx
            end
            nSpikesAttOutPerTrialNw2use = nSpikesAttOutPerTrialNw;
            nSpikesAttOut = sum(nSpikesAttOutPerTrialNw2use);
            clear nSpikesAttOutPerTrialNw2use 
        end
        epochedSpikeTimesAttOut = epochedSpikeTimesAttOutNw;
        nSpikesAttOutPerTrial = nSpikesAttOutPerTrialNw;
        NEAttInNw = NEAttIn; endTimeAttInNw = endTimeAttIn; startTimeAttInNw = startTimeAttIn;
    elseif nSpikesAttIn == nSpikesAttOut
        NEAttOutNw = NEAttOut; endTimeAttOutNw = endTimeAttOut; startTimeAttOutNw = startTimeAttOut;
        NEAttInNw = NEAttIn; endTimeAttInNw = endTimeAttIn; startTimeAttInNw = startTimeAttIn;
    end
    
    if Poisson
        % create Poisson data for comparison
        [poissonBinaryAttIn,~] = calcPoisson(NEAttInNw,endTimeAttInNw,startTimeAttInNw,Fs,data);
        [poissonBinaryAttOut,~] = calcPoisson(NEAttOutNw,endTimeAttOutNw,startTimeAttOutNw,Fs,data);
    end
    
    if Shuffled
        spikeTimesAllTrials = [];
        for n = 1:NEAttOutNw
            spikeTimesAllTrials = [spikeTimesAllTrials epochedSpikeTimesAttOut(n).spikeTimes];
        end
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for n = 1:NEAttOutNw
            if n == 1
                if nSpikesAttOutPerTrial(n) ~= 0
                    nspike2use = spikeTimesAllTrials2use(1:nSpikesAttOutPerTrial(1));
                    idx = nSpikesAttOutPerTrial(1) + 1;
                else
                    idx = 1;
                    continue
                end
            else
                if nSpikesAttOutPerTrial(n) ~= 0
                    if idx+nSpikesAttOutPerTrial<length(spikeTimesAllTrials2use)
                        nspike2use = spikeTimesAllTrials2use(idx:idx+nSpikesAttOutPerTrial(n)-1);
                    else
                        nspike2use = spikeTimesAllTrials2use(idx:end);
                    end
                    idx = nSpikesAttOutPerTrial(n) + idx;
                else
                    attOutRandomized(n).spikeTimes = [];
                    continue
                end
            end
            attOutRandomized(n).spikeTimes = sort(nspike2use);
            clear nspike2use
        end
        % finally get the interspike intervals of the shuffled data
        datatmpShuffledAttOut = [];
        for n = 1:length(attOutRandomized)
            datatmpShuffledAttOut = [datatmpShuffledAttOut diff(attOutRandomized(n).spikeTimes)];
        end

        clear spikeTimesAllTrials2use nSpikesPerTrial
        spikeTimesAllTrials = [];
        for n = 1:NEAttInNw
            spikeTimesAllTrials = [spikeTimesAllTrials epochedSpikeTimesAttIn(n).spikeTimes];
        end
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for n = 1:NEAttInNw
            if n == 1
                if nSpikesAttInPerTrial(n) ~= 0
                    nspike2use = spikeTimesAllTrials2use(1:nSpikesAttInPerTrial(1));
                    idx = nSpikesAttInPerTrial(1) + 1;
                else
                    idx = 1;
                    continue
                end
            else
                if nSpikesAttInPerTrial(n) ~= 0
                    if idx+nSpikesAttInPerTrial<length(spikeTimesAllTrials2use)
                        nspike2use = spikeTimesAllTrials2use(idx:idx+nSpikesAttInPerTrial(n)-1);
                    else
                        nspike2use = spikeTimesAllTrials2use(idx:end);
                    end
                    idx = nSpikesAttInPerTrial(n) + idx;
                else
                    attInRandomized(n).spikeTimes = [];
                    continue
                end
            end
            attInRandomized(n).spikeTimes = sort(nspike2use);
            clear nspike2use
        end
        % finally get the interspike intervals of the shuffled data
        datatmpShuffledAttIn = [];
        for n = 1:length(attInRandomized)
            datatmpShuffledAttIn = [datatmpShuffledAttIn diff(attInRandomized(n).spikeTimes)];
        end
    end
else
    datatmpShuffledAttIn = 0; datatmpShuffledAttOut = 0;
    poissonBinaryAttIn = 0; poissonBinaryAttOut = 0;
end

spikeFs = 1000;
for j = 1:length(startTimeAttIn)
    spikeTimesInTrial = allSpikeTimes(allSpikeTimes >= startTimeAttIn(j) + 1/spikeFs & ...
            allSpikeTimes <= endTimeAttIn(j));
    spikeTimesInTrialAdj{j} = round((spikeTimesInTrial - startTimeAttIn(j)) * spikeFs);
    spikesInTrialsAttIn(j, spikeTimesInTrialAdj{j}) = 1;
end
for j = 1:length(startTimeAttOut)
    spikeTimesInTrial = allSpikeTimes(allSpikeTimes >= startTimeAttOut(j) + 1/spikeFs & ...
            allSpikeTimes <= endTimeAttOut(j));
    spikeTimesInTrialAdj{j} = round((spikeTimesInTrial - startTimeAttOut(j)) * spikeFs);
    spikesInTrialsAttOut(j, spikeTimesInTrialAdj{j}) = 1;
end