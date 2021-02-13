function [firingRatefun, percWDTwoSpikeBurstfun, percWDThreeSpikeBurstfun, numInterSpikeTimesfun, numBurstLongfun,...
    numBurstShortfun, frAttInPerTrialfun,frAttOutPerTrialfun,...
    percWDTwoSpikeBurstAttInfun, percWDTwoSpikeBurstAttOutfun,...
    percWDTwoSpikeBurstBaselinefun, percWDTwoSpikeBurstAttInBaselinefun, percWDTwoSpikeBurstAttOutBaselinefun,...
    percWDTwoSpikeBurstPostCuefun, percWDTwoSpikeBurstAttInPostCuefun, percWDTwoSpikeBurstAttOutPostCuefun,...
    percWDTwoSpikeBurstPostArrayfun, percWDTwoSpikeBurstAttInPostArrayfun, percWDTwoSpikeBurstAttOutPostArrayfun,...
    cueOnset,arrayOnset,targetDimOnset,percRyanTwoSpikeBurstAttInfun,percRyanTwoSpikeBurstAttOutfun] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut, targetDimTrialAttIn,targetDimTrialAttOut, tempResolve, minTime, UE, startTime, endTime)

      
% [firingRatefun, percClassicTwoSpikeBurstfun, percClassicThreeSpikeBurstfun, numInterSpikeTimesfun, numBurstLongfun,...
%     numBurstShortfun, frAttInPerTrialfun,frAttOutPerTrialfun,...
%     percClassicTwoSpikeBurstAttInfun, percClassicThreeSpikeBurstAttInfun,...
%     numInterSpikeTimesAttInfun, numBurstLongAttInfun, numBurstShortAttInfun,...
%     percClassicTwoSpikeBurstAttOutfun, percClassicThreeSpikeBurstAttOutfun,...
%     numInterSpikeTimesAttOutfun, numBurstLongAttOutfun, numBurstShortAttOutfun] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
%                 endTimeAttIn, endTimeAttOut, tempResolve)
%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;
%maxTime = preQuietPeriod; minTime = withinBurstTime;

% calculate traditional bursting in all spike times
interSpikeTimes = diff(spikeTimes2use);
interSpikeTimes = interSpikeTimes*1000; % in ms
firingRatefun = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));

data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

burstInd = find(interSpikeTimes<=minTime);
if length(burstInd) > 1
    if sum(diff(burstInd)==1) ~= 0
        burstInd(find(diff(burstInd)==1)+1) = [];
    end
end
percWDTwoSpikeBurstfun = numel(interSpikeTimes(burstInd)) * 100 / numel(spikeTimes2use);

kernelSigma = .01;
if tempResolve
    cueOnset.window = [0.8 0.8]; % seconds before, after
    cueOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
    cueOnset = createTimeLockedSpdf(spikeTimes2use(burstInd+1), UE.cueOnset, UE.cueOnsetByLoc, cueOnset, kernelSigma, startTime, endTime);
    arrayOnset.window = [0.8 0.8]; % seconds before, after
    arrayOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
    arrayOnset = createTimeLockedSpdf(spikeTimes2use(burstInd+1), UE.arrayOnset, UE.arrayOnsetByLoc, arrayOnset, kernelSigma, startTime, endTime);
    targetDimOnset.window = [0.8 0.8]; % seconds before, after
    targetDimOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
    targetDimOnset = createTimeLockedSpdf(spikeTimes2use(burstInd+1), UE.targetDim, UE.targetDimBalByLoc, targetDimOnset, kernelSigma, startTime, endTime);    
end

if ~isempty(burstInd)
    if burstInd(end)+1 >= numel(interSpikeTimes); burstInd = burstInd(1:end-1); end
    burstIndLong = find(interSpikeTimes(burstInd+1)<=minTime);
    if length(burstIndLong) > 1
        if sum(diff(burstIndLong)==1) ~= 0
            burstIndLong(find(diff(burstIndLong)==1)+1) = [];
        end
    end
    percWDThreeSpikeBurstfun = numel(burstIndLong) * 100 / numel(spikeTimes2use);
    numInterSpikeTimesfun = numel(interSpikeTimes);
    numBurstLongfun = numel(burstIndLong);
    numBurstShortfun = numel(burstInd);
else
    percWDThreeSpikeBurstfun = 0;
    numInterSpikeTimesfun = numel(interSpikeTimes);
    numBurstLongfun = 0;
    numBurstShortfun = 0;
end

% % binarySpikeTrain = zeros(1,length(data_pts));
% % binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,4),.0000000582)) = 1;
% % sum(binarySpikeTrain) 
% % numel(spikeTimes2use)
% 
% if numel(spikeTimes2use)==sum(binarySpikeTrain)
%     discrep = -1;
% else
%     discrep = 1;
% end

% calculate trial based spiking
% 200ms post-cue until TARGET DIM
data = binarySpikeTrain;
Fs = 1000;
NEAttIn=length(cueOnsetTrialAttIn); NEAttOut=length(cueOnsetTrialAttOut);
nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1; nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
datatmpAttIn=[];datatmpAttOut=[];
datatmpAttInBaseline = []; datatmpAttOutBaseline = [];
datatmpAttInPostCue = []; datatmpAttOutPostCue = [];
datatmpAttInPostArray = []; datatmpAttOutPostArray = [];
isiPerTrialAttOutBaseline = []; isiPerTrialAttOutPostCue = []; isiPerTrialAttOutPostArray = [];

isiPerTrialAttIn = struct(); isiPerTrialAttInBaseline = struct(); isiPerTrialAttInPostCue = struct(); isiPerTrialAttInPostArray = struct();
nSpikesAttInPerTrial = nan(NEAttIn,1); nSpikesAttInBaseline = nan(NEAttIn,1); nSpikesAttInPostCue = nan(NEAttIn,1); nSpikesAttInPostArray = nan(NEAttIn,1);
frAttInPerTrialfun = nan(NEAttIn,1); frAttInPerTrialBaselinefun = nan(NEAttIn,1);frAttInPerTrialPostCuefun = nan(NEAttIn,1);frAttInPerTrialPostArrayfun = nan(NEAttIn,1);
isiPerTrialAttOut = struct(); isiPerTrialAttOutBaseline = struct(); isiPerTrialAttOutPostCue = struct(); isiPerTrialAttOutPostArray = struct();
nSpikesAttOutPerTrial = nan(NEAttOut,1); nSpikesAttOutBaseline = nan(NEAttOut,1); nSpikesAttOutPostCue = nan(NEAttOut,1); nSpikesAttOutPostArray = nan(NEAttOut,1);
frAttOutPerTrialfun = nan(NEAttOut,1); frAttOutPerTrialBaselinefun = nan(NEAttOut,1);frAttOutPerTrialPostCuefun = nan(NEAttOut,1);frAttOutPerTrialPostArrayfun = nan(NEAttOut,1);

for n=1:NEAttIn
    nwinr=round(endTimeAttIn(n)*Fs); % arrayOnset
    indx=nEAttIn(n):nwinr-1; % cueOnset:arrayOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
        isiPerTrialAttIn(n).isi = diff(find(data(indx)));
        nSpikesAttInPerTrial(n) = sum(data(indx));
        frAttInPerTrialfun(n) = nSpikesAttInPerTrial(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
    
    % lock to cue to get baseline
    nwinr=round(nEAttIn(n)); % cueOnset
    nwinl=round(nEAttIn(n) - (0.350*Fs)); % baseline period
    if nwinl<1;nwinl=1;end
    indx=nwinl:nwinr-1; % baseline:cueOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttInBaseline=[datatmpAttInBaseline diff(find(data(indx)))];
        isiPerTrialAttInBaseline(n).isi = diff(find(data(indx)));
        nSpikesAttInBaseline(n) = sum(data(indx));
        frAttInPerTrialBaselinefun(n) = nSpikesAttInBaseline(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
    
    % lock to arrayonset to get cue period
    nwinr=round(endTimeAttIn(n)*Fs); % arrayOnset
    nwinl=round((endTimeAttIn(n)-0.350)*Fs); % cue period
    indx=nwinl:nwinr-1; % cueOnset:postCueOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttInPostCue=[datatmpAttInPostCue diff(find(data(indx)))];
        isiPerTrialAttInPostCue(n).isi = diff(find(data(indx)));
        nSpikesAttInPostCue(n) = sum(data(indx));
        frAttInPerTrialPostCuefun(n) = nSpikesAttInPostCue(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
    
    % lock to target dim to get array period
    nwinr=round(targetDimTrialAttIn(n)*Fs); % target dimmin
    nwinl=round((targetDimTrialAttIn(n)-0.350)*Fs); % array period
    indx=nwinl:nwinr-1; % arrayOnset:postArrayOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttInPostArray=[datatmpAttInPostArray diff(find(data(indx)))];
        isiPerTrialAttInPostArray(n).isi = diff(find(data(indx)));
        nSpikesAttInPostArray(n) = sum(data(indx));
        frAttInPerTrialPostArrayfun(n) = nSpikesAttInPostArray(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
end


for n=1:NEAttOut;
    nwinr=round(endTimeAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
        isiPerTrialAttOut(n).isi = diff(find(data(indx)));
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        frAttOutPerTrialfun(n) = nSpikesAttOutPerTrial(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
    
     % lock to cue to get baseline
    nwinr=round(nEAttOut(n)); % cueOnset
    nwinl=round(nEAttOut(n) - (0.350*Fs)); % baseline period
    if nwinl<1;nwinl=1;end
    indx=nwinl:nwinr-1; % baseline:cueOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOutBaseline=[datatmpAttOutBaseline diff(find(data(indx)))];
        isiPerTrialAttOutBaseline(n).isi = diff(find(data(indx)));
        nSpikesAttOutBaseline(n) = sum(data(indx));
        frAttOutPerTrialBaselinefun(n) = nSpikesAttOutBaseline(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
    
    % lock to arrayonset to get cue period
    nwinr=round(endTimeAttOut(n)*Fs); % arrayOnset
    nwinl=round((endTimeAttOut(n)-0.350)*Fs); % cue period
    indx=nwinl:nwinr-1; % cueOnset:postCueOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOutPostCue=[datatmpAttOutPostCue diff(find(data(indx)))];
        isiPerTrialAttOutPostCue(n).isi = diff(find(data(indx)));
        nSpikesAttOutPostCue(n) = sum(data(indx));
        frAttOutPerTrialPostCuefun(n) = nSpikesAttOutPostCue(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
    
    % lock to target dim to get array period
    nwinr=round(targetDimTrialAttOut(n)*Fs); % target dimming
    nwinl=round((targetDimTrialAttOut(n)-0.350)*Fs); % array period
    indx=nwinl:nwinr-1; % arrayOnset:postArrayOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOutPostArray=[datatmpAttOutPostArray diff(find(data(indx)))];
        isiPerTrialAttOutPostArray(n).isi = diff(find(data(indx)));
        nSpikesAttOutPostArray(n) = sum(data(indx));
        frAttOutPerTrialPostArrayfun(n) = nSpikesAttOutPostArray(n) / round(((indx(end) - indx(1))/Fs),3);
        clear indx
    end
end

% for att in
% per trial
burstIndAttInPerTrial = [];
for n = 1:numel(isiPerTrialAttIn)
    burstIndAttInTmp = find(isiPerTrialAttIn(n).isi<=minTime);
    if length(burstIndAttInTmp) > 1
        if sum(diff(burstIndAttInTmp)==1) ~= 0
            burstIndAttInTmp(find(diff(burstIndAttInTmp)==1)+1) = [];
        end
    end
    burstIndAttInPerTrial = [burstIndAttInPerTrial burstIndAttInTmp];
end
if ~isempty(burstIndAttInPerTrial)
    percWDTwoSpikeBurstAttInfun = numel(burstIndAttInPerTrial) * 100 / sum(nSpikesAttInPerTrial);
else
    percWDTwoSpikeBurstAttInfun = 0;
end
    
burstIndAttInPerTrialRyan = [];
for n = 1:numel(isiPerTrialAttIn)
    burstIndAttInTmp = find(isiPerTrialAttIn(n).isi<=10);
    if length(burstIndAttInTmp) > 1
        if sum(diff(burstIndAttInTmp)==1) ~= 0
            burstIndAttInTmp(find(diff(burstIndAttInTmp)==1)+1) = [];
        end
    end
    burstIndAttInPerTrialRyan = [burstIndAttInPerTrialRyan burstIndAttInTmp];
end
if ~isempty(burstIndAttInPerTrialRyan)
    percRyanTwoSpikeBurstAttInfun = numel(burstIndAttInPerTrialRyan) * 100 / sum(nSpikesAttInPerTrial);
else
    percRyanTwoSpikeBurstAttInfun = 0;
end

if ~isempty(isiPerTrialAttInBaseline) && ~isempty(isiPerTrialAttInPostArray) && ~isempty(isiPerTrialAttInPostCue)
    % att in baseline
    burstIndAttInBaseline = [];
    for n = 1:numel(isiPerTrialAttInBaseline)
        burstIndAttInTmp = find(isiPerTrialAttInBaseline(n).isi<=minTime);
        if length(burstIndAttInTmp) > 1
            if sum(diff(burstIndAttInTmp)==1) ~= 0
                burstIndAttInTmp(find(diff(burstIndAttInTmp)==1)+1) = [];
            end
        end
        burstIndAttInBaseline = [burstIndAttInBaseline burstIndAttInTmp];
    end
    burstIndAttInPostCue = [];
    for n = 1:numel(isiPerTrialAttInPostCue)
        burstIndAttInTmp = find(isiPerTrialAttInPostCue(n).isi<=minTime);
        if length(burstIndAttInTmp) > 1
            if sum(diff(burstIndAttInTmp)==1) ~= 0
                burstIndAttInTmp(find(diff(burstIndAttInTmp)==1)+1) = [];
            end
        end
        burstIndAttInPostCue = [burstIndAttInPostCue burstIndAttInTmp];
    end
    burstIndAttInPostArray = [];
    for n = 1:numel(isiPerTrialAttInPostArray)
        burstIndAttInTmp = find(isiPerTrialAttInPostArray(n).isi<=minTime);
        if length(burstIndAttInTmp) > 1
            if sum(diff(burstIndAttInTmp)==1) ~= 0
                burstIndAttInTmp(find(diff(burstIndAttInTmp)==1)+1) = [];
            end
        end
        burstIndAttInPostArray = [burstIndAttInPostArray burstIndAttInTmp];
    end
    if ~isempty(burstIndAttInBaseline) || ~isempty(burstIndAttInPostCue) || ~isempty(burstIndAttInPostArray)
        percWDTwoSpikeBurstAttInBaselinefun = numel(burstIndAttInBaseline) * 100 / sum(nSpikesAttInBaseline);
        percWDTwoSpikeBurstAttInPostCuefun = numel(burstIndAttInPostCue) * 100 / sum(nSpikesAttInPostCue);
        percWDTwoSpikeBurstAttInPostArrayfun = numel(burstIndAttInPostArray) * 100 / sum(nSpikesAttInPostArray);
    else
        percWDTwoSpikeBurstAttInBaselinefun = 0;
        percWDTwoSpikeBurstAttInPostCuefun = 0;
        percWDTwoSpikeBurstAttInPostArrayfun = 0;
    end
    isiFoundAttIn = 1;
else
    burstIndAttInBaseline = [];
    burstIndAttInPostCue = [];
    burstIndAttInPostArray = [];
    percWDTwoSpikeBurstAttInBaselinefun = 0;
    percWDTwoSpikeBurstAttInPostCuefun = 0;
    percWDTwoSpikeBurstAttInPostArrayfun = 0;
    isiFoundAttIn = 0;
end

% for att out
burstIndAttOutPerTrial = [];
for n = 1:numel(isiPerTrialAttOut)
    burstIndAttOutTmp = find(isiPerTrialAttOut(n).isi<=minTime);
    if length(burstIndAttOutTmp) > 1
        if sum(diff(burstIndAttOutTmp)==1) ~= 0
            burstIndAttOutTmp(find(diff(burstIndAttOutTmp)==1)+1) = [];
        end
    end
    burstIndAttOutPerTrial = [burstIndAttOutPerTrial burstIndAttOutTmp];
end
if ~isempty(burstIndAttOutPerTrial)
    percWDTwoSpikeBurstAttOutfun = numel(burstIndAttOutPerTrial) * 100 / sum(nSpikesAttOutPerTrial);
else
    percWDTwoSpikeBurstAttOutfun = 0;
end

burstIndAttOutPerTrialRyan = [];
for n = 1:numel(isiPerTrialAttOut)
    burstIndAttOutTmp = find(isiPerTrialAttOut(n).isi<=10);
    if length(burstIndAttOutTmp) > 1
        if sum(diff(burstIndAttOutTmp)==1) ~= 0
            burstIndAttOutTmp(find(diff(burstIndAttOutTmp)==1)+1) = [];
        end
    end
    burstIndAttOutPerTrialRyan = [burstIndAttOutPerTrialRyan burstIndAttOutTmp];
end
if ~isempty(burstIndAttOutPerTrialRyan)
    percRyanTwoSpikeBurstAttOutfun = numel(burstIndAttOutPerTrialRyan) * 100 / sum(nSpikesAttOutPerTrial);
else
    percRyanTwoSpikeBurstAttOutfun = 0;
end

if ~isempty(isiPerTrialAttOutBaseline) && ~isempty(isiPerTrialAttOutPostArray) && ~isempty(isiPerTrialAttOutPostCue)
    % att out baseline
    burstIndAttOutBaseline = [];
    for n = 1:numel(isiPerTrialAttOutBaseline)
        burstIndAttOutTmp = find(isiPerTrialAttOutBaseline(n).isi<=minTime);
        if length(burstIndAttOutTmp) > 1
            if sum(diff(burstIndAttOutTmp)==1) ~= 0
                burstIndAttOutTmp(find(diff(burstIndAttOutTmp)==1)+1) = [];
            end
        end
        burstIndAttOutBaseline = [burstIndAttOutBaseline burstIndAttOutTmp];
    end
    burstIndAttOutPostCue = [];
    for n = 1:numel(isiPerTrialAttOutPostCue)
        burstIndAttOutTmp = find(isiPerTrialAttOutPostCue(n).isi<=minTime);
        if length(burstIndAttOutTmp) > 1
            if sum(diff(burstIndAttOutTmp)==1) ~= 0
                burstIndAttOutTmp(find(diff(burstIndAttOutTmp)==1)+1) = [];
            end
        end
        burstIndAttOutPostCue = [burstIndAttOutPostCue burstIndAttOutTmp];
    end
    burstIndAttOutPostArray = [];
    for n = 1:numel(isiPerTrialAttOutPostArray)
        burstIndAttOutTmp = find(isiPerTrialAttOutPostArray(n).isi<=minTime);
        if length(burstIndAttOutTmp) > 1
            if sum(diff(burstIndAttOutTmp)==1) ~= 0
                burstIndAttOutTmp(find(diff(burstIndAttOutTmp)==1)+1) = [];
            end
        end
        burstIndAttOutPostArray = [burstIndAttOutPostArray burstIndAttOutTmp];
    end
    if ~isempty(burstIndAttOutBaseline) || ~isempty(burstIndAttOutPostCue) || ~isempty(burstIndAttOutPostArray)
        percWDTwoSpikeBurstAttOutBaselinefun = numel(burstIndAttOutBaseline) * 100 / sum(nSpikesAttOutBaseline);
        percWDTwoSpikeBurstAttOutPostCuefun = numel(burstIndAttOutPostCue) * 100 / sum(nSpikesAttOutPostCue);
        percWDTwoSpikeBurstAttOutPostArrayfun = numel(burstIndAttOutPostArray) * 100 / sum(nSpikesAttOutPostArray);
    else
        percWDTwoSpikeBurstAttOutBaselinefun = 0;
        percWDTwoSpikeBurstAttOutPostCuefun = 0;
        percWDTwoSpikeBurstAttOutPostArrayfun = 0;
    end
    isiFoundAttOut = 1;
else
    burstIndAttOutPostCue = [];
    burstIndAttOutPostArray = [];
    burstIndAttOutBaseline = [];
    percWDTwoSpikeBurstAttOutBaselinefun = 0;
    percWDTwoSpikeBurstAttOutPostCuefun = 0;
    percWDTwoSpikeBurstAttOutPostArrayfun = 0;
    isiFoundAttOut = 0;
end

if isiFoundAttIn && isiFoundAttOut
    percWDTwoSpikeBurstBaselinefun = (numel(burstIndAttOutBaseline) + numel(burstIndAttInBaseline)) * 100 / (sum(nSpikesAttOutBaseline) + sum(nSpikesAttInBaseline));
    percWDTwoSpikeBurstPostCuefun = (numel(burstIndAttOutPostCue) + numel(burstIndAttInPostCue)) * 100 / (sum(nSpikesAttOutPostCue) + sum(nSpikesAttInPostCue));
    percWDTwoSpikeBurstPostArrayfun = (numel(burstIndAttOutPostArray) + numel(burstIndAttInPostArray)) * 100 / (sum(nSpikesAttOutPostArray)+ sum(nSpikesAttInPostArray));
else
    percWDTwoSpikeBurstBaselinefun = 0;
    percWDTwoSpikeBurstPostCuefun = 0;
    percWDTwoSpikeBurstPostArrayfun = 0;
end

%% function [firingRatefun, percWDTwoSpikeBurstfun, percWDThreeSpikeBurstfun,...
%     numInterSpikeTimesfun, numBurstLongfun,...
%     numBurstShortfun, frAttInPerTrialfun,frAttOutPerTrialfun,...
%     percWDTwoSpikeBurstAttInfun,percWDTwoSpikeBurstAttOutfun] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
%                 endTimeAttIn, endTimeAttOut, tempResolve, minTime)
% 
%             
% % [firingRate, percWDTwoSpikeBurst, percWDThreeSpikeBurst, numInterSpikeTimes, numBurstLong,...
% %     numBurstShort, frAttInPerTrial,frAttOutPerTrial,...
% %     percWDTwoSpikeBurstAttIn, percWDThreeSpikeBurstAttIn,...
% %     numInterSpikeTimesAttIn, numBurstLongAttIn, numBurstShortAttIn,...
% %     percWDTwoSpikeBurstAttOut, percWDThreeSpikeBurstAttOut,...
% %     numInterSpikeTimesAttOut, numBurstLongAttOut, numBurstShortAttOut] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
% %                 endTimeAttIn, endTimeAttOut)
% %endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;
% 
% % calculate traditional bursting in all spike times
% interSpikeTimes = diff(spikeTimes2use);
% interSpikeTimes = interSpikeTimes*1000; % in ms
% firingRatefun = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));
% 
% burstInd = find(interSpikeTimes<=minTime);
% percWDTwoSpikeBurstfun = numel(burstInd) * 100 / numel(interSpikeTimes);
% 
% if ~isempty(burstInd)
%     if sum(burstInd+1 >= numel(interSpikeTimes)); burstInd = burstInd(1:end-1); end
%     burstIndLong = find(interSpikeTimes(burstInd+1)<=minTime);
%     percWDThreeSpikeBurstfun = numel(burstIndLong) * 100 / numel(interSpikeTimes);
%     numInterSpikeTimesfun = numel(interSpikeTimes);
%     numBurstLongfun = numel(burstIndLong);
%     numBurstShortfun = numel(burstInd);
% else
%     percWDThreeSpikeBurstfun = 0;
%     numInterSpikeTimesfun = numel(interSpikeTimes);
%     numBurstLongfun = 0;
%     numBurstShortfun = 0;
% end
% 
% % repeat for the task times
% data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
% binarySpikeTrain = zeros(1,length(data_pts));
% binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;
% 
% % % binarySpikeTrain = zeros(1,length(data_pts));
% % % binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,4),.0000000582)) = 1;
% % % sum(binarySpikeTrain) 
% % % numel(spikeTimes2use)
% % 
% % if numel(spikeTimes2use)==sum(binarySpikeTrain)
% %     discrep = -1;
% % else
% %     discrep = 1;
% % end
% 
% % calculate trial based spiking
% % 200ms post-cue until TARGET DIM
% data = binarySpikeTrain;
% Fs = 1000;
% NEAttIn=length(cueOnsetTrialAttIn); NEAttOut=length(cueOnsetTrialAttOut);
% nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1; nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
% datatmpAttIn=[];datatmpAttOut=[];
% for n=1:NEAttIn
%     %     nwinl=round(0.200*Fs);
%     nwinr=round(endTimeAttIn(n)*Fs);
%     indx=nEAttIn(n):nwinr-1;
%     if length(indx) >1 && indx(end) < size(data,2)
%         datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
%         isiPerTrialAttIn(n).isi = diff(find(data(indx)));
%         nSpikesAttInPerTrial(n) = sum(data(indx));
%         frAttInPerTrialfun(n) = nSpikesAttInPerTrial(n) / ((indx(end) - indx(1))/Fs);
% %         attIn(n).spikeTimes=find(data(indx));
%         clear indx
%     end
% end
% for n=1:NEAttOut;
%     %     nwinl=round(0.200*Fs);
%     nwinr=round(endTimeAttOut(n)*Fs);
%     indx=nEAttOut(n):nwinr-1;
%     if length(indx) >1 && indx(end) < size(data,2)
%         datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
%         isiPerTrialAttOut(n).isi = diff(find(data(indx)));
%         nSpikesAttOutPerTrial(n) = sum(data(indx));
%         frAttOutPerTrialfun(n) = nSpikesAttOutPerTrial(n) / ((indx(end) - indx(1))/Fs);
% %         attOut(n).spikeTimes=find(data(indx));
%         clear indx
%     end
% end
% 
% % for att in
% % per trial
% burstIndAttInPerTrial = [];
% for n = 1:numel(isiPerTrialAttIn)
%     burstIndAttInTmp = find(isiPerTrialAttIn(n).isi<=minTime);
%     burstIndAttInPerTrial = [burstIndAttInPerTrial burstIndAttInTmp];
% end
% if ~isempty(burstIndAttInPerTrial)
%     percWDTwoSpikeBurstAttInfun = numel(burstIndAttInPerTrial) * 100 / sum(nSpikesAttInPerTrial);
% else
%     percWDTwoSpikeBurstAttInfun = 0;
% end
% % 
% % burstIndAttIn = find(datatmpAttIn<=5);
% % percWDTwoSpikeBurstAttIn = numel(burstIndAttIn) * 100 / sum(nSpikesAttInPerTrial);
% % 
% % if ~isempty(burstIndAttIn)
% %     if sum(burstIndAttIn+1 >= numel(datatmpAttIn)); burstIndAttIn = burstIndAttIn(1:end-1); end
% %     burstIndLongAttIn = find(datatmpAttIn(burstIndAttIn+1)<=5);
% %     percWDThreeSpikeBurstAttIn = numel(burstIndLongAttIn) * 100 / sum(nSpikesAttInPerTrial);
% %     numInterSpikeTimesAttIn = numel(datatmpAttIn);
% %     numBurstLongAttIn = numel(burstIndLongAttIn);
% %     numBurstShortAttIn = numel(burstIndAttIn);
% % else
% %     percWDThreeSpikeBurstAttIn = 0;
% %     numInterSpikeTimesAttIn = numel(datatmpAttIn);
% %     numBurstLongAttIn = 0;
% %     numBurstShortAttIn = 0;
% % end
% 
% % for att out
% burstIndAttOutPerTrial = [];
% for n = 1:numel(isiPerTrialAttOut)
%     burstIndAttOutTmp = find(isiPerTrialAttOut(n).isi<=minTime);
%     burstIndAttOutPerTrial = [burstIndAttOutPerTrial burstIndAttOutTmp];
% end
% if ~isempty(burstIndAttOutPerTrial)
%     percWDTwoSpikeBurstAttOutfun = numel(burstIndAttOutPerTrial) * 100 / sum(nSpikesAttOutPerTrial);
% else
%     percWDTwoSpikeBurstAttOutfun = 0;
% end
% % burstIndAttOut = find(datatmpAttOut<=5);
% % percWDTwoSpikeBurstAttOut = numel(burstIndAttOut) * 100 / sum(nSpikesAttOutPerTrial);
% % 
% % if ~isempty(burstIndAttOut)
% %     if sum(burstIndAttOut+1 >= numel(datatmpAttOut)); burstIndAttOut = burstIndAttOut(1:end-1); end
% %     burstIndLongAttOut = find(datatmpAttOut(burstIndAttOut+1)<=5);
% %     percWDThreeSpikeBurstAttOut = numel(burstIndLongAttOut) * 100 / sum(nSpikesAttOutPerTrial);
% %     numInterSpikeTimesAttOut = numel(datatmpAttOut);
% %     numBurstLongAttOut = numel(burstIndLongAttOut);
% %     numBurstShortAttOut = numel(burstIndAttOut);
% % else
% %     percWDThreeSpikeBurstAttOut = 0;
% %     numInterSpikeTimesAttOut = numel(datatmpAttOut);
% %     numBurstLongAttOut = 0;
% %     numBurstShortAttOut = 0;
% % end