function [firingRatefun, percClassicTwoSpikeBurstfun, percClassicThreeSpikeBurstfun, numInterSpikeTimesfun, numBurstLongfun,...
    numBurstShortfun, frAttInPerTrialfun,frAttOutPerTrialfun,...
    percClassicTwoSpikeBurstAttInfun, percClassicTwoSpikeBurstAttOutfun,...
    percClassicTwoSpikeBurstBaselinefun, percClassicTwoSpikeBurstAttInBaselinefun, percClassicTwoSpikeBurstAttOutBaselinefun,...
    percClassicTwoSpikeBurstPostCuefun, percClassicTwoSpikeBurstAttInPostCuefun, percClassicTwoSpikeBurstAttOutPostCuefun,...
    percClassicTwoSpikeBurstPostArrayfun, percClassicTwoSpikeBurstAttInPostArrayfun, percClassicTwoSpikeBurstAttOutPostArrayfun,...
    cueOnsetfun,arrayOnset,targetDimOnset] = histClassicAttCondi_homepc(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut, targetDimTrialAttIn,targetDimTrialAttOut, tempResolve, minTime, maxTime, UE, startTime, endTime)

      
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

% indices of 100ms quiet periods
quietInd = find(interSpikeTimes>maxTime);
% find times after 100ms quiet that are <= to 4ms
if quietInd(end) == numel(interSpikeTimes); quietInd = quietInd(1:end-1); end
burstInd = find(interSpikeTimes(quietInd+1)<=minTime);
percClassicTwoSpikeBurstfun = numel(interSpikeTimes(quietInd(burstInd))) * 100 / numel(spikeTimes2use);

if ~isempty(burstInd)
    if quietInd(burstInd(end))+2 >= numel(interSpikeTimes); burstInd = burstInd(1:end-1); end
    burstIndLong = find(interSpikeTimes(quietInd(burstInd)+2)<=minTime);
    percClassicThreeSpikeBurstfun = numel(burstIndLong) * 100 / numel(spikeTimes2use);
    numInterSpikeTimesfun = numel(interSpikeTimes);
    numBurstLongfun = numel(burstIndLong);
    numBurstShortfun = numel(burstInd);
else
    percClassicThreeSpikeBurstfun = 0;
    numInterSpikeTimesfun = numel(interSpikeTimes);
    numBurstLongfun = 0;
    numBurstShortfun = 0;
end

% repeat for the task times
data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

kernelSigma = .001;
if tempResolve
    cueOnsetfun.window = [0.8 0.8]; % seconds before, after
    cueOnsetfun.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
    cueOnsetfun = createTimeLockedSpdf(spikeTimes2use(quietInd(burstInd)+1), UE.cueOnset, UE.cueOnsetByLoc, cueOnsetfun, kernelSigma, startTime, endTime);
    arrayOnset.window = [0.8 0.8]; % seconds before, after
    arrayOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
    arrayOnset = createTimeLockedSpdf(spikeTimes2use(quietInd(burstInd)+1), UE.arrayOnset, UE.arrayOnsetByLoc, arrayOnset, kernelSigma, startTime, endTime);
    targetDimOnset.window = [0.8 0.8]; % seconds before, after
    targetDimOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
    targetDimOnset = createTimeLockedSpdf(spikeTimes2use(quietInd(burstInd)+1), UE.targetDim, UE.targetDimBalByLoc, targetDimOnset, kernelSigma, startTime, endTime);    
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

isiPerTrialAttInfun = struct(); isiPerTrialAttInBaseline = struct(); isiPerTrialAttInPostCue = struct(); isiPerTrialAttInPostArray = struct();
nSpikesAttInPerTrial = nan(NEAttIn,1); nSpikesAttInBaseline = nan(NEAttIn,1); nSpikesAttInPostCue = nan(NEAttIn,1); nSpikesAttInPostArray = nan(NEAttIn,1);
frAttInPerTrialfun = nan(NEAttIn,1); frAttInPerTrialBaselinefun = nan(NEAttIn,1);frAttInPerTrialPostCuefun = nan(NEAttIn,1);frAttInPerTrialPostArrayfun = nan(NEAttIn,1);
isiPerTrialAttOutfun = struct(); isiPerTrialAttOutBaseline = struct(); isiPerTrialAttOutPostCue = struct(); isiPerTrialAttOutPostArray = struct();
nSpikesAttOutPerTrial = nan(NEAttOut,1); nSpikesAttOutBaseline = nan(NEAttOut,1); nSpikesAttOutPostCue = nan(NEAttOut,1); nSpikesAttOutPostArray = nan(NEAttOut,1);
frAttOutPerTrialfun = nan(NEAttOut,1); frAttOutPerTrialBaselinefun = nan(NEAttOut,1);frAttOutPerTrialPostCuefun = nan(NEAttOut,1);frAttOutPerTrialPostArrayfun = nan(NEAttOut,1);

for n=1:NEAttIn
    nwinr=round(endTimeAttIn(n)*Fs); % arrayOnset
    indx=nEAttIn(n):nwinr-1; % cueOnset:arrayOnset
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
        isiPerTrialAttInfun(n).isi = diff(find(data(indx)));
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
        isiPerTrialAttOutfun(n).isi = diff(find(data(indx)));
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
for n = 1:numel(isiPerTrialAttInfun)
    qPeriod = isiPerTrialAttInfun(n).isi>maxTime;
    if sum(qPeriod)>=1
        if sum(find(qPeriod)+1<=numel(isiPerTrialAttInfun(n).isi)) < sum(qPeriod)
            qPeriod = qPeriod(1:end-1);
        end
        if isiPerTrialAttInfun(n).isi(find(qPeriod)+1) <= minTime
            burstIndAttInTmp = [find(qPeriod) find(qPeriod)+1];
            burstIndAttInPerTrial = [burstIndAttInPerTrial burstIndAttInTmp];
        end
    end
end
if ~isempty(burstIndAttInPerTrial)
    percClassicTwoSpikeBurstAttInfun = numel(burstIndAttInPerTrial) * 100 / sum(nSpikesAttInPerTrial);
else
    percClassicTwoSpikeBurstAttInfun = 0;
end
    

if ~isempty(isiPerTrialAttInBaseline) && ~isempty(isiPerTrialAttInPostArray) && ~isempty(isiPerTrialAttInPostCue)
    % att in baseline
    burstIndAttInBaseline = [];
    for n = 1:numel(isiPerTrialAttInBaseline)
        qPeriod = isiPerTrialAttInBaseline(n).isi>maxTime;
        if sum(qPeriod)>=1
            if sum(find(qPeriod)+1<=numel(isiPerTrialAttInBaseline(n).isi)) < sum(qPeriod)
                qPeriod = qPeriod(1:end-1);
            end
            if isiPerTrialAttInBaseline(n).isi(find(qPeriod)+1) <= minTime
                burstIndAttInTmp = [find(qPeriod) find(qPeriod)+1];
                burstIndAttInBaseline = [burstIndAttInBaseline burstIndAttInTmp];
            end
        end
    end
    burstIndAttInPostCue = [];
    for n = 1:numel(isiPerTrialAttInPostCue)
        qPeriod = isiPerTrialAttInPostCue(n).isi>maxTime;
        if sum(qPeriod)>=1
            if sum(find(qPeriod)+1<=numel(isiPerTrialAttInPostCue(n).isi)) < sum(qPeriod)
                qPeriod = qPeriod(1:end-1);
            end
            if isiPerTrialAttInPostCue(n).isi(find(qPeriod)+1) <= minTime
                burstIndAttInTmp = [find(qPeriod) find(qPeriod)+1];
                burstIndAttInPostCue = [burstIndAttInPostCue burstIndAttInTmp];
            end
        end
    end
    burstIndAttInPostArray = [];
    for n = 1:numel(isiPerTrialAttInPostArray)
        qPeriod = isiPerTrialAttInPostArray(n).isi>maxTime;
        if sum(qPeriod)>=1
            if sum(find(qPeriod)+1<=numel(isiPerTrialAttInPostArray(n).isi)) < sum(qPeriod)
                qPeriod = qPeriod(1:end-1);
            end
            if isiPerTrialAttInPostArray(n).isi(find(qPeriod)+1) <= minTime
                burstIndAttInTmp = [find(qPeriod) find(qPeriod)+1];
                burstIndAttInPostArray = [burstIndAttInPostArray burstIndAttInTmp];
            end
        end
    end
    if ~isempty(burstIndAttInBaseline) || ~isempty(burstIndAttInPostCue) || ~isempty(burstIndAttInPostArray)
        percClassicTwoSpikeBurstAttInBaselinefun = numel(burstIndAttInBaseline) * 100 / sum(nSpikesAttInBaseline);
        percClassicTwoSpikeBurstAttInPostCuefun = numel(burstIndAttInPostCue) * 100 / sum(nSpikesAttInPostCue);
        percClassicTwoSpikeBurstAttInPostArrayfun = numel(burstIndAttInPostArray) * 100 / sum(nSpikesAttInPostArray);
    else
        percClassicTwoSpikeBurstAttInBaselinefun = 0;
        percClassicTwoSpikeBurstAttInPostCuefun = 0;
        percClassicTwoSpikeBurstAttInPostArrayfun = 0;
    end
    isiFoundAttIn = 1;
else
    burstIndAttInBaseline = [];
    burstIndAttInPostCue = [];
    burstIndAttInPostArray = [];
    percClassicTwoSpikeBurstAttInBaselinefun = 0;
    percClassicTwoSpikeBurstAttInPostCuefun = 0;
    percClassicTwoSpikeBurstAttInPostArrayfun = 0;
    isiFoundAttIn = 0;
end

% for att out
burstIndAttOutPerTrial = [];
for n = 1:numel(isiPerTrialAttOutfun)
    qPeriod = isiPerTrialAttOutfun(n).isi>maxTime;
    if sum(qPeriod)>=1
        if sum(find(qPeriod)+1<=numel(isiPerTrialAttOutfun(n).isi)) < sum(qPeriod)
            qPeriod = qPeriod(1:end-1);
        end
        if isiPerTrialAttOutfun(n).isi(find(qPeriod)+1) <= minTime
            burstIndAttOutTmp = [find(qPeriod) find(qPeriod)+1];
            burstIndAttOutPerTrial = [burstIndAttOutPerTrial burstIndAttOutTmp];
        end
    end
end
if ~isempty(burstIndAttOutPerTrial)
    percClassicTwoSpikeBurstAttOutfun = numel(burstIndAttOutPerTrial) * 100 / sum(nSpikesAttOutPerTrial);
else
    percClassicTwoSpikeBurstAttOutfun = 0;
end

if ~isempty(isiPerTrialAttOutBaseline) && ~isempty(isiPerTrialAttOutPostArray) && ~isempty(isiPerTrialAttOutPostCue)
    % att out baseline
    burstIndAttOutBaseline = [];
    for n = 1:numel(isiPerTrialAttOutBaseline)
        qPeriod = isiPerTrialAttOutBaseline(n).isi>maxTime;
        if sum(qPeriod)>=1
            if sum(find(qPeriod)+1<=numel(isiPerTrialAttOutBaseline(n).isi)) < sum(qPeriod)
                qPeriod = qPeriod(1:end-1);
            end
            if isiPerTrialAttOutBaseline(n).isi(find(qPeriod)+1) <= minTime
                burstIndAttOutTmp = [find(qPeriod) find(qPeriod)+1];
                burstIndAttOutBaseline = [burstIndAttOutBaseline burstIndAttOutTmp];
            end
        end
    end
    burstIndAttOutPostCue = [];
    for n = 1:numel(isiPerTrialAttOutPostCue)
        qPeriod = isiPerTrialAttOutPostCue(n).isi>maxTime;
        if sum(qPeriod)>=1
            if sum(find(qPeriod)+1<=numel(isiPerTrialAttOutPostCue(n).isi)) < sum(qPeriod)
                qPeriod = qPeriod(1:end-1);
            end
            if isiPerTrialAttOutPostCue(n).isi(find(qPeriod)+1) <= minTime
                burstIndAttOutTmp = [find(qPeriod) find(qPeriod)+1];
                burstIndAttOutPostCue = [burstIndAttOutPostCue burstIndAttOutTmp];
            end
        end
    end
    burstIndAttOutPostArray = [];
    for n = 1:numel(isiPerTrialAttOutPostArray)
        qPeriod = isiPerTrialAttOutPostArray(n).isi>maxTime;
        if sum(qPeriod)>=1
            if sum(find(qPeriod)+1<=numel(isiPerTrialAttOutPostArray(n).isi)) < sum(qPeriod)
                qPeriod = qPeriod(1:end-1);
            end
            if isiPerTrialAttOutPostArray(n).isi(find(qPeriod)+1) <= minTime
                burstIndAttOutTmp = [find(qPeriod) find(qPeriod)+1];
                burstIndAttOutPostArray = [burstIndAttOutPostArray burstIndAttOutTmp];
            end
        end
    end
    if ~isempty(burstIndAttOutBaseline) || ~isempty(burstIndAttOutPostCue) || ~isempty(burstIndAttOutPostArray)
        percClassicTwoSpikeBurstAttOutBaselinefun = numel(burstIndAttOutBaseline) * 100 / sum(nSpikesAttOutBaseline);
        percClassicTwoSpikeBurstAttOutPostCuefun = numel(burstIndAttOutPostCue) * 100 / sum(nSpikesAttOutPostCue);
        percClassicTwoSpikeBurstAttOutPostArrayfun = numel(burstIndAttOutPostArray) * 100 / sum(nSpikesAttOutPostArray);
    else
        percClassicTwoSpikeBurstAttOutBaselinefun = 0;
        percClassicTwoSpikeBurstAttOutPostCuefun = 0;
        percClassicTwoSpikeBurstAttOutPostArrayfun = 0;
    end
    isiFoundAttOut = 1;
else
    burstIndAttOutPostCue = [];
    burstIndAttOutPostArray = [];
    burstIndAttOutBaseline = [];
    percClassicTwoSpikeBurstAttOutBaselinefun = 0;
    percClassicTwoSpikeBurstAttOutPostCuefun = 0;
    percClassicTwoSpikeBurstAttOutPostArrayfun = 0;
    isiFoundAttOut = 0;
end

if isiFoundAttIn && isiFoundAttOut
    percClassicTwoSpikeBurstBaselinefun = (numel(burstIndAttOutBaseline) + numel(burstIndAttInBaseline)) * 100 / (sum(nSpikesAttOutBaseline) + sum(nSpikesAttInBaseline));
    percClassicTwoSpikeBurstPostCuefun = (numel(burstIndAttOutPostCue) + numel(burstIndAttInPostCue)) * 100 / (sum(nSpikesAttOutPostCue) + sum(nSpikesAttInPostCue));
    percClassicTwoSpikeBurstPostArrayfun = (numel(burstIndAttOutPostArray) + numel(burstIndAttInPostArray)) * 100 / (sum(nSpikesAttOutPostArray) + sum(nSpikesAttInPostArray));
else
    percClassicTwoSpikeBurstBaselinefun = 0;
    percClassicTwoSpikeBurstPostCuefun = 0;
    percClassicTwoSpikeBurstPostArrayfun = 0;
end

% % indices of 100ms quiet periods
% quietIndAttIn = find(datatmpAttIn>100);
% % find times after 100ms quiet that are <= to 4ms
% if ~isempty(quietIndAttIn)
%     if quietIndAttIn(end) == numel(datatmpAttIn); quietIndAttIn = quietIndAttIn(1:end-1); end
%     burstIndAttIn = find(datatmpAttIn(quietIndAttIn+1)<=4);
%     percClassicTwoSpikeBurstAttInfun = numel(datatmpAttIn(quietIndAttIn(burstIndAttIn))) * 100 / sum(nSpikesAttInPerTrial);
% else
%     burstIndAttIn = [];
%     percClassicTwoSpikeBurstAttInfun = NaN;
% end

% if ~isempty(burstIndAttIn)
%     if quietIndAttIn(burstIndAttIn(end))+2 >= numel(datatmpAttIn); burstIndAttIn = burstIndAttIn(1:end-1); end
%     burstIndLongAttIn = find(datatmpAttIn(quietIndAttIn(burstIndAttIn)+2)<=4);
%     percClassicThreeSpikeBurstAttInfun = numel(burstIndLongAttIn) * 100 / sum(nSpikesAttInPerTrial);
%     numInterSpikeTimesAttInfun = numel(datatmpAttIn);
%     numBurstLongAttInfun = numel(burstIndLongAttIn);
%     numBurstShortAttInfun = numel(burstIndAttIn);
% else
%     percClassicThreeSpikeBurstAttInfun = 0;
%     numInterSpikeTimesAttInfun = numel(datatmpAttIn);
%     numBurstLongAttInfun = 0;
%     numBurstShortAttInfun = 0;
% end




% % indices of 100ms quiet periods
% quietIndAttOut = find(datatmpAttOut>100);
% % find times after 100ms quiet that are <= to 4ms
% if ~isempty(quietIndAttOut)
%     if quietIndAttOut(end) == numel(datatmpAttOut); quietIndAttOut = quietIndAttOut(1:end-1); end
%     burstIndAttOut = find(datatmpAttOut(quietIndAttOut+1)<=4);
%     percClassicTwoSpikeBurstAttOutfun = numel(datatmpAttOut(quietIndAttOut(burstIndAttOut))) * 100 / sum(nSpikesAttOutPerTrial);
% else
%     burstIndAttOut = [];
%     percClassicTwoSpikeBurstAttOutfun = NaN;
% end
% 
% if ~isempty(burstIndAttOut)
%     if quietIndAttOut(burstIndAttOut(end))+2 >= numel(datatmpAttOut); burstIndAttOut = burstIndAttOut(1:end-1); end
%     burstIndLongAttOut = find(datatmpAttOut(quietIndAttOut(burstIndAttOut)+2)<=4);
%     percClassicThreeSpikeBurstAttOutfun = numel(burstIndLongAttOut) * 100 / sum(nSpikesAttOutPerTrial);
%     numInterSpikeTimesAttOutfun = numel(datatmpAttOut);
%     numBurstLongAttOutfun = numel(burstIndLongAttOut);
%     numBurstShortAttOutfun = numel(burstIndAttOut);
% else
%     percClassicThreeSpikeBurstAttOutfun = 0;
%     numInterSpikeTimesAttOutfun = numel(datatmpAttOut);
%     numBurstLongAttOutfun = 0;
%     numBurstShortAttOutfun = 0;
% end
