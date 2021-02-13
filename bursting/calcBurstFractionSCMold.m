function [percentageBurst,percBurstWMAttWhole,percBurstWMAttIn,percBurstWMAttOut] = calcBurstFractionSCMold(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut)

if ~exist('cueOnsetTrialAttIn','var')
    % third parameter does not exist, so default it to something
    param3 = 42;
end

%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;
%endTimeAttIn = targetDimTrialAttIn; endTimeAttOut =targetDimTrialAttOut;
% calculate percentage of spikes >= 200 Hz firing rate
percentageBurst =  sum(diff(spikeTimes2use)<(1/200)) * 100 / size(spikeTimes2use,2);

data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

% calculate trial based spiking
% 200ms post-cue until TARGET DIM
data = binarySpikeTrain;
Fs = 1000;
NEAttIn=length(cueOnsetTrialAttIn); NEAttOut=length(cueOnsetTrialAttOut);
nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1; nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
datatmpAttIn=[];datatmpAttOut=[];
for n=1:NEAttIn
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttIn(n)*Fs);
    indx=nEAttIn(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
        nSpikesAttInPerTrial(n) = sum(data(indx));
        attIn(n).spikeTimes=find(data(indx));
        clear indx
    end
end
for n=1:NEAttOut;
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        attOut(n).spikeTimes=find(data(indx));
        clear indx
    end
end
if exist('attIn','var') && exist('attOut','var')
    NEAttIn = length(attIn); NEAttOut = length(attOut);

    % Perform spike count matching
    clear datatmpAttOut2use datatmpAttIn2use
    nSpikesAttIn = length(datatmpAttIn);
    nSpikesAttOut = length(datatmpAttOut);
    counter = 0;
    if nSpikesAttIn > nSpikesAttOut
        datatmpAttInOld = datatmpAttIn;
        trialSel = randperm(NEAttIn);
        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
            counter = counter + 1;
            NEAttInNw = NEAttIn - counter;
            nEAttInNw = nEAttIn(trialSel(1:NEAttInNw)); endTimeAttInNw = endTimeAttIn(trialSel(1:NEAttInNw));
            for n=1:NEAttInNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(endTimeAttInNw(n)*Fs);
                indx=nEAttInNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                    attInNw(n).spikeTimes=find(data(indx));
                end
                clear indx
            end
            datatmpAttIn2use = datatmpAttIn;
            nSpikesAttIn = length(datatmpAttIn2use);
            datatmpAttIn = [];
        end
    elseif nSpikesAttOut > nSpikesAttIn
        datatmpAttOutOld = datatmpAttOut; datatmpAttOut = [];
        trialSel = randperm(NEAttOut);
        while nSpikesAttOut + mean(nSpikesAttOutPerTrial) / 2 > nSpikesAttIn
            counter = counter + 1;
            NEAttOutNw = NEAttOut - counter;
            nEAttOutNw = nEAttOut(trialSel(1:NEAttOutNw)); endTimeAttOutNw = endTimeAttOut(trialSel(1:NEAttOutNw));
            for n=1:NEAttOutNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(endTimeAttOutNw(n)*Fs);
                indx=nEAttOutNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                    attOutNw(n).spikeTimes=find(data(indx));
                end
                clear indx
            end
            datatmpAttOut2use = datatmpAttOut;
            nSpikesAttOut = length(datatmpAttOut2use);
            datatmpAttOut = [];
        end
    end

    if exist('datatmpAttIn2use','var') == 0 && exist('datatmpAttOut2use','var') == 1
        datatmpAttIn2use = datatmpAttIn;
        endTimeAttInNw = endTimeAttIn;
        nEAttInNw = nEAttIn; NEAttInNw = NEAttIn;
        if exist('trialSel','var') == 0
            cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(1:NEAttOutNw);
        else
            cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(trialSel(1:NEAttOutNw));
        end
        cueOnsetTrialAttInNw = cueOnsetTrialAttIn;
        attInNw = attIn;
    elseif exist('datatmpAttIn2use','var') == 1 && exist('datatmpAttOut2use','var') == 0
        datatmpAttOut2use = datatmpAttOut;
        endTimeAttOutNw = endTimeAttOut;
        nEAttOutNw = nEAttOut; NEAttOutNw = NEAttOut;
        if exist('trialSel','var') == 0
            cueOnsetTrialAttInNw = cueOnsetTrialAttIn(1:NEAttInNw);
        else
            cueOnsetTrialAttInNw = cueOnsetTrialAttIn(trialSel(1:NEAttInNw));
        end
        cueOnsetTrialAttOutNw = cueOnsetTrialAttOut;
        attOutNw = attOut;
    else
        datatmpAttIn2use = datatmpAttIn;
        endTimeAttInNw = endTimeAttIn;
        nEAttInNw = nEAttIn; NEAttInNw = NEAttIn;
        datatmpAttOut2use = datatmpAttOut;
        endTimeAttOutNw = endTimeAttOut;
        nEAttOutNw = nEAttOut; NEAttOutNw = NEAttOut;
        if exist('trialSel','var') == 0
            cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(1:NEAttOutNw);
        else
            cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(trialSel(1:NEAttOutNw));
        end
        cueOnsetTrialAttInNw = cueOnsetTrialAttIn;
        attInNw = attIn;
        if exist('trialSel','var') == 0
            cueOnsetTrialAttInNw = cueOnsetTrialAttIn(1:NEAttInNw);
        else
            cueOnsetTrialAttInNw = cueOnsetTrialAttIn(trialSel(1:NEAttInNw));
        end
        cueOnsetTrialAttOutNw = cueOnsetTrialAttOut;
        attOutNw = attOut;
    end
    totalAttConds = [datatmpAttIn2use datatmpAttOut2use];
    percBurstWMAttWhole.four = sum(totalAttConds <= 4)  / size(totalAttConds,2);
    percBurstWMAttWhole.five = sum(totalAttConds <= 5)  / size(totalAttConds,2);
    percBurstWMAttWhole.ten = sum(totalAttConds <= 10)  / size(totalAttConds,2);
    percBurstWMAttIn.four = sum(datatmpAttIn2use <= 4)  / size(datatmpAttIn2use,2);
    percBurstWMAttIn.five = sum(datatmpAttIn2use <= 5)  / size(datatmpAttIn2use,2);
    percBurstWMAttIn.ten = sum(datatmpAttIn2use <= 10)  / size(datatmpAttIn2use,2);
    percBurstWMAttOut.four = sum(datatmpAttOut2use <= 4)  / size(datatmpAttOut2use,2);
    percBurstWMAttOut.five = sum(datatmpAttOut2use <= 5)  / size(datatmpAttOut2use,2);
    percBurstWMAttOut.ten = sum(datatmpAttOut2use <= 10)  / size(datatmpAttOut2use,2);
else
    percBurstWMAttWhole.four = 0;
    percBurstWMAttWhole.five = 0;
    percBurstWMAttWhole.ten = 0;
    percBurstWMAttIn.four = 0;
    percBurstWMAttIn.five = 0;
    percBurstWMAttIn.ten = 0;
    percBurstWMAttOut.four = 0;
    percBurstWMAttOut.five = 0;
    percBurstWMAttOut.ten = 0;
end