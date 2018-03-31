function timeLockStruct = alignLfpToEvents(data, events, Fs, timeLockStruct)
% data: nChannel x nSample

nwinl = round(timeLockStruct.windowOffset(1) * Fs);
nwinr = round(timeLockStruct.windowOffset(2) * Fs) - 1;
eventIndex = round(events * Fs);
relIndex = nwinl:nwinr;
timeLockStruct.t = relIndex / Fs;
timeLockStruct.Fs = Fs;

% nChannel x nEvent x nSample
timeLockStruct.aligned = nan(size(data, 1), numel(events), numel(timeLockStruct.t));
for n = 1:numel(events)
    indx = eventIndex(n) + relIndex;
    timeLockStruct.lfp(:,n,:) = data(:,indx);
end
