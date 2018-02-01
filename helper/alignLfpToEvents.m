function [aligned,t] = alignLfpToEvents(data, events, Fs, windowOffset)
% data: nChannel x nSample

nwinl = round(windowOffset(1) * Fs);
nwinr = round(windowOffset(2) * Fs) - 1;
eventIndex = round(events * Fs);
relIndex = nwinl:nwinr;
t = relIndex / Fs;

% nChannel x nEvent x nSample
aligned = nan(size(data, 1), numel(events), numel(t));
for n = 1:numel(events)
    indx = eventIndex(n) + relIndex;
    aligned(:,n,:) = data(:,indx);
end
