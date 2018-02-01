function adjLfp = padNaNsToAdjustLfpOffset(lfp, lfp_ts, lfp_ind, Fs)

% add NaNs at start of AD01 and within AD01 so that
% AD01(round(timestamp*Fs)) will yield the actual LFP signal at that time
% also overwrite data points if the next block begins before all of the
% data in the current block has been tagged with a timestamp.

adjLfp = nan(size(lfp));
% add dummy value as if there were another block at the end
lfp_ind(end+1) = numel(lfp) + 1;




for i = 1:numel(lfp_ts)
    startBlockIndex = round(lfp_ts(i)*Fs);
    % special case: lfp_ts(1) = 0 -- just remove the first data point t=0
    if startBlockIndex == 0 && lfp_ind(i) == 1
        startBlockIndex = 1;
        lfp_ind(i) = 2;
    end
    numValuesInBlock = lfp_ind(i+1) - lfp_ind(i) - 1;
    adjLfp(startBlockIndex:startBlockIndex+numValuesInBlock) = lfp(lfp_ind(i):(lfp_ind(i+1)-1));
end