function d = createnonemptydatamatpt(sig,times,window)
% Wrapper to chronux's createdatamatpt() that returns datamat with no NaN/Inf in
% it (in case that ever happens...) and some min spikes per window (currently 0)
% The whole d may be empty though!!! (check with isempty())
% 
% Inputs & Output same as for createdatamatpt()

min_spikes_per_window = 0;

datamat = createdatamatpt(sig,times,window);
d = struct([]); c = 1;
for j = 1:numel(datamat)
    if length(datamat(j).times) >= min_spikes_per_window && ...
                ~any(isnan(datamat(j).times)) && ...
                ~any(isinf(datamat(j).times))
        d(c).times = datamat(j).times; c = c + 1;
    elseif any(isnan(datamat(j).times)) || any(isinf(datamat(j).times))
        fprintf('createdatamatpt() is acting funky: %d=%d ',j,datamat(j).times);
    end
end