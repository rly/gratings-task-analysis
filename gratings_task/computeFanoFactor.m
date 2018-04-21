function ff = computeFanoFactor(countStruct, loc)
% compute fano factor given a averageFiringRatesByCount struct
% that struct is already in units of Hz. convert back to spikes by
% multiplying by the time window
ff = countStruct.byLocSD(loc)^2 / countStruct.byLoc(loc) * diff(countStruct.windowOffset);