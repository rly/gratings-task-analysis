function t = computeTForSpdf(eventTime, analysisWindowOffset, kernelSigma)

analysisWindow = eventTime + analysisWindowOffset;
totalWindowTime = analysisWindow(2) - analysisWindow(1);

% create evenly distributed time vector with spacing kernelSigma / 10.
% this way, when each spike is replaced by a gaussian at the sampled time
% points, a spike that occurs not at a sampled time will generate a peak
% value that is 98% of the peak value generated from a spike that occurs
% exactly at a sampled time. replace 10 with a higher value to reduce this
% sampling frequency effect, but at the cost of increasing computation time
% 
% in the original chronux psth(), if kernelSigma is 0.01, then time points
% are separated by about 2ms. here time points are separated by 1ms.
nTime = fix(10 * totalWindowTime / kernelSigma) + 1; % number of time steps
t = linspace(analysisWindow(1), analysisWindow(2), nTime);