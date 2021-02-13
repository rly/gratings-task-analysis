function [s] = makeSDF(spikes,sigma)

mirror = fliplr(spikes);
padded = [mirror spikes mirror];
%%
% convolve psth with kernel
kwidth = -sigma*3:3*sigma;
kernel = normpdf(kwidth,0,sigma);
s = conv(padded,kernel);
s(1:length(mirror)) = [];
s(length(s)-length(mirror)+1:length(s)) = [];
s(1:floor(length(kernel)/2)) = [];
s(length(spikes)+1:length(s)) = [];