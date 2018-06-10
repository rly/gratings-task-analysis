function [b,Hd] = alphaLowBetaBandpassFilter
%ALPHALOWBETABANDPASSFILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.1 and the Signal Processing Toolbox 7.3.
% Generated on: 30-May-2018 17:49:50

% FIR Window Bandpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fstop1 = 9;                 % First Stopband Frequency
Fpass1 = 10;                % First Passband Frequency
Fpass2 = 20;                % Second Passband Frequency
Fstop2 = 21;                % Second Stopband Frequency
Dstop1 = 0.001;             % First Stopband Attenuation
Dpass  = 0.00057564620966;  % Passband Ripple
Dstop2 = 0.001;             % Second Stopband Attenuation
flag   = 'noscale';         % Sampling Flag

% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                             1 0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dfilt.dffir(b);

% [EOF]