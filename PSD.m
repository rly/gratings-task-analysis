function [freq_vector,power_norm]=PSD(lfp)
% Calculates LFP power across multicontact probe depth
% Input: LFP in SAMPLES x CHANNELS
% Output: Normalized LFP power across frequencies, contacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs       = 1000; % Hz
chanN    = size(lfp,2);
% loop through channels
for ch = 1:chanN
    clear x n Spec
    x        = lfp(:,ch);
    n        = size(lfp,1); % Number of data points
    % prep for psd
    nfft     = 2048; %512;
    window   = hanning(nfft);
    nwind    = length(window);
    if n < nwind    % zero-pad x if it has length less than the window length
        x(nwind)=0;  n=nwind;
    end
    noverlap = 1;
    k        = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
    index    = 1:nwind;
    % compute PSD
    %x        = lfp(:,ch);
    Spec     = zeros(nfft,1);
    % loop through windows
    for i=1:k
        xw    = window.*(x(index));
        index = index + (nwind - noverlap);
        Xx    = abs(fft(xw,nfft)).^2;
        Spec  = Spec + Xx;
    end
    % Select first half
    if ~any(any(imag(x)~=0))   % check if x is complex
        if rem(nfft,2)    % nfft odd
            select = (1:(nfft+1)/2)';
        else
            select = (1:nfft/2+1)';
        end
        Spec = Spec(select);
    else
        select = (1:nfft)';
    end
    freq_vector = (select - 1)*Fs/nfft;
    if ch == 1
        power = nan(size(Spec,1),chanN);
    end
    power(:,ch) = Spec;
end
% normalize power @ each frequency relative to power across contacts
power_norm = nan(size(power));
for ch = 1:size(power,2)
    for f = 1:size(power,1)
     %   power_norm(f,ch) = (power(f,ch)*100) / mean(power(f,:));
     %power_norm(f,ch) = power(f,ch)/mean(power(f,:)); 
     power_norm(f,ch) = (( power(f,ch) - mean(power(f,:)) )./ mean(power(f,:)) ) *100;
    end
end

