function t = computeTForLfp(analysisWindowOffset, Fs)

t = round(analysisWindowOffset(1) * Fs):round(analysisWindowOffset(2) * Fs);