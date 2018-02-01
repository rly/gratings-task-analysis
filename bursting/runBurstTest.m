function runBurstTest(spikesInTrials, N, maxLag)

disp('-----------------------------');
disp('spikesInTrials:');
disp(spikesInTrials);
disp('SPNA:');
SPNA = computeSPNA(spikesInTrials, N, maxLag);
disp([0:maxLag; SPNA(maxLag+1:end)]')

end