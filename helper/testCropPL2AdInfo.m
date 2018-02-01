tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 5.1/40000, 9.1/40000);
assert(isequal(newAdInfo.Values, [6 9]));
assert(all(abs(newAdInfo.FragTs - [6 9]'/40000) < tol));
assert(isequal(newAdInfo.FragCounts, [1 1]'));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 5.1/40000, 9/40000);
assert(isequal(newAdInfo.Values, [6 9]));
assert(all(abs(newAdInfo.FragTs - [6 9]'/40000) < tol));
assert(isequal(newAdInfo.FragCounts, [1 1]'));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 14/40000);
assert(isequal(newAdInfo.Values, adInfo.Values));
assert(all(abs(newAdInfo.FragTs - adInfo.FragTs) < tol));
assert(isequal(newAdInfo.FragCounts, adInfo.FragCounts));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 15/40000);
assert(isequal(newAdInfo.Values, adInfo.Values));
assert(all(abs(newAdInfo.FragTs - adInfo.FragTs) < tol));
assert(isequal(newAdInfo.FragCounts, adInfo.FragCounts));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 13.9/40000);
assert(isequal(newAdInfo.Values, adInfo.Values(1:end-1)));
assert(all(abs(newAdInfo.FragTs - adInfo.FragTs) < tol));
assert(isequal(newAdInfo.FragCounts, [2 3 2 1]'));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 0.9/40000, 14/40000);
assert(isequal(newAdInfo.Values, adInfo.Values));
assert(all(abs(newAdInfo.FragTs - adInfo.FragTs) < tol));
assert(isequal(newAdInfo.FragCounts, adInfo.FragCounts));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 1/40000);
assert(isequal(newAdInfo.Values, 1));
assert(all(abs(newAdInfo.FragTs - 1/40000) < tol));
assert(isequal(newAdInfo.FragCounts, 1));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 1.1/40000);
assert(isequal(newAdInfo.Values, 1));
assert(all(abs(newAdInfo.FragTs - 1/40000) < tol));
assert(isequal(newAdInfo.FragCounts, 1));
assert(newAdInfo.ADFreq == 40000);

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
try
    newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 0.9/40000);
catch ME
    assert(strcmp(ME.identifier,'MATLAB:assertion:failed'));
end

tol = 1e8;
adInfo.Values = [1 2 4 5 6 9 10 NaN 13 14];
adInfo.FragTs = [1 4 9 13]'/40000;
adInfo.FragCounts = [2 3 2 2]';
adInfo.ADFreq = 40000;
try
    newAdInfo = cropPL2AdInfo(adInfo, 1/40000, 0.9/40000);
catch ME
    assert(strcmp(ME.identifier,'MATLAB:assertion:failed'));
end