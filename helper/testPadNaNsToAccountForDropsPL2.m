adInfo.Values = [1 2 4 5 6 9 10];
adInfo.FragTs = [1 4 9]'*0.000025;
adInfo.FragCounts = [2 3 2]';
adInfo.ADFreq = 40000;
assert(isequaln(padNaNsToAccountForDropsPL2(adInfo), [1 2 NaN 4 5 6 NaN NaN 9 10]'));

adInfo.Values = [2 3 4 6];
adInfo.FragTs = [2 6]'*0.000025;
adInfo.FragCounts = [3 1]';
adInfo.ADFreq = 40000;
assert(isequaln(padNaNsToAccountForDropsPL2(adInfo), [NaN 2 3 4 NaN 6]'));

adInfo.Values = [-1 1 2 3 4];
adInfo.FragTs = 0;
adInfo.FragCounts = 5;
adInfo.ADFreq = 40000;
assert(isequaln(padNaNsToAccountForDropsPL2(adInfo), [1 2 3 4]'));

adInfo.Values = [-1 1 2 3 4 6];
adInfo.FragTs = [0 6]'*0.000025;
adInfo.FragCounts = [5 1]';
adInfo.ADFreq = 40000;
assert(isequaln(padNaNsToAccountForDropsPL2(adInfo), [1 2 3 4 NaN 6]'));
