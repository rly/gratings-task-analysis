function dataOut=createEventLockedGAKS(data,E,Fs,win)
% adjusted createdatamatc from chronux here to make the formula spit out
% event locked data in trial format

% data = allSpikesGaussSep;
% E = UE.cueOnset-spikeTimes(1);
% Fs = D.directFs;
% win = allSpikesGaussSepCueLocked.window;
NE=length(E);
nwinl=round(win(1)*Fs);
nwinr=round(win(2)*Fs);
nE=floor(E*Fs)+1;
datatmp=[];
for n=1:NE;
    indx=nE(n)-nwinl:nE(n)+nwinr-1;
    datatmp=[datatmp; data(indx)];
end
dataOut.eventLockedGAKS = datatmp;
dataOut.window = win;