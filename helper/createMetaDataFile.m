function MD = createMetaDataFile(D, metaDataFilePath)
% save just a few vars for more efficient access later

MD.blockStartTimes = D.blockStartTimes;
MD.blockStopTimes = D.blockStopTimes;
if isfield(D, 'allUnitStructs')
    MD.allUnitStructs = D.allUnitStructs;
end
save(metaDataFilePath, 'MD');