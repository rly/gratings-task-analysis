function [f1,ax1] = grcPlotNoiseCorrDiff(corrAllUnits, goodUnitsDPul, goodUnitsVPul)
corrP3DPul = [];
corrP1DPul = [];
corrP3VPul = [];
corrP1VPul = [];

if any(goodUnitsDPul)
corrP3DPul = squeeze(corrAllUnits(goodUnitsDPul,goodUnitsDPul,3));
corrP1DPul = squeeze(corrAllUnits(goodUnitsDPul,goodUnitsDPul,1));
corrP3DPul = corrP3DPul(:);
corrP1DPul = corrP1DPul(:);
corrP3DPul(isnan(corrP3DPul)) = [];
corrP1DPul(isnan(corrP1DPul)) = [];
end
if any(goodUnitsVPul)
corrP3VPul = squeeze(corrAllUnits(goodUnitsVPul,goodUnitsVPul,3));
corrP1VPul = squeeze(corrAllUnits(goodUnitsVPul,goodUnitsVPul,1));
corrP3VPul = corrP3VPul(:);
corrP1VPul = corrP1VPul(:);
corrP3VPul(isnan(corrP3VPul)) = [];
corrP1VPul(isnan(corrP1VPul)) = [];
end
% D always comes before V per session so don't need the reserve (which should be all blank)
corrP3Cross = [];%squeeze(arrayResponseHoldLateNoiseCorr(goodUnitsDPul,goodUnitsVPul,3));
corrP1Cross = [];%squeeze(arrayResponseHoldLateNoiseCorr(goodUnitsDPul,goodUnitsVPul,1));
corrP3Cross = corrP3Cross(:);
corrP1Cross = corrP1Cross(:);
corrP3Cross(isnan(corrP3Cross)) = [];
corrP1Cross(isnan(corrP1Cross)) = [];

corrP3 = [corrP3DPul; corrP3VPul; corrP3Cross];
corrP1 = [corrP1DPul; corrP1VPul; corrP1Cross];
isDPul = [true(size(corrP3DPul)); false(size(corrP3VPul)); false(size(corrP3Cross))];
isVPul = [false(size(corrP3DPul)); true(size(corrP3VPul)); false(size(corrP3Cross))];

[f1,ax1] = grcPlotMetricDiff(corrP3, corrP1, isDPul, isVPul, false(size(isDPul)) );
