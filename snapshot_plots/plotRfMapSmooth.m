function [h,rfmapSmooth] = plotRfMapSmooth(rfmapSmooth, rfmapCount, ...
        numPixelsPerDegree, mapScale)
% plot smoothed heatmap represented with circles at sampled locations, in
% units of degrees of visual angle
%
% Inputs: 
% - rfmapSmooth - N x N matrix representing smoothed responses over space,
%                 centered at (0,0)
% - rfmapCount - N x N matrix with number of samples at each location
% - numPixelsPerDegree - number of screen pixels per degree visual angle
% - mapScale - scaling factor of number of points, or units, per pixel
    
% rfmapSmooth(1,401) should go to x = -200, y = 200
% transpose because imagesc puts rows into columns and v.v.  
mapDim = size(rfmapSmooth, 1);
mapDimY = size(rfmapSmooth, 2);
assert(mapDim == mapDimY);
mapXOffset = (mapDim + 1)/2;
mapYOffset = (mapDim + 1)/2;

% plot heatmap, centered at (0,0), in degrees of visual angle
h = imagesc(((1:mapDim) - mapXOffset) / numPixelsPerDegree / mapScale, ...
        ((1:mapDim) - mapYOffset) / numPixelsPerDegree / mapScale, rfmapSmooth');

% formatting
xlim(([1 mapDim] - mapXOffset) / numPixelsPerDegree / mapScale);
ylim(([1 mapDim] - mapYOffset) / numPixelsPerDegree / mapScale);
xlabel('Degrees horizontal');
ylabel('Degrees vertical');
colorbar;

% plot open circles at sampled locations
[rfmapX, rfmapY] = find(rfmapCount > 0);
plot((rfmapX - mapXOffset) / numPixelsPerDegree / mapScale, ...
        (rfmapY - mapYOffset) / numPixelsPerDegree / mapScale, ...
        'o', 'Color', [0.5 0.5 0.5]);
