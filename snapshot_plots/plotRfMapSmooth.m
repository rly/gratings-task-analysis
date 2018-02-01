function [h,rfmapSmooth] = plotRfMapSmooth(rfmapSmooth, rfmapCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset)

% rfmapSmooth(1,401) should go to x = -200, y = 200
% transpose because imagesc puts rows into columns and v.v.  
h = imagesc(((1:mapDim) - mapXOffset) / numPixelsPerDegree / mapScale, ...
        ((1:mapDim) - mapYOffset) / numPixelsPerDegree / mapScale, rfmapSmooth');
xlim(([1 mapDim] - mapXOffset) / numPixelsPerDegree / mapScale);
ylim(([1 mapDim] - mapYOffset) / numPixelsPerDegree / mapScale);
xlabel('Degrees horizontal');
ylabel('Degrees vertical');
colorbar;

[rfmapX, rfmapY] = find(rfmapCount > 0);
plot((rfmapX - mapXOffset) / numPixelsPerDegree / mapScale, ...
        (rfmapY - mapYOffset) / numPixelsPerDegree / mapScale, ...
        'o', 'Color', [0.5 0.5 0.5]);
