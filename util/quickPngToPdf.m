function pdfFilePath = quickPngToPdf(pngFilePath, pdfFilePath)
% converts a png file into a pdf file (this loses vector graphics)

assert(strcmp(pngFilePath(end-3:end), '.png'));
if nargin < 2 || isempty(pdfFilePath)
    pdfFilePath = [pngFilePath(1:end-4) '.pdf'];
end

img = imread(pngFilePath);
fh = figure('Visible', 'off');
imshow(img, 'border', 'tight');
export_fig(pdfFilePath, '-nocrop');
close(fh);
