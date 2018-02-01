function deleteAndAppendPdfs(fileName, pdfs)
% pdfs is a cell matrix of pdf file names

if ~isempty(fileName)
    if exist(fileName, 'file') == 2
        delete(fileName);
    end
    append_pdfs(fileName, pdfs{:});
end
