targetFunction = 'CostRailway'; % Nombre de la función a buscar
dirPath = pwd;  % Ruta del directorio con los archivos

% Lista todos los archivos en el directorio
fileList = dir(fullfile(dirPath, '*.m'));
fileNames = {fileList.name};

% Verificar dependencias para cada archivo
for i = 1:length(fileNames)
    filePath = fullfile(dirPath, fileNames{i});
    dependencies = matlab.codetools.requiredFilesAndProducts(filePath);
    if any(contains(dependencies, targetFunction))
        fprintf('La función "%s" es utilizada en: %s\n', targetFunction, fileNames{i});
    end
end