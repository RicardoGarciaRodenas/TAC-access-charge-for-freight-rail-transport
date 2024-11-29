function exploreStruct(structure, prefix)
    % Recursively explore the fields of a structure, including nested substructures.
    % Arguments:
    % - structure: The main structure to explore.
    % - prefix: The prefix for the field names (used for recursion).
    
    if nargin < 2
        prefix = ''; % Empty prefix for the initial call
    end
    
    fields = fieldnames(structure); % Get the fields of the structure
    for i = 1:numel(fields)
        fieldName = fields{i};
        fullFieldName = [prefix, '.', fieldName]; % Construct the full field name
        fieldValue = structure.(fieldName);
        
        % Display the field name
        if isstruct(fieldValue)
            disp(['Structure found at: ', fullFieldName]);
            exploreStruct(fieldValue, fullFieldName); % Recursive call for nested structures
        else
            disp(['Field: ', fullFieldName, ' - Type: ', class(fieldValue)]);
        end
    end
end