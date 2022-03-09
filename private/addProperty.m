function array = addProperty(array, propertyName, values)
% This function adds the property 'propertyName' to the input array of
% structures, 'array', using the values given. If values are given as a 
% matrix, each row will be added to each element of the array. If 'values'
% is a cell array, each element is added as the value to the specified
% field in each element of the input array.

% Check if a matrix was given, loop over each row if so
if isa(values,'double')
    
    for k = 1:length(array)
        array{k}.(propertyName) = values(k,:);
    end
    
% If not a matrix that was given, add each element of the values array as
% a propety of the input array
else
    
    for k = 1:length(values)
        array{k}.(propertyName) = values{k};
    end
    
end