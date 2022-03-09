function output_mat = getProperty( input_array, property_name )
% This function extracts the property named <property_name> from each
% element in the array <input_array> and compiles them into a matrix, where
% each row is the property for that element. If the property being
% retrieved is matrix-valued, then this matrix is reshaped into a vector.

output_mat = cell2mat( cellfun( @(x) x.(property_name)(:)', input_array(:), 'UniformOutput', false ) );

end