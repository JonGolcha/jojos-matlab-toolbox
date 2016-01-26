function linear_index = subv2ind(siz, subscript_vector)
%SUBV2IND return linear indexes into an ND-array given by a vector positions
%
% linear_index = subv2ind(siz, subscript_vectors)
%
% Inputs:
%                   siz Dx1 or 1xD
%      subscript_vector    NxD       each row specifies an element in an array
%                                    with dimensions siz
%
% Outputs:
%         linear_index     Nx1       linear indexes specifying the same elements
%
% Matlab's sub2ind returns the scalar linear indexes corresponding to an
% argument list of subscripts giving the rows, cols, 3rd-dims, ... of each
% element. In code that deals with multi-dimensional arrays, it is common to
% have the subscripts stored in a rows of an array. This function will take the
% subscripts in a single array and give the corresponding linear indexes for a
% given array size, siz.
%
% See also: ind2subv, ind2sub, sub2ind

% Iain Murray, November 2007, March 2011

if (size(subscript_vector, 2) ~= numel(siz))
    error('subscripts are wrong length for stated array dimensions.');
end

subscripts = num2cell(subscript_vector, 1);
linear_index = sub2ind(siz(:)', subscripts{:});
