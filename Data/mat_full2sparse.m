function sparse_vec = mat_full2sparse(full_mat)
% sparse_vec = mat_full2sparse(full_mat)
%
% ************************************************************************
% Convert a symetric matrix (full) to a sparse representation (vector)
%
% Input:
% ------
% full_mat      --> symetric matrix (may be more 3 dimensions, in
%                   which case a set of vectors is returned)
%
% Output:
% -------
% sparse_vec    --> vector with only upper triangular part of full_mat
%
% Example
% -------
%   full_mat = reshape(1:16,[4 4])'
%   sparse_vec = mat_full2sparse(full_mat)
%
% History:
% --------
% *** 2016-04-08
% Created by J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com

% from full to sparse representation
N = size(full_mat,1);
if ismatrix(full_mat)
    sparse_vec = full_mat(tril(ones(N),-1)~=0);
elseif ndims(full_mat) == 3
    K = size(full_mat,3);
    sparse_vec = zeros(.5*N*(N-1),K);
    for k = 1:K
        sparse_vec(:,k) = mat_full2sparse(squeeze(full_mat(:,:,k)));
    end
    
else
   error('[mat_full2sparse] can only convert 2D or 3D matrices.'); 
end


