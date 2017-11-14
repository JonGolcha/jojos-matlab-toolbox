function full_mat = mat_sparse2full(sparse_vec,halfOnly)
% full_mat = mat_sparse2full(sparse_vec)
%
% ************************************************************************
% Convert a vector (sparse representation) to a symetric matrix (full 
% representation) 
%
% Input:
% ------
% sparse_vec    --> vector with only upper triangular part of full_mat Lx1
%                   (may be LxM dimensions, in which case a set of M
%                   matrices is returned)
% halfOnly      --> (optionnal) fills only half of full_mat
%                  	'1' or '-1 ' to fill only upper/lower triangular part
%
% Output:
% -------
% full_mat      --> symetric matrix 

% Example
% -------
%   sparse_vec = [1 2 3 4 5 6]
%   full_mat = mat_sparse2full(sparse_vec)
%
% History:
% --------
% *** 2016-04-08
% Created by J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com

if nargin<2 | isempty(halfOnly)
    halfOnly = 0;
end
if isrow(sparse_vec)
    sparse_vec = sparse_vec';
end

% from sparse to full representation
L = size(sparse_vec,1);
N = .5*(1+(1+8*L)^.5); % find back size of full_mat from length of sparse_vec (solution of quadratic equation)
if ~isvector(sparse_vec)
    M = size(sparse_vec,2);
    full_mat = zeros(N,N,M);
    for m = 1:M
        full_mat(:,:,m) = mat_sparse2full(sparse_vec(:,m),halfOnly);
    end
else
    full_mat = zeros(N,N);
    full_mat(find(tril(ones(N),-1))) = sparse_vec;
    if ~halfOnly
        full_mat = full_mat + full_mat';
    elseif halfOnly == -1
        full_mat = full_mat';
    end
end
