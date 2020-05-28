function map = tutte_map(v,f,bdy_index,bdy)
% Compute the Tutte map with given boundary condition.
% Invoked only if the harmonic map fails due to very bad triangulations.
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% bdy_index: nb x 1 boundary vertex indices
% bdy: nb x 1 complex coordinates of the target boundary position
% 
% Output:
% map: nv x 2 vertex coordinates of the disk Tutte map
% 
% If you use this code in your own work, please cite the following paper:
% [1] T. W. Meng, G. P.-T. Choi and L. M. Lui, 
%     "TEMPO: Feature-Endowed Teichmüller Extremal Mappings of Point Clouds."
%     SIAM Journal on Imaging Sciences, 9(4), pp. 1922-1962, 2016.
%
% Copyright (c) 2015-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nf = length(f);

% Construct the Tutte Laplacian
I = reshape(f',nf*3,1);
J = reshape(f(:,[2 3 1])',nf*3,1);
V = ones(nf*3,1)/2;
W = sparse([I;J],[J;I],[V;V]);
M = W + sparse(1:nv,1:nv,-diag(W)-(sum(W)'), nv, nv);

[mrow,mcol,mval] = find(M(bdy_index,:));
M = M - sparse(bdy_index(mrow), mcol, mval, nv, nv) + ...
        sparse(bdy_index, bdy_index, ones(length(bdy_index),1), nv, nv);
    
b = zeros(nv,1);

b(bdy_index) = bdy;

z = M\b;
map = [real(z),imag(z)];

end