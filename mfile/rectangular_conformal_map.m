function map = rectangular_conformal_map(v,f,corner)
% Compute the rectangular conformal mapping using the fast method in [1].
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% (optional) corner: 4 x 1 vertex indices for the four corners of the 
% rectangle, with anti-clockwise orientation
% 4 - 3
% |   |
% 1 - 2
% 
% Output:
% map: nv x 2 vertex coordinates of the rectangular conformal parameterization
% 
% Remarks:
% 1. Please make sure that the input mesh does not contain any 
%    unreferenced vertices/non-manifold vertices/non-manifold edges.
% 2. Please remove all valence 1 boundary vertices (i.e. vertices with 
%    only 1 face attached to them) before running the program.
% 3. Please make sure that the input triangulations f are with 
%    anti-clockwise orientation.
% 4. The output rectangular domain will always have width = 1, while the
%    height depends on the choice of the corners and may not be 1.
%    (The Riemann mapping theorem guarantees that there exists a conformal  
%    map from any simple-connected open surface to the unit square, but if 
%    four vertices on the surface boundary are specified to be the four 
%    corners of the planar domain, the theorem is no longer applicable.)
% 
% If you use this code in your own work, please cite the following paper:
% [1] T. W. Meng, G. P.-T. Choi and L. M. Lui, 
%     "TEMPO: Feature-Endowed Teichmüller Extremal Mappings of Point Clouds."
%     SIAM Journal on Imaging Sciences, 9(4), pp. 1922-1962, 2016.
%
% Copyright (c) 2015-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%% Initial setup

nv = length(v);

if size(v,1) < size(v,2)
    v = v';
end
if size(f,1) < size(f,2)
    f = f';
end
if size(v,2) == 2
    v = [v, zeros(nv,1)];
end

% TR = TriRep(f,v);     % for older versions of MATLAB
% B = freeBoundary(TR); % for older versions of MATLAB
TR = triangulation(f,v);
B = TR.freeBoundary;
bdy_index = B(:,1);
% Remark: The above approach for getting the surface boundary may not work
% well in case the boundary contains vertices with valence 1
% In that case, use some other method to obtain the boundary

if nargin < 3
    % if the corner indices are not specified, pick 4 boundary points as corner
    corner = bdy_index(1:ceil(length(bdy_index)/4):end);
end

if sum(ismember(corner,bdy_index)) ~= 4
    error('The input corner indices do not lie on the surface boundary!');
end

% rearrange the boundary indices
id = find(bdy_index == corner(1));
bdy_index = bdy_index([id:end,1:id-1]);

id1 = 1;
id2 = find(bdy_index == corner(2));
id3 = find(bdy_index == corner(3));
id4 = find(bdy_index == corner(4));

if id2 > id3
    % the boundary orientation is wrong, meaning the input f has wrong orientation
    % we correct the orientation and start over
    warning('The input triangulations are with clockwise orientation!');
    f = fliplr(f);
    % TR = TriRep(f,v);     % for older versions of MATLAB
    % B = freeBoundary(TR); % for older versions of MATLAB
    TR = triangulation(f,v);
    B = TR.freeBoundary;
    bdy_index = B(:,1);
    id = find(bdy_index == corner(1));
    bdy_index = bdy_index([id:end,1:id-1]);
    id1 = 1;
    id2 = find(bdy_index == corner(2));
    id3 = find(bdy_index == corner(3));
    id4 = find(bdy_index == corner(4));
end

%% Step 1: Mapping the input mesh onto the unit disk

bdy_length = sqrt((v(bdy_index,1) - v(bdy_index([2:end,1]),1)).^2 + ...
            (v(bdy_index,2) - v(bdy_index([2:end,1]),2)).^2 + ...
            (v(bdy_index,3) - v(bdy_index([2:end,1]),3)).^2);
partial_edge_sum = zeros(length(bdy_length),1);

% arc-length parameterization boundary constraint
for i = 2:length(bdy_length)
    for j = 1:i-1
        partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length);
bdy = exp(theta*1i);

% disk harmonic map
M = cotangent_laplacian(v,f);
[mrow,mcol,mval] = find(M(bdy_index,:));
M = M - sparse(bdy_index(mrow),mcol,mval,nv, nv) + ...
        sparse(bdy_index,bdy_index,ones(length(bdy_index),1),nv, nv);
c = zeros(nv,1); 
c(bdy_index) = bdy;
z = M\c;
disk = [real(z),imag(z)]; 

if sum(sum(isnan(disk))) ~= 0
    % use tutte embedding instead
    disk = tutte_map(v,f,bdy_index,bdy); 
end

%% Step 2: Mapping the unit disk to the unit square

% compute the generalized Laplacian
mu = beltrami_coefficient(disk,f,v);
Ax = generalized_laplacian(disk,f,mu); Ay = Ax;

% set the boundary constraints
bottom = bdy_index(id1:id2);
right = bdy_index(id2:id3);
top = bdy_index(id3:id4);
left = bdy_index([id4:end,id1]);

bx = zeros(nv,1); by = bx;
Ax([left;right],:) = 0;
Ax([left;right],[left;right]) = diag(ones(length([left;right]),1));
bx(right) = 1;
Ay([top;bottom],:) = 0;
Ay([top;bottom],[top;bottom]) = diag(ones(length([top;bottom]),1));
by(top) = 1;

% solve the generalized Laplace equation
square_x = Ax\bx;
square_y = Ay\by;

%% Step 3: Optimize the height of the square to achieve a conformal map

h_opt = fminbnd(@(h) sum(abs(beltrami_coefficient([square_x,h*square_y],f,v)).^2),0,5);
map = [square_x, h_opt*square_y];

end