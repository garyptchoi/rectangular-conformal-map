% Compute the rectangular conformal mapping using the fast method in [1].
%
% Usage:
% map = rectangular_conformal_map(v,f,corner)
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

addpath('mfile')

%% Example 1: Human Face

load('human_face.mat');
corner = [611;107;2304;14924];

% visualize the input mesh with the specified corners
plot_mesh(v,f); 
view([0 90]); hold on;
plot3(v(corner(1),1),v(corner(1),2),v(corner(1),3),'ro','MarkerFaceColor','r');
plot3(v(corner(2),1),v(corner(2),2),v(corner(2),3),'go','MarkerFaceColor','g');
plot3(v(corner(3),1),v(corner(3),2),v(corner(3),3),'bo','MarkerFaceColor','b');
plot3(v(corner(4),1),v(corner(4),2),v(corner(4),3),'yo','MarkerFaceColor','y');

% compute the rectangular conformal map
map = rectangular_conformal_map(v,f,corner);

% visualize the rectangular conformal map with the specified corners
plot_mesh(map,f); 
hold on;
plot(map(corner(1),1),map(corner(1),2),'ro','MarkerFaceColor','r');
plot(map(corner(2),1),map(corner(2),2),'go','MarkerFaceColor','g');
plot(map(corner(3),1),map(corner(3),2),'bo','MarkerFaceColor','b');
plot(map(corner(4),1),map(corner(4),2),'yo','MarkerFaceColor','y');

% evaluate the angle distortion
angle_distortion(v,f,map);

%% A simple application: checkerboard texture mapping using the rectangular parameterization 

% create a checkerboard texture on the rectangle
[x,y] = meshgrid(0:0.0025:1,0:0.0025:max(map(:,2)));
I = zeros(size(x,1),size(x,2),3);
I(:,:,1) = 0.55;
I(:,:,2) = 0.25;
I(:,:,3) = 0.2;
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if (mod(round(x(i,j)*16),2) == 1 && mod(round(y(i,j)*16),2) == 1) || (mod(round(x(i,j)*16),2) == 0 && mod(round(y(i,j)*16),2) == 0)
            I(i,j,1) = 1;
            I(i,j,2) = 0.85;
            I(i,j,3) = 0.65;
        end
    end
end
figure;
imshow(I);

% map the texture onto the surface using the parameterization result
% F1 = TriScatteredInterp(map,v(:,1),'natural'); % for older versions of MATLAB
% F2 = TriScatteredInterp(map,v(:,2),'natural'); % for older versions of MATLAB
% F3 = TriScatteredInterp(map,v(:,3),'natural'); % for older versions of MATLAB
F1 = scatteredInterpolant(map,v(:,1),'natural');
F2 = scatteredInterpolant(map,v(:,2),'natural');
F3 = scatteredInterpolant(map,v(:,3),'natural');
X = F1(x,y);
Y = F2(x,y);
Z = F3(x,y);

% visualize the texture mapping result
figure;
surf(X,Y,Z,'FaceColor','texturemap','Cdata',I);
shading interp
set(gcf,'color','w'); 
axis equal tight off
ax = gca; ax.Clipping = 'off';
view([0 90])

%% Example 2: Chinese Lion

load('chinese_lion.mat');
corner = [7013;6582;12652;5716];

% visualize the input mesh with the specified corners
% plot_mesh(v,f); 
% one may also include the third input in plot_mesh to visualize an additional quantity defined on vertices
% here mean_curv is a precomputed nv x 1 vector representing the mean curvature of the input mesh
plot_mesh(v,f,mean_curv); 
view([0 80]); hold on;
plot3(v(corner(1),1),v(corner(1),2),v(corner(1),3),'ro','MarkerFaceColor','r');
plot3(v(corner(2),1),v(corner(2),2),v(corner(2),3),'go','MarkerFaceColor','g');
plot3(v(corner(3),1),v(corner(3),2),v(corner(3),3),'bo','MarkerFaceColor','b');
plot3(v(corner(4),1),v(corner(4),2),v(corner(4),3),'yo','MarkerFaceColor','y');

% compute the rectangular conformal map
map = rectangular_conformal_map(v,f,corner);

% visualize the rectangular conformal map with the specified corners
% plot_mesh(map,f);
plot_mesh(map,f,mean_curv); 
hold on;
plot(map(corner(1),1),map(corner(1),2),'ro','MarkerFaceColor','r');
plot(map(corner(2),1),map(corner(2),2),'go','MarkerFaceColor','g');
plot(map(corner(3),1),map(corner(3),2),'bo','MarkerFaceColor','b');
plot(map(corner(4),1),map(corner(4),2),'yo','MarkerFaceColor','y');

% evaluate the angle distortion
angle_distortion(v,f,map);

%% Example 3: Human Brain (without specifying the corners)

load('human_brain.mat')

% visualize the input mesh
% plot_mesh(v,f);
plot_mesh(v,f,mean_curv); 
view([90 0]); 

% compute the rectangular conformal map without specifying the corners 
map = rectangular_conformal_map(v,f);

% visualize the rectangular conformal map 
% plot_mesh(map,f); 
plot_mesh(map,f,mean_curv); 

% evaluate the angle distortion
angle_distortion(v,f,map);
