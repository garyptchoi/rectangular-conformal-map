function L = cotangent_laplacian(v,f)
% Compute the cotagent Laplacian of a mesh.
%
% If you use this code in your own work, please cite the following paper:
% [1] T. W. Meng, G. P.-T. Choi and L. M. Lui, 
%     "TEMPO: Feature-Endowed Teichm�ller Extremal Mappings of Point Clouds."
%     SIAM Journal on Imaging Sciences, 9(4), pp. 1922-1962, 2016.
%
% Copyright (c) 2015-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);

f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

l1 = sqrt(sum((v(f2,:) - v(f3,:)).^2,2));
l2 = sqrt(sum((v(f3,:) - v(f1,:)).^2,2));
l3 = sqrt(sum((v(f1,:) - v(f2,:)).^2,2));

s = (l1 + l2 + l3)*0.5;
area = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
 
cot12 = (l1.^2 + l2.^2 - l3.^2)./area/2;
cot23 = (l2.^2 + l3.^2 - l1.^2)./area/2; 
cot31 = (l1.^2 + l3.^2 - l2.^2)./area/2; 
diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;

II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3];
V = [cot12; cot12; cot23; cot23; cot31; cot31; diag1; diag2; diag3];
L = sparse(II,JJ,V,nv,nv);

end