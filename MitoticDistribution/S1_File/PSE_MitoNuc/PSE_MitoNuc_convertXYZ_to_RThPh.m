function [radii, th, ph, pos_refManual]= PSE_MitoNuc_convertXYZ_to_RThPh(pos, FileName,sphDataName)

%% description
% helper function for the script PSE_MitoNuc_mitosisDensityEstimation.m
% this script build transforms euclidean position (pos) in a referential
% oriented along the Pseudo Stratified Epithelium (PSE) anatomy. In particular the
% referential is centered on eye lens, its x axis is oriented along the
% PSE junction and it x and y axis defines a plane "closing" the cup formed
% by the PSE
%
% input:
%   pos : a nx3 list of positions
%   Filename : used as an identifier to fetch the right dataset information from
% the excel file
%   sphDataName : an excel file containing information to to build PSE
% referential of multiple datasets
%
% output:
%   radii: a nx1 array with the first spherical cordinate (radius) for each
%   of the position in the anatomy oriented referential
%   th: a nx1 array with the 2nd spherical cordinate (theta) for each
%   of the position in the anatomy oriented referential
%   ph: a nx1 array with the 3rd spherical cordinate (phi) for each
%   of the position in the anatomy oriented referential
%   pos_refManual: a nx3 array with the euclidean coordinates for each
%   of the position in the anatomy oriented referential


% Authors:
% Marija Matejcic, Norden lab, MPI-CBG, Dresden, Germany
% Benoit Lombardot, Scientific Computing Facility, MPI-CBG, Dresden, Germany
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% License: BSD 3, see the LICENSE.txt file coming with the script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%sphDataName = 'E:\project_data\Marija(Norden)\update matlab distance distribution\Object_defining_vectors.xlsx'

fnChunks = strsplit(FileName, '_');% used to identify the data corresponding tp the current dataset
sample = strjoin(fnChunks(1:5), '_');

[num, txt] = xlsread (sphDataName, fnChunks{3});
[row,col] = find(ismember(txt,sample)==1);
row = row-1;

A = [num(row,1:3)]; % center of the ocular sphere (=lens center)
B = [num(row,4:6)]; % first junction of PSE external ring
C = [num(row,7:9)]; % second junction of PSE external ring
D = [num(row,10:12)]; % point 1 in between the junctions 
E = [num(row,13:15)]; % point 2 in between the 2 junctions


%% definition of the direct orthonormal referential (u,v,n)
% n is roughly (because it is manually determined) the symetrie axis of the PSE
% u is pointing toward the tissue junction/stitch
% v completes the referential such that u,v,n is direct

% defintion of the symetry axis unit vector with the cross product of BC and DE
BC = C-B;
DE = E-D;
N = cross(BC,DE);
n = N/sum(N.^2).^0.5;

% defintion of the second and third unit vector of the referential (u,v)
%  (n,u,v) is a referential with direct orientation
u = B-A;
u = u-n*(sum(u.*n));
u = u/norm(u,2); % u is perpendicular to n and belong to the plane A,B,F
v = cross(n,u); 
v = v/norm(v,2); % v is perpendicular to both n and u

%% projection of nuclei position in the new referential
% a) spherical coordinates
% mesure of nuclei position in the spherical coordinate defined by the base u,v,n
% radius,th, ph 
% ncell = number of spots
ncell = numel(pos)/3;
ph= zeros(ncell,1);
th = zeros(ncell,1);
radii = zeros(ncell,1);
for i=1:ncell
    % phi : pi/2 - latitude
    AX = pos(i,:)-A;
    ax = AX/sum(AX.^2).^0.5;
    ph(i)= acos(sum(ax.*n)); % belongs to [0, pi] 
    
    % theta : longitude
    AXn = n*(sum(AX.*n));
    AXuv = AX-AXn;
    axuv = AXuv/norm(AXuv,2);
    cos_th = sum(axuv.*u);
    th(i) = acos(cos_th)*sign(sum(axuv.*v)); % belongs to [-pi pi]
    
    % d(center-spot)
    radii(i) = norm(AX,2);
end

% b) euclidean coordinates

pos_aux = pos - ones(ncell,1)*reshape(A,[1,3]);
pos_refManual = zeros(size(pos)); 
pos_refManual(:,1) = sum( pos_aux .* (ones(ncell,1)*reshape(u,[1,3])) , 2 );
pos_refManual(:,2) = sum( pos_aux .* (ones(ncell,1)*reshape(v,[1,3])) , 2 );
pos_refManual(:,3) = sum( pos_aux .* (ones(ncell,1)*reshape(n,[1,3])) , 2 );



