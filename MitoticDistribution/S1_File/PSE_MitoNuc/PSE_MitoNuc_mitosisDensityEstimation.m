%% description:
% The script builds a density map of mitosis in the pseudostratified
% epithelium (PSE). the mitosis position are expressed in a referential 
% adapted to the anatomy of the retinal PSE and project on a unit sphere. the
% density map is build from that projection by associating a gaussian
% kernel to each position and accumulating the projection of each mitosis
% at every position of the sphere.
% 
% the script produces two figures:
%  - a spherical map that can be interactively rotated:
%   the mitosisposition are shown with a with cross the sphere equator with a black dashed line 
%   and the first meridian with a dashed white line. the density map is show
%   with matlab built-in 'pink' colormap
%  - a projection of the north emisphere of the density map realized with
%   the m_map package from Rich Pawlowicz ( https://www.eoas.ubc.ca/~rich/map.html )
%   in particular we used the 'Azimuthal Equal-area' projection.
%
% usage:
% Before starting the script:
%   - copy the PSE_MitoNuc  script to your disc. There should be 4 of them.
%   this script is the master script whereas the 3 others are helper to
%   load data and do a few calculations
%   - copy the m_map package to you disc. it can be found at https://www.eoas.ubc.ca/~rich/map.html
%   - Set the paremeters in the parameters section below
%       - pseMitoNucFolder : folder where the PSE_MitoNuc scripts can be found
%       - m_map_Folder : folder where the m_map package was copied
%       - PathName : folder containing a file with a list of positions
%       if set to and empty string (i.e. '') the ui will popup to request
%       the file location
%       - FileName : name of a comma separated file with a of positions.
%       the file is expected to 3 lines header, one line with the 3 comma
%       separated column title and n lines with the comma separated value
%       of position coordinates (the format correspond to the spot position
%       we exported from Imaris).
%       if set to and empty string (i.e. '') the ui will popup to request
%       the file location
%       - sphDataName: the complete path to an excel file containing the
%       information to build a anatomy orientation referential of the PSE.
%       the excel file is organized in multiple sheet named after each
%       dataset name. this dataset name is deduced from the FileName
%       parameter. inside each spreadsheet the data are organize on one line.
%           - 3 first values are the coordinates of the center of the ocular sphere (=lens center)
%           - value  4 to  6 are the coordinates of the first junction of PSE external ring
%           - value  7 to  9 are the coordinates of the second junction of PSE external ring
%           - value 10 to 12 are the coordinates of the first point in between the junctions 
%           - value 13 to 15 are the coordinates of the second point in between the junctions 
% 
%   once these values are set you can run the script that will output the 2
%   density map representation describe above (sphere and north emisphere projection).


% Authors:
% Marija Matejcic, Norden lab, MPI-CBG, Dresden, Germany
% Benoit Lombardot, Scientific Computing Facility, MPI-CBG, Dresden, Germany
% the script was written between March and September 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% License: BSD 3, see the LICENSE.txt file coming with the script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% script path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pseMitoNucFolder = 'E:\centos_home_copy_2014-09-02\lombardo\Workspace\matlab\Marija\PSE_MitoNuc';
m_map_Folder = 'E:\centos_home_copy_2014-09-02\lombardo\Workspace\matlab\Marija\m_map';

% data location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PathName='';
FileName='';
sphDataName = '/Users/matejcic/Dropbox/OpticCup/Revisions/ActionPoints/Supporting/data_repository/MitoticDistribution/Object_defining_vectors.xlsx';



%% import
addpath(pseMitoNucFolder);
addpath(m_map_Folder);


%% processing
% load nuclei position
[pos, feat_name, FileName]= PSE_MitoNuc_loadNucData(PathName,FileName);
n_cell = numel(pos)/3;

% project cell position to a referential representative of the PSE
[r, th, ph, pos_centered] = PSE_MitoNuc_convertXYZ_to_RThPh(pos, FileName,sphDataName);

% projection of position on the unit sphere
pos_centered2 = pos_centered./(reshape(r, [numel(r) 1])* ones(1,3));

%% density estimation on a sphere
% define a sphere mesh
[Xs,Ys,Zs] = sphere(60);
sphPos = [Xs(:),Ys(:),Zs(:)];

% calculate density for each point on the sphere mesh
h=10*pi/180; % this is the standard deviation of the gaussian used for the density estimation (angle of 10 degree was used)

n_sphPos = size(sphPos,1);
dens = zeros(n_sphPos,1);
for i=1:n_sphPos
    for j = 1:n_cell
        d = PSE_MitoNuc_getdistance(sphPos(i,:),pos_centered2(j,:),'Euclidean');
        % d = PSE_MitoNuc_getdistance(X1,X2,'GreatCircle');
        dens(i) = dens(i) + exp(-(d/h)^2/2);
    end
end
%dens = dens/(sum(dens));
dens = reshape(dens,size(Xs));



%% display 

figure
hold on
% plot the density map
surf(Xs,Ys,Zs,dens,'EdgeColor','None','FaceColor','interp');

%plot the original points
plot3(pos_centered2(:,1),pos_centered2(:,2),pos_centered2(:,3),'w+')

% plot equator
theta_eq = (0:360) * pi/180;
x_eq = cos(theta_eq);
y_eq = sin(theta_eq);
z_eq = 0*theta_eq;
plot3(x_eq,y_eq,z_eq,'--k','linewidth',3)


% plot zeros meridian
phi_0 = (0:180)*pi/180;
x_mer0 = cos(0)*sin(phi_0);
y_mer0 = sin(0)*sin(phi_0);
z_mer0 = cos(phi_0);
plot3(x_mer0,y_mer0,z_mer0,'--w','linewidth',3)
colorbar
colormap(pink)
caxis([0, 20])

axis vis3d
grid on
% axis off



%% project the data with m_map package

% for documentation on the m_map package see:
% http://www2.ocgy.ubc.ca/~rich/private/mapug.html

% set a projection type
m_proj('Azimuthal Equal-area','longitudes',0,'latitudes',90,'radius',170);
% m_coast; m_grid;

% convert (th,ph) position to map coordinate and plot them
%  long = th*180/pi;
%  lat = 90-ph*180/pi;
[long,lat,r] = cart2sph(pos_centered2(:,1),pos_centered2(:,2),pos_centered2(:,3));
long = long*180/pi;
lat  = lat*180/pi;
% sel = lat>-pi;
[xNuc_map, yNuc_map] = m_ll2xy( long, lat );

% convert surf to patch
[f,v,c] = surf2patch(Xs(2:end-1,:),Ys(2:end-1,:),Zs(2:end-1,:),dens(2:end-1,:));
Xp = reshape( v(f,1) , size(f) );
Yp = reshape( v(f,2) , size(f) );
Zp = reshape( v(f,3) , size(f) );
[longP,latP,r] = cart2sph(Xp,Yp,Zp);
longP = longP*180/pi;
latP  = latP*180/pi;

[Xdens_map, Ydens_map] =m_ll2xy( longP, latP );
Cp = reshape( c(f)   , size(f) );

% display the projected data and map
figure
hold on
%m_plot(long,lat,'+w');
% m_patch(longP, latP, Cp,'EdgeColor','None','facecolor','flat')
% plot(yNuc_map, -xNuc_map,'+k')
plot(xNuc_map, yNuc_map,'+w')
colorbar
colormap(pink)
patch(Xdens_map', Ydens_map', Cp','EdgeColor','None');
m_grid;
caxis([0, 20])

   

% figure
% polar(th, ph*180/pi,'+r')
clear all
clc
   
   
   
   
   
   
   
   
   