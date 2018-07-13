function d = PSE_MitoNuc_getdistance(X1,X2,distType)

%% description:
% helper function for the script PSE_MitoNuc_mitosisDensityEstimation.m
% Measure the distance between two points.
% X1: an array indicating a position
% X2: an array indication a position
% the type of distance is indicated with a string set in the variable
% distType. It is either 'Euclidean' for the euclidean distance or
% 'GreatCircle' for the length of arc between to two points belonging to
% the same sphere
% the script returns the distance d between X1 and X2

% Authors:
% Marija Matejcic, Norden lab, MPI-CBG, Dresden, Germany
% Benoit Lombardot, Scientific Computing Facility, MPI-CBG, Dresden, Germany
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% License: BSD 3, see the LICENSE.txt file coming with the script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if strcmp(distType,'Euclidean')
    d = sum((X1-X2).^2).^0.5;
    
elseif strcmp(distType,'GreatCircle')
    % expect data belonging to the same sphere
    % X1 and X2 are vector from the center of the sphere to the point
    % if R1!=R2 the mean radius is used
    R1 = sum(X1.^2).^0.5;
    R2 = sum(X2.^2).^0.5;
    X1 = X1/R1;
    X2 = X2/R2;
    R = (R1+R2)/2;
    d = R*real(acos( sum( X1.* X2 ) ) );    
end
    