README
PSE_MitoNuc tool
————————
————————

This file contains the following sections: 
Introduction, Input data, MATLAB Workflow

//Introduction: 
This MATLAB tool package is used in the analysis of mitotic cell distribution, to project 3D coordinates onto a circular 2D density heat map. 
It was built for the analysis of tissue-wide 3D images of immunostained zebrafish retinal samples. It might, however, be adapted to any nearly-spherical image for which specific 3D coordinates wish to be projected into 2D.


//Input data: 
1. 5 object-defining points (3D coordinates; see below)
2. List of 3D coordinates (see below)

1. 5 object-defining points
In order to use the “PSE MitoNuc” tool, the user should define five points on eachsample to build tissue object-defining vectors. For the retina, these were: 
1) center of the retinal cup (~widest point of the cup, close to the center of the lens), remaining 4 points were defined on this plane; 
2) outermost point of the optic fissure; 
3) point opposite to 2); 
4) nasal-most point;5) temporal-most point. 
It is essential here to correctly recognize nasal and temporal points. For this, it is useful to note in advance for each sample whether it is the left, or the right eye. 

X, y, z coordinates of these object-defining points and respective sample names should all be collected in a single matrix, formatted in columns as: 1x, 1y, 1z, 2x, 2y, 2z, 3x, …, to be read-in by the MATLAB script. 

In addition, the 3D coordinates of mitotic cells should be in distinct files for each sample, formatted as: 
1x, 1y, 1z
2x, 2y, 2z
…


//MATLAB Workflow:
 To run the “PSE MitoNuc” tool, the main script “PSE_MitoNuc_mitosisDensityEstimation.m” should be started. 
The value of the variable sphDataName should be changed to the path to the file that contains the object-defining points. When prompted, the user should also choose the file containing the mitotic cell coordinates. 

The script “PSE_MitoNuc_convertXYZ_to_RThPh.m” will then estimate thesymmetry axis of the tissue, as well as 2 perpendicular vectors, to define a newcoordinate system. 
X, y, z nuclear coordinates will be transformed into new Cartesian and spherical (r, theta, phi) coordinates, by projecting them onto a unit sphere defined in the new system. 
Euclidean distance will be calculated for each such mitotic cell to a point in the defined spherical mesh, and all these distances then used to calculate mitotic cell densities on each point of the mesh. The result will be plotted as a heatmap both on a sphere, and on a 2D projection. To generate this, a MATLAB “m_map” package and its “Azimuthal Equal-area” projection is used (see: http://www2.ocgy.ubc.ca/~rich/private/mapug.html).