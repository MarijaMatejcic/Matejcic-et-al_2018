function [pos, feat_name, FileName]= PSE_MitoNuc_loadNucData(PathName,FileName)

%% description:
% helper function for the script PSE_MitoNuc_mitosisDensityEstimation.m
% PahtName: path to the position file
% FileName: name of the position file
%  the file read are are orginaly obtained from Imaris. The script expect an ascii file
%  with 3 header lines, one line with comma separated column title and  
%  n lines  with 3 coordinates values separated by a comma.
%
% the script returns
% pos : an nx3 array of the position read in the file
% feat_name : a string of the column title
% FileName : the name of the opened file

% Authors:
% Marija Matejcic, Norden lab, MPI-CBG, Dresden, Germany
% Benoit Lombardot, Scientific Computing Facility, MPI-CBG, Dresden, Germany
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% License: BSD 3, see the LICENSE.txt file coming with the script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load nuclei data
if isempty(FileName) || isempty(PathName) 
    [FileName,PathName,FilterIndex] = uigetfile('/Volumes/IKNM_TissueMechanics/PSE_neurogenesis/01-Analysis/pH3_detection_n_OC_size/02-NeuroDiv_Positions/*.csv','Choose mitotic nuclei positions file','');
end
% display ('loading data script runs!')

fnChunks = strsplit(FileName, '_');
f = fopen([PathName FileName]);
dump = fgetl(f);
dump = fgetl(f);
dump = fgetl(f);
feat_name = fgetl(f);

% generating array pos, containing xyz coords extracted from file 
pos = [];
line = fgetl(f);
while ischar(line)
    chunks = strsplit(line,',');
    % array nx3 with x, y, z coordinates
    aux = [str2num(chunks{1}),str2num(chunks{2}),str2num(chunks{3})];
    % concatenation of arrays pos and aux along their dimension 1 = filling up empty pos array
    pos = cat(1,pos, aux);
    line = fgetl(f);
end
fclose(f);

