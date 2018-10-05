%
% reads in a parotid gland binary mesh file
% 05/10/2018 J.Rugis
%

fname = '4sim_out_N4_p3-p2-p4-1tet.bin';

fileID = fopen(fname);

vertices_count = int32(fread(fileID, 1, 'int32'));
vertices_data = transpose(single(fread(fileID, [4 vertices_count], 'single')));
vertices = vertices_data(:, 1:3);
vertices_dfa = vertices_data(:, 4);

surface_tris_count = int32(fread(fileID, 1, 'int32'));
surface_tris = transpose(int32(fread(fileID, [3 surface_tris_count], 'int32')));

tets_count = int32(fread(fileID, 1, 'int32'));
tets_data = transpose(int32(fread(fileID, [6 tets_count], 'int32')));
tets = tets_data(:, 1:4);
tets_dfa = typecast(tets_data(:, 5), 'single');
tets_dfb = typecast(tets_data(:, 6), 'single');

apical_count = int32(fread(fileID, 1, 'int32'));
apical_tris = transpose(int32(fread(fileID, [1 apical_count], 'int32')));

basal_count = int32(fread(fileID, 1, 'int32'));
basal_tris = transpose(int32(fread(fileID, [1 basal_count], 'int32')));

common_tris_count = int32(fread(fileID, 1, 'int32'));
common_tris = transpose(int32(fread(fileID, [3 common_tris_count], 'int32')));

fclose(fileID);
