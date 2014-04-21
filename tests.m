clc
clear all
addpath('tests');

% ele_suite = ElementTest;
% run(ele_suite)
% 
% mesh_suite = MeshTest;
% run(mesh_suite)

femcase_suite = FemCaseTest;
run(femcase_suite)