clc
clear all
addpath('tests');

%% Should commit new tests to the list after they are passing locally

ele_suite = ElementTest;
run(ele_suite)

mesh_suite = MeshTest;
run(mesh_suite)

femcase_suite = FemCaseTest;
run(femcase_suite)

piezo_suite = PiezoTest;
run(piezo_suite)

mass_suite = DynamicTest;
run(mass_suite)