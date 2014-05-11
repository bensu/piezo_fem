classdef MeshTest < matlab.unittest.TestCase
    methods (Test)
        function AssemblyTest(testCase)
            % Tests the assembly function using on 2 elements checking
            % where are the dofs coupled
            mesh = Factory.ShellMesh('Q4',[2,1],[1,1,0.1]);
            material = Material(1,0.3,1);
            dofs_per_ele = 0;
            dofs_per_node = 5;
            f = @(ele) Physics.K_Shell(ele,material,2);
            S = mesh.assembly_matrix(dofs_per_node,dofs_per_ele,f);
            % Nodes 4 and 6 should not have any coupling
            dofs_4 = index_range(dofs_per_node,4);
            dofs_6 = index_range(dofs_per_node,6);
            testCase.verifyEqual(zeros(dofs_per_node),S(dofs_4,dofs_6));
            % Nodes 1 and 3 should not have any coupling
            dofs_1 = index_range(dofs_per_node,1);
            dofs_3 = index_range(dofs_per_node,3);
            testCase.verifyEqual(zeros(dofs_per_node),S(dofs_1,dofs_3));
            % Nodes 1 and 4 should be coupled
            testCase.verifyEqual(true,any(any(S(dofs_1,dofs_4))));
            % Nodes 3 and 6 should be coupled
            testCase.verifyEqual(true,any(any(S(dofs_3,dofs_6))));
        end
        function AssemblyBasicTest(testCase)
            % Tests the assembly function using only one element
            mesh = Factory.ShellMesh('Q4',[1,1],[1,1,0.1]);
            material = Material(1,0.3,1);
            dofs_per_ele = 0;
            dofs_per_node = 5;
            f = @(ele) Physics.K_Shell(ele,material,2);
            S = mesh.assembly_matrix(dofs_per_node,dofs_per_ele,f);
            K = f(mesh.ele(1));
            testCase.verifyEqual(true,near(norm(diag(K)),norm(diag(S))));
        end
        function ShellMeshTest(testCase)
            mesh = Factory.ShellMesh('Q4',[1,1],[1,1,0.1]);
            % Check basic properties
            testCase.verifyEqual(4,mesh.n_nodes);
            testCase.verifyEqual(1,mesh.n_ele);
            testCase.verifyEqual(4,mesh.nodes_per_ele);
            % Check Element creation
            ele = mesh.ele(1);
            testCase.verifyEqual(true,near(ele.coords(1,1),mesh.coords(1,1)));
        end
        function ShellFactoryTest(testCase)
            [coords, connect, node_normals] = Factory.Shell([1,1],[1,1,0.1]);
            % Extract basic properties.
            [n_nodes, dim] = size(coords);
            [n_ele, nodes_per_ele] = size(connect);
            % Compare them to expected results.
            testCase.verifyEqual(4,n_nodes);
            testCase.verifyEqual(1,n_ele);
            testCase.verifyEqual(3,dim);
            testCase.verifyEqual(4,nodes_per_ele);
            testCase.verifyEqual(size(node_normals),size(coords));
        end
        function BoxFactoryTest(testCase)
            % Generate simplest box
            [coords, connect] = Factory.Box([1,1,1],[1,1,1]);
            % Extract basic properties.
            [n_nodes, dim] = size(coords);
            [n_ele, nodes_per_ele] = size(connect);
            % Compare them to expected results.
            testCase.verifyEqual(8,n_nodes);
            testCase.verifyEqual(1,n_ele);
            testCase.verifyEqual(3,dim);
            testCase.verifyEqual(8,nodes_per_ele);
        end
    end
end