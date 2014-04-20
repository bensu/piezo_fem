classdef MeshTestSuite < matlab.unittest.TestCase
    methods (Test)
         function Assembly_Test(testCase)
            [coords, connect, node_normals] = Factory.Shell([1,1],[1,1,0.1]);
            mesh = Mesh(coords,connect,node_normals);
            material = Material(1,0.3,1);
            dofs_per_ele = 0;
            dofs_per_node = 5;
            f = @(ele) Physics.K(ele,material,2);
            S = mesh.assembly(dofs_per_node,dofs_per_ele,f);
            K = f(mesh.ele(1));
            testCase.verifyEqual(K,S);
        end
        function ShellMeshTest(testCase)
            [coords, connect, node_normals] = Factory.Shell([1,1],[1,1,0.1]);
            mesh = Mesh(coords,connect,node_normals);
            % Check basic properties
            testCase.verifyEqual(4,mesh.n_nodes);
            testCase.verifyEqual(1,mesh.n_ele);
            testCase.verifyEqual(4,mesh.nodes_per_ele);
            % Check Element creation
            ele = mesh.ele(1);
            testCase.verifyEqual(ele.coords,mesh.coords);
            testCase.verifyEqual(ele.normals,mesh.normals);
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