classdef TestSuite < matlab.unittest.TestCase
    methods (Test)
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
        function N_Q4Tests(testCase)
            % Loop through all the corners:
            % (-1,-1,-1),(1,-1,-1),(-1,1,-1), etc.
            % evaluate the shape functions and check if they yield one.
            magnitude = zeros(8,1);
            value = zeros(8,1);
            count = 1;
            for xi = -1:2:1
                for eta = -1:2:1
                    for mu = -1:2:1
                        results = Element.N_Q4(xi,eta);
                        value(count) = any(results==1);
                        magnitude(count) = norm(results);
                        count = count + 1;
                    end
                end
            end
            % We check that there is only one 1 in the series
            testCase.verifyEqual(true,all(magnitude==1));
            testCase.verifyEqual(true,all(value));
        end
        function ShellQ4_Test(testCase)
            % Since the element is 'cubic'
            a = 1;
            % we expect a diagonal jacobian with a 1/2 value
            [coords, ~, node_normals] = Factory.Shell([1,1],[a,a,a]);
            ele = Element(coords,node_normals);
            results = zeros(8,3);
            is_diag = zeros(8,1);
            count = 1;
            for xi = -1:2:1
                for eta = -1:2:1
                    for mu = -1:2:1
                        jac = ele.jacobian(xi,eta,mu);
                        results(count,:) = diag(jac);
                        is_diag(count) = isequal(diag(diag(jac)),jac);
                        count = count + 1;
                    end
                end
            end
            % Are all the diagonal values == 1/2?
            testCase.verifyEqual(true,all(all(results==1/2)));
            % Were all the jacobians diagonal?
            testCase.verifyEqual(true,all(is_diag));
        end
    end
end