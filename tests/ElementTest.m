classdef ElementTest < matlab.unittest.TestCase
    methods (Test)
        function K_Mech(testCase)
            a = 1;
            [coords, ~, node_normals] = Factory.Shell([1,1],[a,a,a]);
            ele = Element(coords,node_normals);
            material = Material(1,0.3,1);
            K = Physics.K_Shell(ele,material,2);
            is_symmetric = @(x) all(all(abs(x-x.')<1e-5));
            testCase.verifyEqual(true,is_symmetric(K));
            % There should be a decoupling of some dofs for a perfect
            % element
            testCase.verifyEqual(zeros(4),K(1:5:end,3:5:end));
        end
        function B_Mech(testCase)
            a = 1;
            [coords, ~, node_normals] = Factory.Shell([1,1],[a,a,a]);
            ele = Element(coords,node_normals);
            B = Physics.B(ele,0,0,0);
            % B should be [6 x 5*nodes]
            testCase.verifyEqual([6,20],size(B));
            % in a cubic element there should be no coupling between the
            % angles + the z direction and the flat x,y components
            w_coupling = B(1:3,sort([3:5:20 4:5:20 5:5:20]));
            testCase.verifyEqual(zeros(3,12),w_coupling);
        end
        function Ndevs_Mech_ShellQ4Test(testCase)
            a = 1;
            [coords, ~, node_normals] = Factory.Shell([1,1],[a,a,a]);
            ele = Element(coords,node_normals);
            % Ndevs should be [9x20]
            right_size = true;
            % The sum of the first part should be 0 because its pure shape
            % functions
            begins_zero = true;
            % Since the mesh is perfectly square, the last row is always
            % zero
            last_row_zero = true;
            count = 1;
            for xi = [-1,1]
                for eta = [-1,1]
                    for mu = [-1,1]
                        devs = Physics.Ndevs(ele,xi,eta,mu);
                        right_size = right_size && ...
                            all(([9, 20] == size(devs)));
                        last_row_zero = last_row_zero && ...
                            all(zeros(1,20)==devs(9,:));
                        % Sum the matrixs as Cook 360 12.5-9
                        aux = zeros(9,5);
                        for n = 1:4
                            aux = aux + devs(:,index_range(5,n));
                        end
                        begins_zero = begins_zero && ...
                            all(all(zeros(9,3)==aux(:,1:3)));
                        count = count + 1;
                    end
                end
            end
            testCase.verifyEqual(true,right_size);
            testCase.verifyEqual(true,begins_zero);
            testCase.verifyEqual(true,last_row_zero);
        end
        function Jacobian_ShellQ4Test(testCase)
            % Since the element is 'cubic'
            a = 1;
            % we expect a diagonal jacobian with a 1/2 value
            [coords, ~, node_normals] = Factory.Shell([1,1],[a,a,a]);
            ele = Element(coords,node_normals);
            results = zeros(8,3);
            is_diag = true;
            count = 1;
            for xi = [-1,1]
                for eta = [-1,1]
                    for mu = [-1,1]
                        jac = ele.jacobian(xi,eta,mu);
                        results(count,:) = diag(jac);
                        is_diag = is_diag && ...
                            isequal(diag(diag(jac)),jac);
                        count = count + 1;
                    end
                end
            end
            % Are all the diagonal values == 1/2?
            testCase.verifyEqual(true,all(all(results==1/2)));
            % Were all the jacobians diagonal?
            testCase.verifyEqual(true,is_diag);
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
    end
end