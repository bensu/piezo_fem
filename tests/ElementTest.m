classdef ElementTest < matlab.unittest.TestCase
    methods (Test)
        function K_Mech(testCase)
            a = 1;
            mesh = Factory.ShellMesh([1,1],[a,a,a]);
            ele = mesh.ele(1);
            material = Material(1,0.3,1);
            K = Physics.K_Shell(ele,material,2);
            is_symmetric = @(x) all(all(abs(x-x.')<1e-5));
            testCase.verifyEqual(true,is_symmetric(K));
        end
        function B_Mech(testCase)
            a = 1;
            mesh = Factory.ShellMesh([1,1],[a,a,a]);
            ele = mesh.ele(1);
            B = Physics.B_Shell(ele,0,0,0);
            % B should be [6 x 5*nodes]
            testCase.verifyEqual([6,20],size(B));
            % in a cubic element there should be no coupling between the
            % angles + the z direction and the flat x,y components
            w_coupling = B(1:3,sort([3:5:20 4:5:20 5:5:20]));
            testCase.verifyEqual(zeros(3,12),w_coupling);
        end
        function Jacobian_ShellQ4Test(testCase)
            % Since the element is 'cubic'
            % we expect a diagonal jacobian with a 1/2 value
            a = 1;
            mesh = Factory.ShellMesh([1,1],[a,a,a]);
            ele = mesh.ele(1);
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