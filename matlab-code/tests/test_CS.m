classdef test_CS < matlab.unittest.TestCase
% Tests for the CS functions and classes. Run with
%>>> run(test_CS)

    methods (Test)
        function test_MeasEntityMu(testCase)
            x0 = 0.3;
            y0 = -0.2;
            N = 4;
            ind = 1;
            xy_rate = [0.1, 0.2];
            me_mu = MeasEntityMu(x0, y0, N, ind, xy_rate);
            expected = [x0; y0; ind;
                     x0 + xy_rate(1)*1; y0 + xy_rate(2)*1; ind;
                     x0 + xy_rate(1)*2; y0 + xy_rate(2)*2; ind;
                     x0+xy_rate(1)*3; y0+xy_rate(2)*3;      -1];
            
            testCase.assertEqual(expected, me_mu.as_vector(),...
                                'AbsTol', 1e-10)
            
        end
        
        function test_MeasEntityStatic(testCase)
            % Test for MeasEntityStatic
            xr = [1];
            yr = [2];
            N = 4;
            ind = 5
            me = MeasEntityStatic(xr, yr, N, ind);
            expected = [xr; yr; ind;
                        xr; yr; ind;
                        xr; yr; ind;
                        xr; yr; -1];
            testCase.verifyEqual(expected, me.as_vector)                                
        end

        function test_MoveEntityStatic(testCase)
            % Test for MoveEntityStatic
            xr = [1];
            yr = [2];
            N = 4;
            mve = MoveEntityStatic(xr, yr, N);
            expected = [xr; yr; 0;
                        xr; yr; 0;
                        xr; yr; 0;
                        xr; yr; 0];
            testCase.verifyEqual(expected, mve.as_vector)
        end
        
        function test_MasterTrajster(testCase)
            % test the master trajster class using statice ME and MVE.
            N_mve = 3;
            XR = [-.5, 0.25];
            YR = [-.4, 0.15];
            meta_cell = {2, 3};
            MT = MasterTrajster(XR, YR, meta_cell,...
                MoveEntityStatic.factory(N_mve), MeasEntityStatic.factory);

            % The last index of each ME should be -1.
            expected = [XR(1); YR(1); 0;  % move
                        XR(1); YR(1); 0;  % move
                        XR(1); YR(1); 0;  % move
                        XR(1); YR(1); 1;  % meas 1
                        XR(1); YR(1); -1;  % meas 1
                        XR(2); YR(2); 0;  % move 2
                        XR(2); YR(2); 0;  % move 2
                        XR(2); YR(2); 0;  % move 2
                        XR(2); YR(2); 2;  % meas 2
                        XR(2); YR(2); 2;  % meas 2
                        XR(2); YR(2); -1]; % meas 2
                        
            master_data_vec = MT.as_vector(); 
            
            testCase.verifyEqual(expected, MT.as_vector)
        end
    end
    
end

