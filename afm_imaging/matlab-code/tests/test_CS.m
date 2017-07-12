classdef test_CS < matlab.unittest.TestCase
% Tests for the CS functions and classes. Run with
%>>> run(test_CS)

    methods (Test)
        function test_MVE_static(testCase)
            % Test for MoveEntityStatic
            xr = [1];
            yr = [2];
            N = 4;
            mve = MoveEntityStatic(xr, yr, N);
            expected = [xr; yr; 0;
                        xr; yr; 0;
                        xr; yr; 0;
                        xr; yr; 0];
            testCase.verifyEqual(expected, mve.asVector)
        end
        
        function test_ME_static(testCase)
            % Test for MeasEntityStatic
            xr = [1];
            yr = [2];
            N = 4;
            ind = 5
            me = MeasEntityStatic(xr, yr, N, ind);
            expected = [xr; yr; ind;
                        xr; yr; ind;
                        xr; yr; ind;
                        xr; yr; ind];
            testCase.verifyEqual(expected, me.asVector)                                
        end
        
        function test_MasterTrajster(testCase)
            % test the master trajster class using statice ME and MVE.
            N_mve = 3;
            XR = [-.5, 0.25];
            YR = [-.4, 0.15];
            meta_cell = {2, 3};
            MT = MasterTrajster(XR, YR, meta_cell,...
                MoveEntityStatic.factory(N_mve), MeasEntityStatic.factory);

            expected = [XR(1); YR(1); 0; % move
                        XR(1); YR(1); 0; % move
                        XR(1); YR(1); 0; % move
                        XR(1); YR(1); 1; % meas 1
                        XR(1); YR(1); 1; % meas 1
                        XR(2); YR(2); 0; % move 2
                        XR(2); YR(2); 0; % move 2
                        XR(2); YR(2); 0; % move 2
                        XR(2); YR(2); 2; % meas 2
                        XR(2); YR(2); 2; % meas 2
                        XR(2); YR(2); 2;] % meas 2
                        
            master_data_vec = MT.asVector(); 
            
            testCase.verifyEqual(expected, MT.asVector)
        end
    end
    
end

