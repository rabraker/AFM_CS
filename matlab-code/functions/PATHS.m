classdef PATHS
    %PATHS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [ PATH ] = exp()
            % PATH constant to where all experimental data is stored for the
            % MPC journal paper.
            if ispc
                PATH = 'Z:\afm-cs';
            else
                PATH = '/media/labserver/afm-cs';
            end
        end
        function PATH = sim_data()
          % Path to simulation data folder.
          PATH = fullfile(PATHS.MPCJ_root(), 'data');
        end
        function PATH = step_exp()
          PATH = fullfile(PATHS.exp, 'steps');
        end
        
        function [ PATH ] = CS_root()
            % PATH constant to where all experimental data is stored for the 
            % MPC journal paper.
            if ~ispc()
              PATH = fullfile(getMatPath, 'afm-cs');
            else
              PATH = 'C:\Users\arnold\Documents\afm-cs';
            end
        end
        
        function [ PATH ] = labview()
            % PATH constant to where all experimental data is stored for the 
            % MPC journal paper.

            PATH = fullfile(PATHS.CS_root, 'labview');
        end
        
        function path = reconstruction_BP()
          path = fullfile(PATHS.CS_root, 'reconstruction', 'BP');
        end
        function path = reconstruction_SMP1D()
          path = fullfile(PATHS.CS_root, 'reconstruction', 'SMP_1D');
        end
        
%         function [ PATH ] = sim_models()
%             % PATH constant to where all experimental data is stored for the
%             % MPC journal paper.
%             
%             PATH = fullfile(PATHS.MPCJ_root, 'models');
%         end
%         
%         
        function [ PATH ] = sysid()
            % PATH constant to where the system ID data is stored
            
            PATH = fullfile(PATHS.exp, 'sysID');
        end
        
        function PATH = jfig()
           PATH = fullfile(PATHS.MPCJ_root, 'latex', 'figures'); 
        end

    end
    
end

