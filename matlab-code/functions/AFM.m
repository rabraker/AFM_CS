classdef AFM
  % Static methods that provide constants associated with the nPoint stage and
  % z-axis parameters.
  
  methods (Static)
    
    function mic = volts2mic_xy()
      % mic = volts2mic_xy()
      % Unit conversion of volts 2 microns for xy-stage
      mic = 5;
    end
    
    function volt = mic2volt_xy()
      % volt = mic2volt_xy()
      % Unit conversion of microns to volts for xy-stage
      volt = 1/AFM.volts2mic_xy();
    end

    function mic = volts2mic_z()
      % mic = volts2mic_z()
      % Unit conversion of volts 2 microns for z-axis
      mic = 7/20; % 7 micron range in +-20 volts
    end

    function volts = mic2volts_z()
      % volts = mic2volts_z()
      % Unit conversion of volts 2 microns for z-axis
      volts = 1/AFM.volts2mic_z(); 
    end

    function nm = volts2nm_z()
      % nm = volts2nm_z()
      % Unit conversion of volts 2 nanometers for z-axis
      nm= AFM.volts2mic_z()*1000; 
    end

    function volts = nm2volts_z()
    % volts = nm2volts_z()
    % Unit conversion of volts 2 microns for z-axis
      volts = 1/AFM.volts2nm_z(); 
    end
    
  end
  
end