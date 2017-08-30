function [ s ] = metadata2latex(ExpMetaData, Ts)

    state_ticks = ExpMetaData.state_counts;
    state_times = state_ticks*Ts;
    time_total = sum(state_times);

    s_time = sprintf('xy-move & z-down & z-settle & xy-scan & z-up & total\\\\');
    s_dat  = sprintf('%.2f    &  %.2f  &   %.2f   &  %.2f   & %.2f &  %.2f\\\\',...
                        state_times,  sum(state_times));
    s= sprintf('%s\n%s', s_time, s_dat);
    s = sprintf('%s\n\\midrule\n', s);
    s = sprintf('%s $K_{i,z}$ & $\\delta_z$ [v/Ts] & $z_{ref} [v]$ & $z_{up} [v]$& -- &-- \\\\ \n', s);
    s = sprintf('%s    %.3f &     %.3g      &%.2f       & %.2f    & -- & --\\\\ \n',...
        s, ExpMetaData.Ki, ExpMetaData.ramp_rate, ExpMetaData.setpoint, ExpMetaData.z_UP);
    
    
    % s_row = '------------------------------------------------------------------';
    % for fld =fields(ExpMetaData)'
    %     if strcmp(fld{1},'state_counts')
    %         continue
    %     end
    %    s_i = sprintf('%s: %g', fld{1}, ExpMetaData.(fld{1}));
    %    if length(s_row) + length(s_i) < 79-5
    %        s_row = sprintf('%s  |  %s',s_row, s_i);
    %    else
    %        % Row is to long. Tack it onto the whole thing. 
    %        s = sprintf('%s\n%s', s, s_row);
    %        s_row = sprintf('%s', s_i);
    %    end

    % end
    
    

end

