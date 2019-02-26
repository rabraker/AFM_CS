% [ s ] = metadata2text(ExpMetaData, Ts)
%
% Converts the metadata saved by labview as a .mat file into a string that
% can be wither displayed or printed onto a figure. This is specific to
% play-cs-imaging.vi

function [ s ] = metadata2text(ExpMetaData, Ts)

    state_ticks = ExpMetaData.state_counts;
    state_times = state_ticks*Ts;
    time_total = sum(state_times);

    s_time = sprintf('xy-move | z-down | z-settle | xy scan | z-up ||| total');
    s_dat  = sprintf('%.2f      %.2f      %.2f      %.2f    %.2f     %.2f    %g',...
                        state_times,  sum(state_times));
    s= sprintf('%s\n%s', s_time, s_dat);
    s_row = '------------------------------------------------------------------';
%     s_cell = {s_row}
%     s_row = '';
    for fld =fields(ExpMetaData)'
        if strcmp(fld{1},'state_counts')
            continue
        end
        %
        %        s_i = sprintf('%s: %g', fld{1}, ExpMetaData.(fld{1}));
        %        if length(s_row) + length(s_i) < 79-5
        %            s_row = sprintf('%s  |  %s',s_row, s_i);
        %        else
        %            % Row is to long. Tack it onto the whole thing.
        %            s = sprintf('%s\n%s', s, s_row);
        %            s_row = sprintf('%s', s_i);
        %        end
        
        try
          fld_str = fld2str(ExpMetaData.(fld{1}));
          
          s_i = sprintf('%s: %s', fld{1}, fld_str);
          
          if length(s_row) + length(s_i) < 79-5
            s_row = sprintf('%s  |  %s',s_row, s_i);
          else
            % Row is to long. Tack it onto the whole thing.
            s = sprintf('%s\n%s', s, s_row);
            s_row = sprintf('%s', s_i);
          end
          %        if strcmp(fld{1}, 'Den') | strcmp(fld{1}, 'Num')
          %            keyboard
          %        end
          
          
        catch
          fprintf('skipping field %s\n', fld{1});
        end
    end
    s = sprintf('%s\n%s', s, s_row);

end

% Process arrays. Assume they are 1 x n, and trailing zeros are not of
% interest. 
function fld_str = fld2str(fld)

    if isa(fld, 'numeric')
        if length(fld) >1
           k = find(fliplr(fld) ==0, 1, 'last');
        else
            k = 1;
        end
       fld_str =  num2str(fld(1:end-k+1));
    else
        fld_str = sprintf('%s', fld);
    end

end