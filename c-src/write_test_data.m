% write_test_data(fid, var, var_name_c)
%
% This function is designed to write variables into a c file or c header
% file. It can currently accept scarlars and single dimension arrays, of
% type double, int8, int16, and int32. 
%
% Usage:
% var = 23.5
% var_name_c = 'var_ml';
% fid = fopen('somefile.h')
% write_test_data(fid, var, var_name_c);
%
% will produce a file with the content
% double var_ml = 23.5;
%
% If an array is passed, it will be set as
% type var_name_c[] = {val1, .... valn};
%
% The function will also create a variable prepended with LEN_var_name_c
% which is an int containing the length of the array.

function [ output_args ] = write_test_data(fid, var, var_name_c)
% Use map container to translate to c data-types.
%
%   matlab types
    keySet  = {'double',
               'int8',
               'int16',
               'int32'};
  % c-types
   % This could get a lot fancier, to account for different c-types, long, unsigned etc.       
    valSet = { {'double', '%0.15g'};
              {'int', '%d'};
              {'short int', '%d'};
              {'long int', '%d'}
              };

    map = containers.Map(keySet, valSet); % like a dict in python.

    if ~isKey(map, class(var))
        error('unrecognized type');
    end

    ctype = map(class(var));
    if isscalar(var)
        fmt_str = sprintf('%s %s = %s;\n', ctype{1}, var_name_c, ctype{2});
        fprintf(fid, fmt_str, var);
        
    else % account for arrays
        if min(size(var)) > 1
            error('we cant yet handle multidimensional arrays. Consider using var(:)');
        end
        
        scalar_fmt = sprintf('%s, ', ctype{2});
        scalar_fmt_end = sprintf('%s', ctype{2});
        array_fmt = repmat(scalar_fmt, 1, length(var(:))-1);
        
        array_fmt = [array_fmt, scalar_fmt_end];
        fmt_str = sprintf('%s %s[] = {%s};\n', ctype{1}, var_name_c, array_fmt);
        fprintf(fid, fmt_str, var);
        
        len_name = sprintf('LEN_%s', var_name_c);
        fprintf(fid, 'int %s = %d;\n', len_name, length(var(:)));

    end


end

