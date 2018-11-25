classdef WriteHeader < handle
    properties
       fpath_header_name;
       fpath_cglobal_name;
       fname;
       fid_h;
       fid_c;
    end
    
    methods
        function self = WriteHeader(fname)
            fpath_parts = strsplit(fname, '/');
            self.fname = fpath_parts{end};
            self.fpath_header_name = sprintf('%s.h', fname);
            self.fpath_cglobal_name = sprintf('%s.c', fname);
            self.fid_h = 1;
            self.fid_c = 1;
        end
        
        function  open(self)
            self.fid_h = fopen(self.fpath_header_name, 'w+');
            self.fid_c = fopen(self.fpath_cglobal_name, 'w+');
        end
        
        function stat = close(self)
           stat = fclose(self.fid_h);
           stat = fclose(self.fid_c);
        end
        function stat = write_ifndef(self)
            fparts = strsplit(self.fname, '.');
            
            defname = sprintf('_%s_', fparts{1});
            
            ifdef_str = sprintf('#ifndef %s\n#define %s\n', defname,...
                defname);
            fprintf(self.fid_h, '%s', ifdef_str);
        end
        
        function stat = close_ifdef(self)
            fprintf(self.fid_h, '#endif\n');
        end

        function write_test_data(self, var, var_name_c)
        % write_test_data(self, var, var_name_c)
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
        % Use map container to translate to c data-types.
        %
        %   matlab types
            keySet  = {'double',
                       'int8',
                       'int16',
                       'int32'};
          % c-types
          % This could get a lot fancier, to account for different c-types,
          % long, unsigned etc.
          valSet = { {'double', '%0.30g'};
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
                fmt_str_c = sprintf('%s %s = %s;\n', ctype{1}, var_name_c, ctype{2});
                fmt_str_h = sprintf('%s %s;\n', ctype{1}, var_name_c);
                fprintf(self.fid_c, fmt_str_c, var);
                fprintf(self.fid_h, fmt_str_h);
                
            else % account for arrays
                if min(size(var)) > 1
                    error('we cant yet handle multidimensional arrays. Consider using var(:)');
                end

                scalar_fmt_c = sprintf('%s, ', ctype{2});
                scalar_fmt_c_end = sprintf('%s', ctype{2});
                array_fmt_c = repmat(scalar_fmt_c, 1, length(var(:))-1);

                array_fmt_c = [array_fmt_c, scalar_fmt_c_end];
                fmt_str_c = sprintf('%s %s[] = {%s};\n', ctype{1}, var_name_c, array_fmt_c);
                fmt_str_h = sprintf('%s %s[%d];\n', ctype{1}, var_name_c, length(var));
                
                fprintf(self.fid_c, fmt_str_c, var);
                fprintf(self.fid_h, fmt_str_h);

                len_name = sprintf('LEN_%s', var_name_c);
                fprintf(self.fid_c, 'int %s = %d;\n', len_name, length(var(:)));
                fprintf(self.fid_h, 'int %s;\n', len_name);
                

            end


        end
    end
end