function [parent_name ] = get_parent_name(child_name, child_handle)


k = regexp(child_name, child_handle);
parent_name = sprintf('%s.csv', child_name(1:k-1));



end

