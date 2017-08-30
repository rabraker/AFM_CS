function [ hp ] = plotbyindex(ax, t, x, ind_vec, indc)

    hp = [];

    ind_cur = 1;
         stop = 1;
     while stop
         ind_next = find(ind_vec(ind_cur:end) ~= ind_vec(ind_cur), 1, 'first')+ind_cur;
         if isempty(ind_next)
               stop = 0;
               ind_next = length(ind_vec);
         end
         
         t_i = t(ind_cur:ind_next-1);
         x_i = x(ind_cur:ind_next-1);
         try
         if ind_vec(ind_cur) < 0
             col_i = indc{1,abs(ind_vec(ind_cur))};
         else
             col_i = indc{1,end};
         end
         catch
            keyboard 
         end
         
         
         h = plot(ax, t_i, x_i, 'color', col_i, 'LineWidth', 1.25);
         hold on
         hp = [hp; h];
         ind_cur = ind_next;
     end




end

function kjump = jumpindex(indeces)

kjump = find(diff(indeces) ~=1)
% for k = 1:length(indeces) - 1
%     if indeces(k) ~= indeces(k+1) -1
%         kjump =[kjump; k];
%     end
% end



end

