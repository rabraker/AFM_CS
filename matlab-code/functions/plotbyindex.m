% [ hp ] = plotbyindex(ax, t, x, ind_vec, indc)
% Inputs
% ------
%   ax : axis handle
%   t  : time vector
%   x  : x data. (nothing necessarily to do with the x-axis)
%   ind_vec : a vector of the meta indexes. This is the 6th column of the
%   csv file.
%   indc  : a cell array that looks something like this:
%
%       indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
%             'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};
%        Each column corresponds to a state. The first rows are color. The
%        second row is legend names. 
%     
% Outputs
%  ------
%  hp : a vector of plot handles. 
% The goal of this function is to plot different parts of the cs trajector
% with different colors based on the state machine state it was generated
% from. We can do this by using the state machine index. This is
% complicated by the fact that we cant just slice into the vector. Ie,
% doing 
% ind = find(ind_vec==-2)
% x_down = x(ind)
% Doesnt work. Matlab will connect the dots between the first set of down
% and and the second. So that is what this mess of code is trying to solve.
% So many little corner cases here. 


function [ hp ] = plotbyindex(ax, t, x, ind_vec, indc)
  hold(ax, 'on')
    hp = []; % For the plot handles. 

    ind_cur = 1; 
    go = 1;
    while go
         % Find the first spot where the state_index changes.
         if ind_cur > length(ind_vec)
             go = 0;
             continue
         end
            
         ind_next_tmp = find(ind_vec(ind_cur:end) ~= ind_vec(ind_cur), 1, 'first');

         % At the end, ind_next_tmp this will be empty. 
         if isempty(ind_next_tmp)
               go = 0;
               ind_next = length(ind_vec);
         else
             ind_next = ind_next_tmp + ind_cur -1;
         end
         
         if ind_next > ind_cur
         
             t_i = t(ind_cur:ind_next-1);
             x_i = x(ind_cur:ind_next-1);
             try
                 if ind_vec(ind_cur) < 0
                     color_i = indc{1,abs(ind_vec(ind_cur))};
                     name_i  =  indc{2,abs(ind_vec(ind_cur))};
                 else
                     color_i = indc{1,end};
                     name_i  = indc{2,end};
                 end
             catch
                keyboard 
             end

             h = plot(ax, t_i, x_i, 'color', color_i, 'LineWidth', 1.25);
             h.DisplayName = name_i;

             hold on
             hp = [hp; h];
             ind_cur = ind_next+1;
         end 
     end
    
    
    
%      while go
%          % Find the first spot where the state_index changes.
%          ind_next_tmp = find(ind_vec(ind_cur:end) ~= ind_vec(ind_cur), 1, 'first');
%          % At the end, this will be empty. 
%          if isempty(ind_next_tmp)
%                go = 0;
%                ind_next = length(ind_vec);
%          else
%              ind_next = ind_next_tmp + ind_cur;
%          end
%          
%          if ind_next > ind_cur
%          
%              t_i = t(ind_cur:ind_next-1);
%              x_i = x(ind_cur:ind_next-1);
%              try
%                  if ind_vec(ind_cur) < 0
%                      color_i = indc{1,abs(ind_vec(ind_cur))};
%                      name_i  =  indc{2,abs(ind_vec(ind_cur))};
%                  else
%                      color_i = indc{1,end};
%                      name_i  = indc{2,end};
%                  end
%              catch
%                 keyboard 
%              end
% 
%              h = plot(ax, t_i, x_i, 'color', color_i, 'LineWidth', 1.25);
%              h.DisplayName = name_i;
% 
%              hold on
%              hp = [hp; h];
%              ind_cur = ind_next;
%          end 
%      end




end

function kjump = jumpindex(indeces)

kjump = find(diff(indeces) ~=1)
% for k = 1:length(indeces) - 1
%     if indeces(k) ~= indeces(k+1) -1
%         kjump =[kjump; k];
%     end
% end



end

