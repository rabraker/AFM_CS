function [hands, legs] = plot_all_cycles(self, ax1, ax2, ax3, ax4, bounds, to_pix)
%   [hands] = plot_all_cycles(self, ax1, ax2, ax3, ax4, bounds, to_pix)
% Plot the color-ized cs cycles to the provided axes.
%
% Usage
% -----
% You can call this method in several ways
%
%  hands = plot_all_cycles(ax1) : plot uz to ax1
%  hands =  plot_all_cycles(..., ax2) : also plot ze to ax2
%  hands = plot_all_cycles(..., ax3) : also plot x to ax3
%  hands = plot_all_cycles(..., ax4) : also plot y to ax4
%
% Alternatively,
%   hands = plot_all_cycles('all') : plot uz, ze, x, and y to the figures
%                                    generated by repeated
%                                    calls to figure();
%  bounds = [upper, lower] for z_error
% TODO: Implement what follows
%   hands = plot_all_cycles(spec) : Plot the data specified by
%                                   spec to created figures.
%                                   Spec should be a string
%                                   such as 'xyz' or 'ux' which
%                                   would result in plots of x,
%                                   y, and ze, or u_z and x, respectivily.
%
  
if ~exist('to_pix', 'var')
  to_pix=false;
end
% First, parse the input. 
  if isa(ax1, 'char') && strcmp(ax1, 'all')
    plotx=true;
    ploty=true;
    plotuz=true;
    plotze=true;
    
    figure()
    ax1 = gca();
    figure();
    ax2 = gca();
    figure();
    ax3 = gca();
    figure();
    ax4 = gca();

  else % assume first arg is ax1
    plotuz = true;
  end
  % The case, unless we find out otherwise.
  plotze = false;
  plotx = false;
  ploty = false;
  
  if nargin  > 2
    plotze = true;
  end
  if nargin > 3
    plotx = true;
  end
  if nargin > 4
    ploty=true;
  end
    
  axs = [];
  % Setup the figures
  if plotuz
    hold(ax1, 'on')
    title(ax1, 'uz')
    grid(ax1, 'on')
    axs = [axs, ax1];
  end
  if plotze
    hold(ax2, 'on')
    title(ax2, 'z-err')
    grid(ax2, 'on')
    axs = [axs, ax2];
  end
  if plotx
    hold(ax3, 'on')
    title(ax3, 'x')
    grid(ax3, 'on')
    axs = [axs, ax3];
  end
  if ploty
    hold(ax4, 'on')
    title(ax4, 'y')
    grid(ax4, 'on')
    axs = [axs, ax4];
  end
  
  % figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
  % axs = {ax1, ax2, ax3, ax4};        

  indc = {   'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
          'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};


  state_seq = {'move', 'tdown', 'tsettle', 'scan', 'tup', 'connect'};

  %   Num_cycles = min([length(self.idx_state_s.move), length(self.idx_state_s.tdown),...
  %                     length(self.idx_state_s.tsettle), length(self.idx_state_s.scan),...
  %                     length(self.idx_state_s.tup)]);
  h_x = gobjects(length(state_seq),1);
  h_y = gobjects(length(state_seq),1);
  h_ze = gobjects(length(state_seq),1);
  h_uz = gobjects(length(state_seq),1);
  
  
  
  for k=1:length(state_seq)
    
    %   try
    %     idx_state = self.idx_state_s.(state_seq{k}){idx_cs_seq};
    %   catch
    %     keyboard
    %   end
    idx_cell = self.idx_state_s.(state_seq{k});
    if isempty(idx_cell)
      h_x(k) = [];
      h_y(k) = [];
      h_ze(k) = [];
      h_uz(k) = [];
      continue;
    end
    if plotuz
      h_uz(k) = plot_cell_of_timeseries(ax1, idx_cell, self.uz, 'color', indc{1, k});
      h_uz(k).DisplayName = indc{2,k};
    end
    if plotze
      h_ze(k) = plot_cell_of_timeseries(ax2, idx_cell, self.ze, 'color', indc{1, k});
      h_ze(k).DisplayName = indc{2,k};
    end
    if plotx
      h_x(k) = plot_cell_of_timeseries(ax3, idx_cell, self.x, 'color', indc{1, k});
      h_x(k).DisplayName = indc{2,k};
    end
    if ploty
      h_y(k) = plot_cell_of_timeseries(ax4, idx_cell, self.y, 'color', indc{1, k});
      %     y_k = self.y(idx_state);
      %     h_y(k) = plot(ax4, t_k, y_k, 'color', indc{1, k});
      h_y(k).DisplayName = indc{2,k};
    end
  end
  hands = [h_x, h_y, h_ze, h_uz];
legs = [];
  if plotuz
    legs(end+1) = legend(h_uz);
  end

  if plotze
    legs(end+1) = legend(h_ze);
  end
  if plotx
    legs(end+1) = legend(h_x);
  end
  if ploty
    legs(end+1) = legend(h_y);
  end

  if exist('bounds', 'var') && ~isempty(bounds)
    plot(ax2, [self.t(1), self.t(end)], [bounds(2), bounds(2)], '--k')
    plot(ax2, [self.t(1), self.t(end)], [bounds(1), bounds(1)], '--k')
  end
  linkaxes(axs, 'x');
end

function hand = plot_cell_of_timeseries(ax, idx_cell, data, varargin)
  hand = gobjects(1);
  for k=1:length(idx_cell)
    t = idx_cell{k}*AFM.Ts;
    y = data(idx_cell{k});
    hand = plot(ax, t, y, varargin{:});
  end

end