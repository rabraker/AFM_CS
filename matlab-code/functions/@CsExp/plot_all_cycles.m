function [hands] = plot_all_cycles(self, ax1, ax2, ax3, ax4)
%% Plot the color-ized cs cycles to the provided axes.
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
% 
% TODO: Implement what follows
%   hands = plot_all_cycles(spec) : Plot the data specified by
%                                   spec to created figures.
%                                   Spec should be a string
%                                   such as 'xyz' or 'ux' which
%                                   would result in plots of x,
%                                   y, and ze, or u_z and x, respectivily.
%
  

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

  indc = {   'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], ;
          'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up',};


  state_seq = {'move', 'tdown', 'tsettle', 'scan', 'tup'};

  Num_cycles = min([length(self.idx_state_s.move), length(self.idx_state_s.tdown),...
                    length(self.idx_state_s.tsettle), length(self.idx_state_s.scan),...
                    length(self.idx_state_s.tup)]);
h_x = gobjects(length(state_seq),1);
h_y = gobjects(length(state_seq),1);
h_ze = gobjects(length(state_seq),1);
h_uz = gobjects(length(state_seq),1);
  for idx_cs_seq = 1:Num_cycles

    for k=1:length(state_seq)
      try
        idx_state = self.idx_state_s.(state_seq{k}){idx_cs_seq};
      catch
        keyboard
      end
      
      t_k  = idx_state*self.Ts;
      if plotuz
        uz_k = self.uz(idx_state);
        h_uz(k) = plot(ax1, t_k, uz_k, 'color', indc{1, k});
        h_uz(k).DisplayName = indc{2,k};
      end
      if plotze
        ze_k = self.ze(idx_state);
        h_ze(k) = plot(ax2, t_k, ze_k, 'color', indc{1, k});
        h_ze(k).DisplayName = indc{2,k};
      end
      if plotx
        x_k = self.x(idx_state);
        h_x(k) = plot(ax3, t_k, x_k, 'color', indc{1, k});
        h_x(k).DisplayName = indc{2,k};
      end
      if ploty
        y_k = self.y(idx_state);
        h_y(k) = plot(ax4, t_k, y_k, 'color', indc{1, k});
        h_y(k).DisplayName = indc{2,k};
      end
    end

  end

  if plotuz
    legend(h_uz);
  end

  if plotze
    legend(h_ze);
  end
  if plotx
    legend(h_x);
  end
  if ploty
    legend(h_y);
  end

  plot(ax2, [self.t(1), self.t(end)], [.05, .05], '--k')
  plot(ax2, [self.t(1), self.t(end)], -[.05, .05], '--k')

  linkaxes(axs, 'x');
end