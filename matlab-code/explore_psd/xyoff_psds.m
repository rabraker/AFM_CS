% add sysID to path
clc
addpath('~/gradschool/sysID/matlab/functions')

% load the PSD data for z-disengaged
fname_dis = ['/media/labserver/afm-cs/sysID/multi-axis/' ...
             'z_PSD_disengaged.json'];

fname_en = ['/media/labserver/afm-cs/sysID/multi-axis/' ...
             'z_PSD_xyoff_zengaged.json'];

psd_dis = LvPSD(fname_dis);
psd_en = LvPSD(fname_en);


figure(1); clf
ax1 = subplot(3,2,1);

plot_psd_and_harmonic(psd_dis, 60, 20, ax1);
title(ax1, 'z-disengaged, 60-hz')

ax2 = subplot(3,2,3);
plot_psd_and_harmonic(psd_dis, 120, 20, ax2);
title(ax2, 'z-disengaged, 120-hz')
ax3 = subplot(3,2,5);
plot_psd_and_harmonic(psd_dis, 300, 20, ax3);
title(ax3, 'z-disengaged, 300-hz')

ax21 = subplot(3,2,2);

plot_psd_and_harmonic(psd_en, 60, 20, ax21);
title(ax21, 'z-enengaged, 60-hz')

ax22 = subplot(3,2,4);
plot_psd_and_harmonic(psd_en, 300, 20, ax22);
title(ax22, 'z-engaged, 120-hz')

ax23 = subplot(3,2,6);
plot_psd_and_harmonic(psd_en, 300, 20, ax23);
title(ax23, 'z-engaged, 300-hz')

figure(2); clf
clrs = {'k', '--g', ':r'};
ax1 = subplot(1,2,1);

hands1 = plot_psd_and_harmonic_s(psd_dis, [60, 120, 300], 20, ax1, clrs);
title(ax1, 'z-disengaged, 60,120,300-hz')
legend(hands1)
ax21 = subplot(1,2,2);

hands2 = plot_psd_and_harmonic_s(psd_en, [60,120,300], 20, ax21, clrs);
title(ax21, 'z-enengaged, 60-hz')

legend(hands2)




function hands = plot_psd_and_harmonic_s(psd, harm_freq_s, N, ax, harm_clrs)

  
  h_dis = semilogx(ax, psd.freqs, psd.psd);
  hold(ax, 'on')
  xlim(ax, [50, psd.freqs(end)]);

  ylm = ylim(ax);

  % compute and plot 60Hz harmonics.
  hands = gobjects(length(harm_freq_s),1);
  for j=1:length(harm_freq_s)
    harm_freq = harm_freq_s(j);
    clr = harm_clrs{j};
    harm = harm_freq*[1:N];

    for k=1:length(harm)
      h = plot(ax, [harm(k), harm(k)], ylm, clr);
    end
    h.DisplayName = sprintf('%.1f harmonic', harm_freq);
    hands(j) = h;
  end
  h_dis = semilogx(ax, psd.freqs, psd.psd);
  
end


function plot_psd_and_harmonic(psd, harm_freq, N, ax)

  h_dis = semilogx(ax, psd.freqs, psd.psd);
  hold(ax, 'on')

  xlim(ax, [50, psd.freqs(end)]);

  ylm = ylim(ax);

  % compute and plot 60Hz harmonics.

  harm = harm_freq*[1:N];

  for k=1:length(harm)
    plot(ax, [harm(k), harm(k)], ylm, ':r')
  end

end
