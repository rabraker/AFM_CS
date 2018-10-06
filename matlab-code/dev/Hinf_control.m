clear, clc, close all

matroot = getMatPath();
notepath = '/home/arnold/gradschool/proposals/nsf_afm_sept2017/latex/figs';
model_path = fullfile(matroot, 'publications', 'afmMPC', 'mimo_modelData.mat');
load(model_path);

PLANT_preT_x.InputDelay = 0;
G_x = PLANT_preT_x;
Ts = G_x.Ts;

We = tf([4, -0.9841], 100*[1 -1], Ts);
Wu = tf(0.1, [1], Ts);
W = tf([2 -1.595], 2.024*[1 0], Ts);


bode(G_x, We, Wu, W )
grid on
legend('G_x', 'We', 'Wu', 'W')