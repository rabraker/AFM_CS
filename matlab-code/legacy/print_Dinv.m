
Ts = 40e-6;

D_inv_old = tf([1.02370858192, -2.04344248772, 1.02259731293], [1, -1.99611163139, 0.998974978924], Ts);

fname = 'Z:\afm-cs\sysID\z-axis_sines_info_quick_firstResFourierCoef_11-26-2018-03.mat'

modelFit = load(fname)

G1 = modelFit.g_1st_res;
G1.IOdelay = 0;
D_inv_new = 1/G1;
D_inv_new = D_inv_new/dcgain(D_inv_new);

figure
bode(D_inv_new, D_inv_old)


format long
[num, den] = tfdata(D_inv_new, 'v')