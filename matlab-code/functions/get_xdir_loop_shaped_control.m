function xdirControl = get_xdir_loop_shaped_control(verbose)
    if nargin < 1
        verbose = true;
    end
    
    verbose = true;
    addpath(fullfile(getCsRoot(), 'matlab-code', 'functions'));
    
    % ------- Load Plants -----
    [plants, frf_data] = CanonPlants.plants_ns14(9, '5micron');
    
    newest_xfrf_path = fullfile(PATHS.sysid, 'x-axis_sines_info_first_resFourierCoef_5-31-2019-03.json');
    sysx_new = SweptSinesOnline(newest_xfrf_path);
    
    
    [z, p, k] = zpkdata(plants.Gvib, 'v');
    p1 = p(end-1:end);
    z1 = z(end-3:end-2);
    p(end-1:end) = [];
    z(end-3:end-2) = [];
    
    g_tmp = zpk(z, p, k, plants.Gvib.Ts)*plants.gdrift;
    g_tmp = g_tmp*dcgain(plants.Gvib*plants.gdrift)/dcgain(g_tmp);
    
    Gx_bend_ = sysx_new.FRF_from_FC(1, [2]);
    Gx_bend = Gx_bend_/g_tmp;
    
    % Initial guess.
    gx_bend0 = zpk(z1, p1, 1, plants.Gvib.Ts);
    
    sos_fos_xbend = SosFos(gx_bend0, 'iodelay', 9);
    lg = LogCostZPK(squeeze(Gx_bend.Response), Gx_bend.Frequency*2*pi, sos_fos_xbend);
    lg.solve_lsq(1);
    gx_bend_lg = lg.sos_fos.realize();
    
    if verbose
        figure(100)
        opts = bodeoptions();
        opts.FreqUnits = 'Hz';
        bodeplot(Gx_bend, opts)
        grid
        hold on
        bodeplot(gx_bend_lg, Gx_bend.Frequency*2*pi);
    end
%     Dinv = 1/gx_bend0;
    Dinv = 1/gx_bend_lg;
    Dinv = Dinv/dcgain(Dinv);
    
    D_notch= zpk([0.979 + 1i*0.163, 0.979  - 1i*0.163], [0.8, 0.8], 1, AFM.Ts);
    D_notch = D_notch/dcgain(D_notch);
    

    
    Ki = 0.03;
    Dki = zpk(0, 1, Ki, AFM.Ts);
    % THe gain crossover is nominally at 615 rad/s. The real pole at 400 hz and the
    % integrator collide and exit unit circle.
    wz = 2.1e3;
    wp = 7*wz;
    kk = wp/wz;
    D_lead_ = kk*tf([1, wz], [1, wp]);
    D_lead = c2d(D_lead_, AFM.Ts, 'matched');

    Dx = balreal(Dinv)*balreal(D_notch); %*balreal(D_lead);
    Loop = balreal(minreal(absorbDelay(ss(plants.Gvib))*Dx * Dki));
    
    Hy_rprime = minreal(balreal(feedback(Loop, 1)));
    Dx_ff = make_dx_ff();
    Hyr = Hy_rprime * Dx_ff;
    if verbose
        figure(1)
        step(Hy_rprime, Hyr)
        figure(2); clf
        bode(Loop, Hy_rprime)
        grid on
        figure(3)
        rlocus(Loop)
        xlim([0.7, 1.1])
        ylim([-0.6, 0.6])
    end
    
    bw = bandwidth(Hy_rprime)/2/pi;
    [gm, pm, wcg, wcp] = margin(Loop);
    
    fprintf('Bandwidth %f Hz\n', bw);
    fprintf('Gain margin: %f\n', 20*log10(gm));
    fprintf('Phase margin: %f\n', pm);
    
    Sens = minreal(1/(1+Loop));
    
    g_static = zpk([], [], 1, AFM.Ts);
    xdirControl = struct('Sens', Sens,...
        'Hyr', Hyr,...
        'Hy_rprime', Hy_rprime,...
        'Loop', Loop,...
        'plants', plants,...
        'Ki', Ki,...
        'D_ki', Dki,...
        'D', Dx,...
        'D_ff', make_dx_ff(),...
        'D1_loop', g_static,...
        'D2_ss_ff', g_static);
    
    

end


function Dx_ff = make_dx_ff()
    wz1 = 214 * 2 * pi;
    wz2 = 505 * 2 * pi;
    zz = 0.01;
    
    
    g = zpk(tf([1, 2*zz*wz1, wz1^2], conv([1, 350*2*pi],  [1, 350*2*pi])));
    g2 = zpk(tf([1, 2*zz*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi])));
    w_lpf = 2 * pi * 500;
    g3 = zpk([], [-w_lpf], w_lpf);
    
    g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched') * c2d(g3, AFM.Ts, 'matched');
    Dx_ff = g * (1/dcgain(g));
end
