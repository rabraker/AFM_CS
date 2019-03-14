

root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
files = {...
'first_res_fit-3-11-2019-1.json',...
'first_res_fit-3-4-2019.json',...
'first_res_fit-3-5-2019.json',...
'first_res_fit-3-7-2019.json',...
'first_res_fit-3-9-2019-2.json',...
};


h = gobjects(length(files), 1);
for k=1:length(files)
  name = strrep(files{k}, 'first_res_fit-', '')
  name = strrep(name, '-1.json', '')
  name = strrep(name, '-2.json', '')
  name = strrep(name, '.json', '')
  dat = loadjson(fullfile(root, files{k}))

  Fig = figure(1);
  gz_k = dat.SOS_frf_data.resp_real + dat.SOS_frf_data.resp_imag*1i;
  h(k) = frf_bode_mag(gz_k, dat.SOS_frf_data.omegas/2/pi, Fig, 'Hz')
  h(k).DisplayName = name;
end

legend(h)