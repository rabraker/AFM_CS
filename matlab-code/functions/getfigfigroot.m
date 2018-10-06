
function fig_figroot = getfigfigroot()

if ispc()
    fig_figroot = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
else
    fig_figroot = '/home/arnold/gradschool/afm-cs/afm_imaging/matlab-codefigures';
end

end