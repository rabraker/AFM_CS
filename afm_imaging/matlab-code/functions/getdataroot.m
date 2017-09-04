function [dataroot] = getdataroot(  )
if ispc
    dataroot = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
else
    dataroot = '/home/arnold/gradschool/afm-cs/afm_imaging/acc2018-data';
end


end

