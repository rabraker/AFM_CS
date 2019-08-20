
fp_data = '/media/labserver/afm-cs/imaging/cs-imaging/5microns/6-5-2019/cs-traj-512pix-12perc-500nm-5mic-01Hz-150prescan-notconnect_out_6-5-2019-01.mat'

fp_caller = [mfilename('fullpath'), '.m']

FI_data = dir(fp_data)
FI_caller = dir(fp_caller)

FI_data.datenum < FI_caller.datenum