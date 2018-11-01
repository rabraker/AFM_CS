

w = 100;

g = tf([1, 2*w*0.3, (w+0.5)^2], [1, 2*w*0.2, w^2]);
gz = c2d(g,.001)
sys_z = ss(gz)
figure, step(sys_z)


[y,t] = step(sys_z);
size(t)


u = y*0+1;


[a,b,c,d] = ssdata(balreal(sys_z))
%%
PHI = [a, b; c d]

PHI_col = [];
for k=1:size(PHI,1);
  PHI_col = [PHI_col; PHI(k,:)'];
end

root = 'C:\Users\arnold\Documents\afm-cs\labview\UnitTests\fpga_harnesses\SS_Compensators';

fname_model = fullfile(root, 'model.csv');
fname_traj = fullfile(root, 'traj.csv');
csvwrite(fname_model, PHI)
csvwrite(fname_traj, [u, y])