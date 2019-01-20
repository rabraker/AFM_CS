clear
clc
load notes/data/cs_sim_CS20NG.mat
whos

%%
cs_sim.Img_sub_sampled = cs_sim.Img_original.*cs_sim.pix_mask;

cs_sim.solve_bp()

%%
figure(1)
thresh0 = min(cs_sim.Img_original(:));
thresh1 = max(cs_sim.Img_original(:));
imshow(cs_sim.Img_original, [thresh0, thresh1])

u = PixelMatrixToVector(cs_sim.Img_original);



%%
% cs_sim.Img_bp = cs_sim.Img_bp - min(cs_sim.Img_bp(:));
Ts = 1/(2*512);
g = zpk([], 0.7, 1-0.5, Ts);
g = g/dcgain(g);


clc
t = (0:length(u)-1)'*Ts;

y = lsim(g, u, t);
y = lsim(g, flipud(y), t);
y = flipud(y);
Img_filt = PixelVectorToMatrix(y, [512, 512]);

figure(3);clf
ImshowDataView.setup(gcf)
ax1 = subplot(3, 1, [1,2]);
ax2 = subplot(3, 1, 3);
cb =  @(event_obj)cs_sim.dataview_callback(event_obj, ax1, ax2);
h = ImshowDataView.imshow(cs_sim.Img_bp, [thresh0, thresh1], ax1, ax2, cb)




%%
ImshowDataView.imshow(img_mat, [0,1], ax1, ax2);
%%

cs_sim_filt = CsSim(Img_filt, cs_sim.pix_mask);


cs_sim_filt.solve_bp()

figure(4);

imshow(cs_sim_filt.Img_bp, [thresh0, thresh1])










