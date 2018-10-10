Fscan = 1;

W = 5e-6;

V = Fscan * 2 * W;

d1 = 1e-6;

Tskip = d1/V


% Time to do zup-xymove,zeng, from slide 42 of comps
Tmove = 0.014; 
Tupdown = 0.066;
Tmu = Tmove + Tupdown;

% Even though Tmove depends on d1, it doesn't change that much. So assume it is
% constant. Then we can solve for d1:

d1_max = Tmu*V
d1_max_mic = d1_max*1e6
%%



G = tf(100, [1, 100]);
D = tf([1], [1 0]);

H = -feedback(G*D, 1)
bode(G*D, H)
grid on