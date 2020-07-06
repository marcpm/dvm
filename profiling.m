%%% DATA INPUT %%%
%%%%%%%%%%%%%%%%%%
data.naca = 2408;
data.alpha = 4.0;
data.geometry = 2;
data.Mpanels = 100;
data.x_h = 1.0;     % flap hinge position
data.eta_f = 0; % flap deflection angle


for M = [60, 100,200]
    data.Mpanels = M;
    f = @() dvm(data);
    t1 = timeit(f,2);
    disp(t1)
end