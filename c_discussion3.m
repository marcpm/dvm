clear;
clc;
data.naca = 2408;
data.alpha = 0.0; % para calculo de alpha_l0 
data.geometry = 2;
data.Mpanels = 100;
data.x_h = 1.0;     % flap hinge position
data.eta_f = 0; % flap deflection angle

cl_alpha = 0.10936;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3.Discussion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 50;  % numero de puntos
alpha_l0_array = zeros(4,4);
Cm0_array = zeros(4,4);
i = 1;
fs = [0, 2, 4, 6];
ps = [1,2,4,6];
for f = fs
    j = 1;
    for p = ps
        data.alpha = 0;
        data.naca = strcat(num2str(f), num2str(p), '00') ;
        [Cl, Cmle,Cp, x, z, pxnorm, pznorm, xvort, zvort] = dvm(data);
%         data.alpha = 3.2;
%         [Cl2, Cmle2,Cp, x, z, pxnorm, pznorm, xvort, zvort] = dvm(data);

%         a_l0 = -Cl/0.10936;

        alpha_l0_array (i, j) = -Cl/cl_alpha;
        Cm0_array(i, j) = Cmle + Cl/4;

        j = j + 1;
    end
    i = i + 1;
end

plot(ps*0.1 , alpha_l0_array', '*-', 'LineWidth', 2)
ylabel('$\alpha_{l0} \; [deg]$', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('Max thickness Position (p)')
legend('f=0', 'f=0.01', 'f=0.02', 'f=0.04', 'f=0.06')
figure
plot(ps*0.1, Cm0_array', '*-', 'LineWidth', 2)
ylabel('$Cm_{0}$', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('Max thickness Position (p)')
legend('f=0', 'f=0.01', 'f=0.02', 'f=0.04', 'f=0.06')