% clear;
clc;
data.naca = 2408;
data.alpha = 4.0;
data.geometry = 2;
data.Mpanels = 100;
data.x_h = 1.0;     % flap hinge position
data.eta_f = 0; % flap deflection angle


% [Cl, Cmle,Cp, x, z, pxnorm, pznorm, xvort, zvort] = dvm(data);
%     fprintf( 'Cl= %.4f \nCmle=%.4f\n', Cl, Cmle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% TAT Compare Values %%%%%%%
Cl_TAT = 0.66644;
Cmle_TAT = -0.21973;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 2 VALIDATION  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FIND CL_a, alpha_0, Cm0, x_ac %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linealidad, solo necesitariamos 2 puntos, pero se adoptan 5 puntos y una
% regresion por ellos.
n = 5;  % numero de puntos
i =1;
Cl_array = zeros(1,n);
Cmle_array = zeros(1,n);

alphas  = linspace(-3,5,n); % rango Mpaneles
for alpha = alphas
    data.alpha = alpha;
    [Cl, Cmle] = dvm(data);
       
    Cl_array(i) = Cl;
    Cmle_array(i) = Cmle;
    i = i + 1;
end

Cl_coefs = polyfit(alphas,Cl_array,1);
Cmle_coefs = polyfit(Cl_array,Cmle_array,1);

fprintf('2pi: %.4f \n', 2*pi);
fprintf('dCl/dalpha: %.5f \nalpha_l0: %.5f\n\n', Cl_coefs(1), -Cl_coefs(2)/Cl_coefs(1));
fprintf('Cm_0: %.4f\nx_ac: %.4f\n\n', Cmle_coefs(2), abs(Cmle_coefs(1)));

    
%  cl_alpha_error = abs(Clalpha_abbot-Clalpha)/Clalpha_abbot;
%  alpha_l0_error = abs(alpha_l0_abbot-alpha_l0Cl)/abs(alpha_l0_abbot);
%  Cm0_error = abs(Cm0_abbot-Cm0)/abs(Cm0_abbot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%COMMENTS%%%%%%%%%%%%%%%%%%%%
% Tomando 5 alphas, obtenemos una pendiente 
% Cl_alpha de 6.2658 rad-1 = 0.10936 deg-1
% alpha_0 = 2.05918 deg
% Cm_0 = -0.0536

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FLAP EFFICENCY%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 15;

fchord_ratios = [0.0,0.005,0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.40];
% fchord_ratios = linspace(0,0.4, 9);
f_efficencies = zeros(size(fchord_ratios));
Cls = zeros(1,n);
alpha_l0s = zeros(1,n);
eta_fs = linspace(0,20, n);
data.alpha = 0.0;
i = 1;
for fchord_ratio = fchord_ratios
    data.x_h = 1 - fchord_ratio;
    j = 1;
    for eta_f = eta_fs
        data.eta_f = eta_f;
        [Cl, Cmle] = dvm(data);
        Cls(j) = Cl;
        j = j + 1;
    end
    alpha_l0s =  Cls / 0.10936;
    temp = polyfit(eta_fs, alpha_l0s,1);
    f_efficencies(i) = temp(1);
    i = i + 1;
end

plot(fchord_ratios, f_efficencies,'Color',[0, 0.4470, 0.7410], 'LineWidth', 2)
hold on
requested_fchord_ratios = [0.15, 0.2, 0.25, 0.3];


% scatter(requested_fchord_ratios, f_efficencies(any(fchord_ratios==requested_fchord_ratios)), 'r')

%%% %%% %%% %%% %%% %%% %%% %%% 
%%% efficency DVM corrected %%% 
corrected_effs080 = f_efficencies*0.8;
corrected_effs083 = f_efficencies*0.83;
corrected_effs090 = f_efficencies*0.9;
plot(fchord_ratios, corrected_effs080, ':' ,'LineWidth', 2)
plot(fchord_ratios, corrected_effs083, ':' ,'LineWidth', 2)
plot(fchord_ratios, corrected_effs090, ':' ,'LineWidth', 2)


%%% %%% %%% %%% %%% %%% 
%%% efficency abbot %%% 
abbot_effs = [0.0, 0.125, 0.250, 0.35, 0.6, 0.65, 0.68];
abbot_chordratios = [0.0, 0.05, 0.105, 0.15, 0.3, 0.35, 0.4];
p = polyfit(abbot_chordratios, abbot_effs, 2);
abbot_effs_extrapol = polyval(p,fchord_ratios);
plot(fchord_ratios, abbot_effs_extrapol, '-', 'LineWidth', 2 )


%%%% dots enunciado 0.15.. 0.30 %%%%%%%
for fchord_ratio = requested_fchord_ratios
    eff = f_efficencies(fchord_ratios == fchord_ratio);
    color = [0.6350 0.0780 0.1840];
    scatter(fchord_ratio, eff,  [], color, 'filled')
    fprintf('flapchord_ratio: %.2f   flap_efficencies: %.4f \n',fchord_ratio,eff); 
    scatter(fchord_ratio, abbot_effs_extrapol(fchord_ratios == fchord_ratio),  [], [0, 0.5, 0], 'filled')
    scatter(fchord_ratio, corrected_effs080(fchord_ratios == fchord_ratio),  [], color, 'filled')
    scatter(fchord_ratio, corrected_effs083(fchord_ratios == fchord_ratio),  [], color, 'filled')
    scatter(fchord_ratio, corrected_effs090(fchord_ratios == fchord_ratio),  [], color, 'filled')
    
end

[hleg, hico] = legend('Flap Efficency DVM', 'Flap Efficency DVM 0.80-Corrected','Flap Efficency DVM 0.83-Corrected', 'Flap Efficency DVM 0.90-Corrected','Flap Efficency Experimental-Abbot', 'Location', 'southeast');
% Delete objects associated with last 3 black lines
istxt = strcmp(get(hico, 'type'), 'text');
hicot = hico(istxt);
hicol = hico(~istxt);
delete(hicot(ismember(get(hicot, 'String'), {'data10','data11','data12'})));
delete(hicol(ismember(get(hicol, 'Tag'),    {'data10','data11','data12'})));

ylabel('Flap Efficency  $\frac{\partial{\alpha_{l0}}}{\partial{\eta_f}}$', 'Interpreter', 'latex', 'FontSize', 20);

xlabel('Flap-Chord Ratio $E$', 'Interpreter', 'latex', 'FontSize', 16);

hold off
