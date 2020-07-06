clear;
clc;

%%%%%%%%%%%%%%%%%%
%%% DATA INPUT %%%
%%%%%%%%%%%%%%%%%%
data.naca = 2408;
data.alpha = 4.0;
data.geometry = 2;
data.Mpanels = 100;
data.x_h = 1.0;     % flap hinge position
data.eta_f = 0; % flap deflection angle

%%%%%%% TAT Analytical Values to Compare %%%%%%%
Cl_TAT = 0.66644;
Cmle_TAT = -0.21973;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1.Verification cl-conv-error, & cmle -conv-error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 50;  % numero de puntos
i =1;

Cl_array = zeros(1,n);
Cl_error_array = zeros(1,n);
Cmle_array = zeros(1,n);
Cmle_error_array = zeros(1,n);

M_array  = floor(linspace(2,200,n)); % rango Mpaneles
for M = M_array
    data.Mpanels = M;
    [Cl, Cmle,Cp, x, z, pxnorm, pznorm, xvort, zvort] = dvm(data);
    
    Cl_error = abs(Cl_TAT-Cl)/Cl_TAT;
    Cmle_error = abs(Cmle_TAT-Cmle)/abs(Cmle_TAT);
    
    Cl_array(i) = Cl;
    Cl_error_array(i) = Cl_error;
    Cmle_array(i) = Cmle;
    Cmle_error_array(i) = Cmle_error;
    
    i = i + 1;
    
end
figure
hold ( 'on')
yyaxis ( 'left')
plot( M_array, Cl_array, '.-', 'LineWidth', 2);
ylabel( '$C_l$', 'Interpreter','latex','FontSize', 20)

yyaxis( 'right')
lightOrange = [212, 122, 95] / 255; 
ax = gca;
ax.YColor = lightOrange;
plot( M_array, Cl_error_array*100, '.-', 'LineWidth', 2, 'Color', lightOrange);
ylabel( '$C_{l} \; error \; [\%]$', 'Interpreter','latex','FontSize', 20)

xlabel( 'Número Paneles (M)','FontSize', 14)
hold('off')

figure
hold ( 'on')
yyaxis ( 'left')
plot( M_array, Cmle_array, '.-', 'LineWidth', 2);
ylabel( '$Cm_{le}$', 'Interpreter','latex','FontSize', 20)

yyaxis ( 'right')
lightOrange = [212, 122, 95] / 255; 
ax = gca;
ax.YColor = lightOrange;
plot( M_array, Cmle_error_array*100, 'g.-', 'LineWidth', 2, 'Color', lightOrange);

ylabel('$Cm_{le} \; error \; [\%]$', 'Interpreter','latex','FontSize', 20)

xlabel( 'Número Paneles (M)','FontSize', 14)
hold('off')


data.Mpanels = 120;
[Cl, Cmle,Cp] = dvm(data);
fprintf( '\n  Cl= %.4f \n  Cmle=%.4f\n', Cl, Cmle);
fprintf( 'Cl_TAT= %.4f \nCmle_TAT=%.4f\n', Cl_TAT, Cmle_TAT);


 

