clear;
clc;
data.naca = 2408;
data.alpha = 4.0;
data.geometry = 2;
data.Mpanels = 100;
data.x_h = 1.0;     % flap hinge position
data.eta_f = 0; % flap deflection angle


[Cl, Cmle,Cp, x, z, pxnorm, pznorm, xvort, zvort] = dvm(data);
fprintf( 'Cl= %.4f \nCmle=%.4f\n', Cl, Cmle);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot camber line + normal vectors %%%%%

quiver(x(1:end-1),z(1:end-1),pxnorm, pznorm)
hold on; 
plot( x, z, 'b')
scatter(xvort, zvort, 'r')
hold off;



 