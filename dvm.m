function [Cl, Cmle,Cp, x, z, pxnorm, pznorm, xvort, zvort] = dvm(data)
    % alias data components
    naca = data.naca;
    alpha = deg2rad(data.alpha);
    geometry = data.geometry;
    M = data.Mpanels;
    x_h = data.x_h;     % flap hinge position
    eta_f = deg2rad(data.eta_f); % flap deflection angle
    
    N = M + 1;
    naca = num2str(naca);
    f = str2num(naca(1)) * 0.01;
    p = str2num(naca(2)) * 0.10;
    t = str2num(naca(3)) * 0.10 + naca(4)*0.01;
    
    if geometry == 1 % discretizacion lineal
       x = linspace(0,1, N); 
    elseif geometry == 2 % discretizacion full cosine
       x = 0.5 * (1- cos((0:N-1)*pi/(N-1)));
    end
    z = zeros(1,N);
    
    % fins 0 < x <p
    z(x<=p) = (f/p^2) * (2*p*x(x<=p) - x(x<=p).^2);
    
    % p < x > 1
    cond = x>p & x<= 1;
    z(cond) = (f/(1-p)^2) * (1 - 2*p + 2*p*x(cond) - x(cond).^2);
    
    % recalcul per flap x_h<x<1
    z(x>=x_h) = z(x>=x_h) -tan(eta_f)*(x(x>=x_h)-x_h);

    
    pchord = zeros(1,M);   % allocate storage space
    pxnorm = zeros(1,M);
    pznorm = zeros(1,M);
    pxtan = zeros(1,M);
    pztan = zeros(1,M);
    xvort = zeros(1,M);
    zvort = zeros(1,M);
    xctrl = zeros(1,M);
    zctrl = zeros(1,M);

    
    u = zeros(M,M);
    w = zeros(M,M);
    for i = 1:M  % calculate panels geometry

        x1 = x(i);  % panel nodes 
        x2 = x(i+1);
        z1 = z(i);
        z2 = z(i+1);

        % calculate panel's chord, normal and tangent vectors

        pchord(i) = sqrt((x2-x1)^2+(z2-z1)^2);
        pxnorm(i) = - (z2-z1)/pchord(i); % components of the normal vector
        pznorm(i) = (x2-x1)/pchord(i);
        pxtan(i) = (x2-x1)/pchord(i) ; % components of the tangent vector
        pztan(i) =  (z2-z1)/pchord(i); 

        % calculate position of the vortex and control point

        xvort(i) = x1 + 0.25*pchord(i)*pxtan(i);
        zvort(i) = z1 + 0.25*pchord(i)*pztan(i);
        xctrl(i) = x1 + 0.75*pchord(i)*pxtan(i);
        zctrl(i) = z1 + 0.75*pchord(i)*pztan(i);

    end
    
    for i = 1:M
        xi = xctrl(i);
        zi = zctrl(i);
        nxi = pxnorm(i);
        nzi = pznorm(i);
        for j = 1:M
           xv = xvort(j);
           zv = zvort(j);
           r_2 = (xi-xv)^2 + (zi-zv)^2;
           u(i,j) = (1/(2*pi*r_2)) * (zi-zv);
           w(i,j) = (-1/(2*pi*r_2)) * (xi-xv);
           A(i,j) = u(i,j)*nxi + w(i,j)*nzi;   
        end
        RHS(i) = -(cos(alpha)*nxi + sin(alpha)*nzi);  
    end
    G = inv(A)*RHS';

    Cl = double( 2 * sum(G));
    Cmle = -2 * sum(G'.*xvort*cos(alpha));
    Cp = 2 * G'./pchord;
    
    
    
    
    
    
    
    
end