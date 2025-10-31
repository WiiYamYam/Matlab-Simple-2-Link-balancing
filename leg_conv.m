function tau = leg_conv()
syms F Tp phi1 phi2
    % lengths
    l1 = 0.25; %thigh
    l2 = 0.3; %calf 

    xc = l1*cos(phi1)+l2*cos(phi1 + phi2);
    yc = l1*sin(phi1)+l2*sin(phi1 + phi2);
    p = [xc;yc];

    p_hat = p/norm(p); %unit vector
    size(p_hat)

    F_vec = F * p_hat; %force vector
    
    %virtual line polar coords
    %l0 = sqrt(xc^2 +yc^2);
    %phi0 = atan2(yc,xc);
    %pos = [l0;phi0];
    
    % Jacobians (in cartesian space)
    J = [-l1*sin(phi1)-l2*sin(phi1+phi2), -l2*sin(phi1+phi2);
          l1*cos(phi1)+l2*cos(phi1+phi2), l2*cos(phi1+phi2)];
    
    F_lqr = pinv(J') * [Tp;0];
    %virtual jacobian (polar coordinate)
    %J_v = [diff(l0,phi1), diff(l0, phi2);
              %diff(phi0,phi1), diff(phi0, phi2)];

    F_tot = F_vec + F_lqr;
    %VMC coversion cartesian
    tau = J' * F_tot;

    % VMC conversion: [T1; T2] = J_virtual' * [F; Tp]
    %conv = @(F, Tp) J' * [F_vec; Tp];
    %tau = conv(F, Tp);

    matlabFunction(tau, 'Vars', {F, Tp, phi1, phi2}, 'File', 'leg_convolution');
end