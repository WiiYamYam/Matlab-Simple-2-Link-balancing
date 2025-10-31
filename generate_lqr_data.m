function generate_lqr_data()
% Parameters
g = 9.81;

% Masses (kg)
mw = 0.25; % wheel
mp = 0.45; % pendulum
mc = 5; % chassis

%links
mp_upper = 0.2;
mp_lower = 0.25;

L0_calc = (0.05:0.01:0.5);

Ks = zeros(2,6,length(L0_calc));
for i = 1:length(L0_calc)
    %lengths
    lwp = (L0_calc(i)*0.5); %length between wheel and pend COM
    lpc = lwp; %length between pend COM and chassis
    lc = 0;   %length between chassis COM and pendulum rotation joint
    r = 0.08;   %wheel radius
    l1 = 0.25;
    l2 = 0.3;
    
    
    %moments of inertia
    Iw_mujo = inertia_cylinder(mw,r,0.03); %wheel moment of inertia
    Ip_mujo = inertia_box(mp,0.04,0.04,L0_calc(i)); %pendulum moment of inertia
    Ic_mujo = inertia_box(mc,0.35,0.03,0.2); %chassis moment of inertia
    
    %pendulum moment of inertia - combine two links
    Ip_thigh = inertia_box(mp_upper,0.04, 0.04, l1);
    Ip_calf = inertia_box(mp_lower, 0.04, 0.04, l2);
    
    %Ip_mujo = Ip_thigh + Ip_calf;
    
    Iw = Iw_mujo(2,2);
    Ip = Ip_mujo(2,2);
    Ic = Ic_mujo(2,2);
    
    %system lol
    syms x xdot xddot theta thetadot thetaddot phi phidot phiddot tauw taup N Nc Pw Pc real
    
    %对驱动轮有
    eq_wheel = xddot == (tauw-N*r)/(mw*r+(Iw/r));
    
    %对摆杆有
    eq_pend_vert = Pw - Pc - (mc*g) == -mp*lwp*(thetaddot*sin(theta)+((thetadot^2)*cos(theta)));
    eq_pend_hori = N - Nc == mp*(xddot + lwp*(thetaddot*cos(theta)-((thetadot^2)*sin(theta))));
    eq_pend_turn = Ip*thetaddot == (Pw*lwp+Pc*lpc)*sin(theta)-(N*lwp+Nc*lpc)*cos(theta)-tauw+taup;
    
    %对机体有
    A = lwp + lpc;
    %C = lpc + lwp;
    
    eq_chassis_vert = Pc - mc*g == -mc*(A*thetaddot*sin(theta)+A*thetadot*thetadot*sin(theta)+lc*phiddot*sin(phi)+lc*phidot*phidot*sin(phi));
    eq_chassis_hori = Nc == mc*(xddot-A*thetaddot*sin(theta)-A*thetadot*thetadot*sin(theta)+lc*phiddot*sin(phi)+lc*phidot*phidot*sin(phi));
    eq_chassis_turn = Ic*phiddot == taup + Nc*lc*cos(phi) + Pc*lc*sin(phi);
    
    % Equations of motion
    eqns = [
        eq_wheel;
        eq_pend_vert;
        eq_pend_hori;
        eq_pend_turn;
        eq_chassis_vert;
        eq_chassis_hori;
        eq_chassis_turn];
    
    vars = [xddot, thetaddot, phiddot, Pw, Pc, N, Nc];
    
    % solve for accelerations, eliminating forces
    sol = solve(eqns, vars);
    
    xddot_sol    = simplify(sol.xddot);
    thetaddot_sol = simplify(sol.thetaddot);
    phiddot_sol   = simplify(sol.phiddot);
    
    % Define state and input vectors
    states = [x; xdot; theta; thetadot; phi; phidot];
    inputs = [tauw; taup];
    
    % Express derivatives in state-space form
    f = [xdot;
        xddot_sol;
        thetadot;
        thetaddot_sol;
        phidot;
        phiddot_sol];
    
    A_sym = simplify(jacobian(f, states));
    B_sym = simplify(jacobian(f, inputs));
    
    %compute A and B and linearise
    A_linear = double(subs(A_sym, {theta,phi,xdot,thetadot,phidot,tauw,taup}, {0,0,0,0,0,0,0}));
    B_linear = double(subs(B_sym, {theta,phi,xdot,thetadot,phidot,tauw,taup}, {0,0,0,0,0,0,0}));

    Q = diag([1/(0.01^2), 1/(1^2), 1/(2^2), 1/(0.01^2), 1/(1^2), 1/(2^2)]);
    R = diag([0.25, 0.25]);

    % Compute K(L0)
    Ks(:,:,i) = lqr(A_linear, B_linear, Q, R);

end

    % Save lookup data
    save('K_lookup_data_good.mat', 'L0_calc', 'Ks');
    disp('Saved lookup table K_lookup_data_good.mat');
end


%% (fun)c ;-;
function I = inertia_box(m,x,y,z)
Ixx = (1/12)*m*(y*y+z*z);
Iyy = (1/12)*m*(x*x+z*z);
Izz = (1/12)*m*(x*x+y*y);

I = diag([Ixx,Iyy,Izz]);
end

function I = inertia_cylinder(m,r,h) %wheel is standing on the y axis!!!
% Cylinder whose axis is along y
Iyy = 0.5*m*r^2;                       % rotation axis
Ixx = (1/12)*m*(3*r^2 + h^2);          % perpendiculars
Izz = Ixx;
I = diag([Ixx,Iyy,Izz]);
end