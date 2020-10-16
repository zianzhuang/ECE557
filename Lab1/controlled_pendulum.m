function Zdot = pendulum(t,Z,parameters,Matrices)
    
    % extract parameters & states
    M = parameters.M;
    g = parameters.g;
    l = parameters.l;
    F = Matrices.F;
    G = Matrices.G;
    H = Matrices.H;
    L = Matrices.L;
    x1 = Z(1); % theta
    x2 = Z(2); % thetadot
    z = Z(3); % theta - controller
    
    % SET control input
    y = x1; % output of pendulum system/controller input
    u = H*z - L*y; % control input of pendulum system/controller output
    
    % construct ODE model
    xdot = [x2;(-g/l)*sin(x1)-1/(M*l)*cos(x1)*u]; % NL pendulum state
    zdot = F*z - G*y; % controller state
    Zdot = [xdot;zdot];
end