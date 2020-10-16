%% Create ODEs for nonlinear systems
function xdot = pendulum(t,x,parameters)

    % extract pendulum parameters & states
    M = parameters.M;
    g = parameters.g;
    l = parameters.l;
    x1 = x(1); % theta
    x2 = x(2); % thetadot
    
    % SET control input
    u = 0;
    
    % construct nonlinear ODE model
    xdot = [x2; -(g/l)*sin(x1)-1/(M*l)*cos(x1)*u]; % note that sin(x) makes it nonlinear
end