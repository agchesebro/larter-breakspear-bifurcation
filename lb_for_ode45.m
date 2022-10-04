%% Larter-Breakspear model for use with ODE45 solver (can be used with other solvers as well)
% This is the implementation of the Larter-Breakspear model used in the
% manuscript. The constants are all from Table 1.

% To call this model on its own, use the command 
% [t, y] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
%
% This will return the integrated time series in y. The structure is:
% [V1, V2, ..., Vn, Z1, Z2, ..., Zn, W1, W2, ..., Wn] for n coupled masses
%
% You'll need to specify the following variables (default values given here):
% delta_V = 0.66; %excitatory/inhibitory voltage threshold variance
% C = 0.1; %coupling constant;
% u = [0, 1; 1, 0]; %matrix to specify connectivity between regions; two in this example
% V_Na = 0.53; %sodium reversal potential
% V_K = -0.7; %potassium reversal potential
% V_Ca = 1.0; %calcium reversal potential
% step_size = 1; %step size to save integrated time series at (in ms)
% sim_len = 2e3; %simulation length (in ms)
% y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)]; %random initial conditions for V, Z, and W

function dydt = lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca)

    %% Constants
    n_roi = size(u, 1);

    T_Ca = -0.01;
    delta_Ca = 0.15;
    g_Ca = 1;
    %V_Ca = 1;
    T_K = 0.0;
    delta_K = 0.3;
    g_K = 2.0;
    %V_K = -0.7;
    T_Na = 0.3;
    delta_Na = 0.15;
    g_Na = 6.7;
    %V_Na = 0.53;
    V_L = -0.5;
    g_L = 0.5;
    V_T = 0.0;
    Z_T = 0.0;
    delta_Z = delta_V;
    Q_Vmax = 1.0;
    Q_Zmax = 1.0;
    I = 0.3;
    a_ee = 0.36;
    a_ei = 2;
    a_ie = 2;
    a_ne = 1;
    a_ni = 0.4;
    b = 0.1;
    psi = 0.7;
    tau_K = 1;
    r_NMDA = 0.25;
    delta = 0;

    m_ion = @(V_i, T_ion, delta_ion) 0.5*(1+tanh((V_i-T_ion)/delta_ion));
    Q_V = @(V_i) 0.5*Q_Vmax*(1+tanh((V_i-V_T)/delta_V));
    Q_Z = @(Z_i) 0.5*Q_Zmax*(1+tanh((Z_i-Z_T)/delta_Z));
    
    % Test is sum(u) == 0
    % Only need in 1-2 ROI cases
    unreal = false;
    if sum(u(:, 1)) == 0 || sum(u(:, 2)) == 0
        unreal = true;
    end
    
    % split y apart into constitutive variables for easier indexing
    v = y(1:n_roi);
    z = y(n_roi+1:2*n_roi);
    w = y(2*n_roi+1:end);
    
    dydt = zeros(n_roi*3, 1);
    
    for i = 1:n_roi
        if ~unreal
           dv = -(g_Ca + ((1-C)*r_NMDA*a_ee*Q_V(v(i))) + (C*r_NMDA*a_ee*sum(u(:, i).*Q_V(v))/sum(u(:, i))))*m_ion(v(i), T_Ca, delta_Ca)*(v(i)-V_Ca)...
                -(g_K*w(i)*(v(i)-V_K))...
                -(g_L*(v(i)-V_L))...
                -(((g_Na*m_ion(v(i), T_Na, delta_Na)) + ((1-C)*a_ee*Q_V(v(i)))+(C*a_ee*sum(u(:, i).*Q_V(v))/sum(u(:, i))))*(v(i)-V_Na))...
                -((a_ie*z(i)*Q_Z(z(i))))...
                +(a_ne*I);
        else
           dv = -(g_Ca + ((1-C)*r_NMDA*a_ee*Q_V(v(i))) + (C*r_NMDA*a_ee*sum(u(:, i).*Q_V(v))))*m_ion(v(i), T_Ca, delta_Ca)*(v(i)-V_Ca)...
                -(g_K*w(i)*(v(i)-V_K))...
                -(g_L*(v(i)-V_L))...
                -(((g_Na*m_ion(v(i), T_Na, delta_Na)) + ((1-C)*a_ee*Q_V(v(i)))+(C*a_ee*sum(u(:, i).*Q_V(v))))*(v(i)-V_Na))...
                -((a_ie*z(i)*Q_Z(z(i))))...
                +(a_ne*I);
        end
         
         dz = b*((a_ni*I) + a_ei*v(i)*Q_V(v(i)));
         
         dw = psi*((m_ion(v(i), T_K, delta_K)-w(i))/tau_K);
         
         dydt(i) = dv;
         dydt(n_roi + i) = dz;
         dydt(2*n_roi + i) = dw;
    end