%% Main code

% Code to accompany Chesebro, Mujica-Parodi, and Weistuch. "Ion
% Gradient-driven Bifurcations of a Multi-Scale Neuronal Model." Pre-print
% available at
% https://www.researchgate.net/publication/364088666_Ion_Gradient-driven_Bifurcations_of_a_Multi-Scale_Neuronal_Model.

% This is the main code file! Run this to generate the results for the
% paper. As you can see commented below, there are some areas with fixed
% values to reproduce the plots in the paper exactly. If you want to vary
% the initial conditions, simply comment out those lines and uncomment the
% lines directly below them to return to random values (which will generate
% similar - but not identical - results to the paper (more useful for your
% own studies).


%% Figure 1 is generated from BioRender.com
% If you want to generate your own figures, visit their site and have fun!

%% Figures 2a, 3a, and 4a
% These figures come from analyses done in MatCont. To recreate these
% analyses and plots, see the accompanying
% "running_matcont_larter_breakspear_bifurcations.md" file for
% information.

%% Generating Figures 2b, 3b and 4b
% The code below is to generate figures 2, 3, and 4. This is broken down 
% into sections in case you want to run one specific one/it takes too long
% to run all. First, run this section to initialize the variables that are 
% constant across all figures.

delta_V = 0.66; %variability of excitatory/inhibitory threshold
C = 0; %coupling constant (for single region, so no self-coupling)
u = 0; %coupling "matrix" (again, only single region)

%% Figure 2b: Na gradient effects
sim_len = 9e5; %ms
step_size = 0.1; %ms
fs = 1e3/step_size;

% Any initial parameters will work to generate a *similar* figure. If you
% want the exact figure from the paper, use these initial conditions.
y_0 = [-0.00827329309156305;-0.0450345569674258;0.902716109915281];

% Random initial conditions (comment line above and uncomment this to have
% your own initial conditions)
%y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];

V_K = -0.7; %default K reversal potential
V_Ca = 1.0; %default Ca reversal potential

range = 1e4:1.1e4;
x = 0:length(range)-1;
x = x/fs*1e3;

V_Na = 0.38;
[t, y1] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_Na = 0.53;
[t, y2] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_Na = 0.68;
[t, y3] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

% Select range of usable figure (100ms once it's past the initial
% oscillations)
range = 1e4:1.1e4;
x = 0:length(range)-1;
x = x/fs*1e3;

% Plot the figure
f = figure;
hold on
f.Position = [100, 100, 900, 500];
plot(x, y2(range, 1), 'LineWidth', 2)
plot(x, y1(range, 1), 'LineWidth', 2)
plot(x, y3(range, 1), 'LineWidth', 2)
xlabel('Time (ms)', 'FontSize', 20)
ylabel('V_i')
ax = gca;
ax.FontSize = 16;
legend({'Between NS and PD', 'Below NS', 'Above PD'})
title('Figure 2b')

%% Figure 2c: Na gradient attractor changes

f = figure;
hold on
f.Position = [100, 100, 700, 500];
lh = plot(y1(1e3:end, 3), y1(1e3:end, 2));
lh.Color = [0.8500, 0.3250, 0.0980, 0.05];
lh = plot(y2(1e3:end, 3), y2(1e3:end, 2));
lh.Color = [0, 0.4470, 0.7410, 0.05];
lh = plot(y3(1e3:end, 3), y3(1e3:end, 2));
lh.Color = [0.9290, 0.6940, 0.1250, 0.2];
xlabel('W', 'FontSize', 20)
xlim([-0.05, 0.9])
ylabel('Z')
ax = gca;
ax.FontSize = 16;
title('Figure 2c')

%% Figure 3b: Ca gradient effects

sim_len = 9e5; %ms
step_size = 0.1; %ms
fs = 1e3/step_size;

% Any initial parameters will work to generate a *similar* figure. If you
% want the exact figure from the paper, use these initial conditions.
y_0 = [0.0422331997796276;0.0270954220673924;0.0426598559350487];

% Random initial conditions (comment line above and uncomment this to have
% your own initial conditions)
%y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];

V_Na = 0.53; %default sodium reversal potential
V_K = -0.7; %default potassium reversal potential

V_Ca = 0.949;
[t, y1] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_Ca = 1.0;
[t, y2] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_Ca = 1.084;
[t, y3] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

% Select range of usable figure (100ms once it's past the initial
% oscillations)
range = (1.75e4:1.85e4);
x = 0:length(range)-1;
x = x/fs*1e3;

% Plot the figure
% Constants are added to the range just to align the peaks to make
% comparison across ion concentrations easier.
f = figure;
hold on
f.Position = [100, 100, 900, 500];
plot(x, y2(range, 1), 'LineWidth', 2)
plot(x, y1(range, 1), 'LineWidth', 2)
plot(x, y3(range-6e1, 1), 'LineWidth', 2)
xlabel('Time (ms)', 'FontSize', 20)
ylabel('V_i')
ax = gca;
ax.FontSize = 16;
legend({'Between NS and PD', 'Below NS', 'Above PD'})
title('Figure 3b')

% Figure 3c: Ca gradient attractor changes

f = figure;
hold on
f.Position = [100, 100, 700, 500];
lh = plot(y1(1e3:end, 3), y1(1e3:end, 2));
lh.Color = [0.8500, 0.3250, 0.0980, 0.05];
lh = plot(y2(1e3:end, 3), y2(1e3:end, 2));
lh.Color = [0, 0.4470, 0.7410, 0.05];
lh = plot(y3(1e3:end, 3), y3(1e3:end, 2));
lh.Color = [0.9290, 0.6940, 0.1250, 0.05];
xlabel('W', 'FontSize', 20)
ylabel('Z')
ax = gca;
ax.FontSize = 16;
title('Figure 3c')

%% Figure 4: K gradient effects
sim_len = 9e5; %ms
step_size = 0.1; %ms
fs = 1e3/step_size;
delta_V = 0.66;
C = 0;
u = 0;

% Random initial conditions
%y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];

% Any initial parameters will work to generate a *similar* figure. If you
% want the exact figure from the paper, use these initial conditions.
y_0 = [-0.0108434788242700;0.0169980902375737;0.381357224224074];

V_Na = 0.53;
V_Ca = 1.0;

range = 1e4:1.1e4;
x = 0:length(range)-1;
x = x/fs*1e3;

V_K = -0.59;
[t, y1] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_K = -0.7;
[t, y2] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_K = -0.9;
[t, y3] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

% Select range of usable figure (100ms once it's past the initial
% oscillations)
range = (0.45e4:0.55e4);
x = 0:length(range)-1;
x = x/fs*1e3;

% Plot the figure
% Constants are added to the range just to align the peaks to make
% comparison across ion concentrations easier.
%close(f)
f = figure;
hold on
f.Position = [100, 100, 900, 500];
plot(x, y2(range+8.1e2, 1), 'LineWidth', 2)
plot(x, y3(range-8e2, 1), 'LineWidth', 2)
plot(x, y1(range-2.2e2, 1), 'LineWidth', 2)
xlabel('Time (ms)', 'FontSize', 20)
ylabel('V_i', 'FontSize', 20)
ax = gca;
ax.FontSize = 16;
legend({'Below PD (physiological)', 'Far below PD (non-physiological)', 'Above PD'})
%title('Figure 4b')

% Figure 4c: K gradient attractor changes

f = figure;
f.Position = [100, 100, 700, 500];
hold on
lh = plot(y3(1e3:end, 3), y3(1e3:end, 2));
lh.Color = [0.8500, 0.3250, 0.0980, 0.008];
grid
lh = plot(y2(1e3:end, 3), y2(1e3:end, 2));
lh.Color = [0, 0.4470, 0.7410, 0.05];
lh = plot(y1(1e3:end, 3), y1(1e3:end, 2));
lh.Color = [0.9290, 0.6940, 0.1250, 0.05];
xlabel('W', 'FontSize', 20)
ylabel('Z')
ax = gca;
ax.FontSize = 16;
%title('Figure 4c')

%% Figure 5 Coherence plots
% As with Figures 2-4, this is split across three sections: one for each
% ion. Please note: running these will take a significant amount of time,
% as there are 50 replicates per point to generate a stable average.

% If you just want to reproduce the plots in Figure 5 directly, simply load
% in the fig5_50replicates.mat file and run the plotting section for each
% ion. That data file contains the data set used to generate the plots in
% the actual manuscript.

% If you want to try Pearson correlations instead of covariance as the
% coherence measure, change all_cov_xx to all_corrs_xx and that will give
% you the Pearson landscape.

%% Set up coupling constant range
% Run this section first! You need the variables here to run all the other
% sections, so run it before any of the others to avoid errors.

delta_V = 0.66; %excitatory/inhibitory threshold variance
c_range = 0.05:0.01:0.15; %coupling constant range
sim_len = 2e3; %ms
step_size = 1; %ms
fs = 1e3/step_size;

% Initialize connections to be weighted as coupled oscillations
u = [0, 1; 1, 0];


%% Computing Na replicates
v_na_range = 0.42:0.005:0.6; %sodium reversal potential range
num_replicates = 50; %number of replicates per grid point
all_corrs_na = zeros(length(c_range), length(v_na_range), num_replicates); %Pearson correlations for time series pairs
all_cov_na = zeros(length(c_range), length(v_na_range), num_replicates); %Covariance for time series pairs

%Included tic/toc functions to help get a sense of how long one replicate
%takes (so multiply by num_replicates to get total run time)
for replicate = 1:num_replicates
    tic
    for i = 1:length(c_range)
        C = c_range(i);
        for j = 1:length(v_na_range)
            V_Na = v_na_range(j);
            
            % Random initial conditions
            y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];
            
            V_K = -0.7;%potassium reversal potential
            V_Ca = 1.0;%calcium reversal potential
            
            [t, y] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
            
            r = corrcoef(y(:, 1), y(:, 2));
            all_corrs_na(i, j, replicate) = r(2);
            c = cov(y(:, 1), y(:, 2));
            all_cov_na(i, j, replicate) = c(2);
        end
    end
    toc
end
disp('Done Na replicates')

%% Plot Na coherence values (Figure 5a)
% Note: this is a 3D surface. To get the angle shown in the figure, simply
% rotate the view to look down on the surface (i.e., projected along the
% coherence plane)
% To visualize the whole surface, comment out the cross line to see 3D
% (cross is above 3D surface to be visible in 2D projection)
v_na_range = 0.42:0.005:0.6; %sodium reversal potential range
na_coherence = mean(all_cov_na, 3);
figure
surf(c_range, v_na_range(1:end), na_coherence(:, 1:end)')
xlabel('c', 'FontSize', 20)
ylabel('V_{Na}', 'FontSize', 20)
zlabel('Mean covariance', 'FontSize', 20)
shading('interp')
ax = gca;
ax.FontSize = 16;
colorbar
hold on
%plot3(0.1, 0.53, 0.8, 'k+', 'MarkerSize', 20, 'LineWidth', 5) %add cross to physiological value

%% Computing Ca replicates

v_ca_range = 0.954:0.005:1.034; %calcium reversal potential range
num_replicates = 50; %number of replicates per grid point
all_corrs_ca = zeros(length(c_range), length(v_ca_range), num_replicates);
all_cov_ca = zeros(length(c_range), length(v_ca_range), num_replicates);

%Included tic/toc functions to help get a sense of how long one replicate
%takes (so multiply by num_replicates to get total run time)
for replicate = 1:num_replicates
    tic
    for i = 1:length(c_range)
        C = c_range(i);
        for j = 1:length(v_ca_range)
            V_Ca = v_ca_range(j);
            % Random initial conditions
            y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];
            
            V_K = -0.7; %potassium reversal potential
            V_Na = 0.53; %sodium reversal potential
            
            [t, y] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
            r = corrcoef(y(:, 1), y(:, 2));
            all_corrs_ca(i, j, replicate) = r(2);
            c = cov(y(:, 1), y(:, 2));
            all_cov_ca(i, j, replicate) = c(2);
        end
    end
    toc
end

disp('done Ca');

%% Plot Ca coherence values (Figure 5b)
% Note: this is a 3D surface. To get the angle shown in the figure, simply
% rotate the view to look down on the surface (i.e., projected along the
% coherence plane)
% To visualize the whole surface, comment out the cross line to see 3D
% (cross is above 3D surface to be visible in 2D projection)
v_ca_range = 0.954:0.005:1.034; %calcium reversal potential range
ca_coherence = mean(all_cov_ca, 3);
figure
surf(c_range, v_ca_range(1:end), ca_coherence(:, 1:end)')
xlabel('c', 'FontSize', 20)
ylabel('V_{Ca}', 'FontSize', 20)
zlabel('Mean covariance', 'FontSize', 20)
shading('interp')
ax = gca;
ax.FontSize = 16;
colorbar
ylim([0.954, 1.034])
hold on
plot3(0.1, 1.0, 0.8, 'k+', 'MarkerSize', 20, 'LineWidth', 5) %add cross to physiological value

%% Computing K replicates

v_k_range = -0.62:-0.005:-0.8; %potassium reversal potential range
num_replicates = 50; %number of replicates per grid point
all_corrs_k = zeros(length(c_range), length(v_k_range), num_replicates);
all_cov_k = zeros(length(c_range), length(v_k_range), num_replicates);

for replicate = 1:num_replicates
    tic
    for i = 1:length(c_range)
        C = c_range(i);
        for j = 1:length(v_k_range)
            V_K = v_k_range(j);
            % Random initial conditions
            y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];
            
            V_Ca = 1.0; %calcium reversal potential
            V_Na = 0.53; %sodium reversal potential
            
            [t, y] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
            r = corrcoef(y(:, 1), y(:, 2));
            all_corrs_k(i, j, replicate) = r(2);
            c = cov(y(:, 1), y(:, 2));
            all_cov_k(i, j, replicate) = c(2);
            
        end
    end
    toc
end

disp('done K');

%% Plot K coherence values (Figure 5c)
% Note: this is a 3D surface. To get the angle shown in the figure, simply
% rotate the view to look down on the surface (i.e., projected along the
% coherence plane)
% To visualize the whole surface, comment out the cross line to see 3D
% (cross is above 3D surface to be visible in 2D projection)
k_coherence = mean(all_cov_k, 3);
figure
surf(c_range, v_k_range(1:end), k_coherence(:, 1:end)')
xlabel('c', 'FontSize', 20)
ylabel('V_{K}', 'FontSize', 20)
zlabel('Mean covariance', 'FontSize', 20)
shading('interp')
ax = gca;
ax.FontSize = 16;
colorbar
ylim([-0.8 -0.62])
hold on
plot3(0.1, -0.7,0.8, 'k+', 'MarkerSize', 20, 'LineWidth', 5) %add cross to physiological value


%% Supplemental Figure Plots
% This is the code to generate the supplemental attractors that support the
% main results. In each case, the first 1e3 timesteps are omitted to remove
% abberantly large trajectories caused by random initial conditions.
clear all; clc;

delta_V = 0.66; %variability of excitatory/inhibitory threshold
C = 0; %coupling constant (for single region, so no self-coupling)
u = 0; %coupling "matrix" (again, only single region)
sim_len = 9e5; %ms
step_size = 0.1; %ms


%% Figure ??? Na gradient attractors
% Any initial parameters will work to generate a *similar* figure. If you
% want the exact figure from the paper, use these initial conditions.
y_0 = [-0.00827329309156305;-0.0450345569674258;0.902716109915281];
% Random initial conditions (comment line above and uncomment this to have
% your own initial conditions)
%y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];
V_K = -0.7; %default K reversal potential
V_Ca = 1.0; %default Ca reversal potential

V_Na = 0.38;
[t, y1] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
V_Na = 0.53;
[t, y2] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
V_Na = 0.68;
[t, y3] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
V_Na = 0.46;
[t, y4] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);
V_Na = 0.58;
[t, y5] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

figure
lh = plot3(y1(3e5:end, 1), y1(3e5:end, 3), y1(3e5:end, 2));
lh.Color = [0, 0, 0, 0.05];
ax = gca;
ax.CameraPosition = [10.0116   -3.6245    1.9525];
ax.FontSize = 16;
fig = gcf;
fig.Position = [100, 100, 900, 700];
grid
xlim([-1 0.5])
ylim([0 0.9])
zlim([0.06 0.4])
title('Depolarized Na Attractor')

figure
lh = plot3(y2(3e5:end, 1), y2(3e5:end, 3), y2(3e5:end, 2));
lh.Color = [0, 0, 0, 0.05];
ax = gca;
ax.CameraPosition = [10.0116   -3.6245    1.9525];
ax.FontSize = 16;
fig = gcf;
fig.Position = [100, 100, 900, 700];
grid
xlim([-1 0.5])
ylim([0 0.9])
zlim([0.06 0.4])
title('Normal Na Attractor')

figure
lh = plot3(y3(3e5:end, 1), y3(3e5:end, 3), y3(3e5:end, 2));
lh.Color = [0, 0, 0, 0.05];
ax = gca;
ax.CameraPosition = [10.0116   -3.6245    1.9525];
ax.FontSize = 16;
fig = gcf;
fig.Position = [100, 100, 900, 700];
grid
xlim([-1 0.5])
ylim([0 0.9])
zlim([0.06 0.4])
title('Hyperpolarized Na Attractor')

%% Figure ??? Ca gradient attractors

% Any initial parameters will work to generate a *similar* figure. If you
% want the exact figure from the paper, use these initial conditions.
y_0 = [0.0422331997796276;0.0270954220673924;0.0426598559350487];

% Random initial conditions (comment line above and uncomment this to have
% your own initial conditions)
%y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];

V_Na = 0.53; %default sodium reversal potential
V_K = -0.7; %default potassium reversal potential

V_Ca = 0.949;
[t, y1] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_Ca = 1.0;
[t, y2] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_Ca = 1.134;
[t, y3] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

figure
lh = plot3(y1(1e3:end, 1), y1(1e3:end, 3), y1(1e3:end, 2));
lh.Color = [0, 0, 0, 0.05];
ax = gca;
ax.CameraPosition = [10.0116   -3.6245    1.9525];
ax.FontSize = 16;
fig = gcf;
fig.Position = [100, 100, 900, 700];
grid
xlim([-1 0.5])
ylim([0 0.9])
zlim([0.06 0.25])
title('Depolarized Ca Attractor')

figure
lh = plot3(y2(1e3:end, 1), y2(1e3:end, 3), y2(1e3:end, 2));
lh.Color = [0, 0, 0, 0.05];
ax = gca;
ax.CameraPosition = [10.0116   -3.6245    1.9525];
ax.FontSize = 16;
fig = gcf;
fig.Position = [100, 100, 900, 700];
grid
xlim([-1 0.5])
ylim([0 0.9])
zlim([0.06 0.25]) 
title('Normal Ca Attractor')

figure
lh = plot3(y3(1e3:end, 1), y3(1e3:end, 3), y3(1e3:end, 2));
lh.Color = [0, 0, 0, 0.05];
ax = gca;
ax.CameraPosition = [10.0116   -3.6245    1.9525];
ax.FontSize = 16;
fig = gcf;
fig.Position = [100, 100, 900, 700];
grid
xlim([-1 0.5])
ylim([0 0.9])
zlim([0.06 0.25])
title('Hyperpolarized Ca Attractor')

%% Figure ??? K gradient attractors
% Random initial conditions
y_0 = [0.1*rand(size(u, 1)*2, 1)-0.05; rand(size(u, 1), 1)];

% Any initial parameters will work to generate a *similar* figure. If you
% want the exact figure from the paper, use these initial conditions.
%y_0 = [-0.0108434788242700;0.0169980902375737;0.381357224224074];

V_Na = 0.53;
V_Ca = 1.0;

V_K = -0.59;
[t, y1] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_K = -0.7;
[t, y2] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);

V_K = -0.9;
[t, y3] = ode45(@(t,y) lb_for_ode45(t, y, delta_V, C, u, V_Na, V_K, V_Ca), [0:step_size:sim_len], y_0);




