clc; clear; close all;

% Parameters
dt = 0.05; % Time step
N = 100; % Number of iterations
L = 60; % Length of board
W = 30; % Width of board
mu = 0.5; % Friction coefficient
g = 9.81; % Acceleration due to gravity

% Board setup
figure('position',[500 150 700 500]);
rectangle('Position', [-L/2 -W/2 L W], 'EdgeColor', 'k', 'LineWidth', 2);
rectangle('Position', [-L/2 -W/2 L W], 'FaceColor', [0.2 0.7 0.2]);

% Coin setup
r = 1; % Coin radius
xc = [10 -10 0 0 0]; % x-coordinate of coins
yc = [0 0 10 -10 0]; % y-coordinate of coins
vc = zeros(2,5); % Velocity of coins
hc = zeros(1,5); % Handle for coins
color = ['w', 'w', 'w', 'w', 'k']; % Color of coins

% Add striker
rs = 1.5*r; % Striker radius
xs = 0; % x-coordinate of striker
ys = 0; % y-coordinate of striker
vs = zeros(2,1); % Velocity of striker
hs = rectangle('Position', [xs-rs, ys-rs, 2*rs, 2*rs], 'Curvature', [1 1], 'FaceColor', 'r'); % Handle for striker

% Create handles for coins
for i = 1:5
    hc(i) = rectangle('Position', [xc(i)-r, yc(i)-r, 2*r, 2*r], 'Curvature', [1 1], 'FaceColor', color(i));
end

% Run simulation
for i = 1:N
    % Update positions
    xc = xc + vc(1,:)*dt;
    yc = yc + vc(2,:)*dt;
    xs = xs + vs(1)*dt;
    ys = ys + vs(2)*dt;
    
    % Handle collisions with board boundaries
    out_of_bounds = xc > L/2-r | xc < -L/2+r | yc > W/2-r | yc < -W/2+r;
    vc(:,out_of_bounds) = -mu*vc(:,out_of_bounds);
    xc(out_of_bounds) = xc(out_of_bounds) - 2*dt*vc(1,out_of_bounds);
    yc(out_of_bounds) = yc(out_of_bounds) - 2*dt*vc(2,out_of_bounds);
    vs(:,xs > L/2-rs | xs < -L/2+rs | ys > W/2-rs | ys < -W/2+rs) = -mu*vs(:,xs > L/2-rs | xs < -L/2+rs | ys > W/2-rs | ys < -W/2+rs);
     % Handle collisions with other coins and striker
    for j = 1:5
        for k = 1:5
            if j ~= k
                dist = norm([xc(j) yc(j)] - [xc(k) yc(k)]);
                if dist < 2*r
                    n = ([xc(j) yc(j)] - [xc(k) yc(k)])/dist;
                    vt1 = dot(vc(:,j), n);
                    vt2 = dot(vc(:,k), n);
                    if vt1 > vt2
                        vc(:,j) = vc(:,j) - (1+mu)*(vt1-vt2)*n;
                    end
                end
            end
        end
    end
end
