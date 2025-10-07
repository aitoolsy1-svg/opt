clc; clear all; close all;

%% Configuration
adapter = 1;
res = 0.01;
yslices = 0.0; % [-0.5:0.1:0.5];

% this is the desired TCP pose
T_b_tcp_desired = [
    0   0 -1  0.5
    0  -1  0  0.0
    -1  0  0  -0.1
    0   0  0  1
];

T_b_tcp_desired(1:3, 1:3) = eul2rotm([0.001, -pi/2+1e-3, 0.001], 'zyx') % !!!change the zero here to small value like 0.001 to avoid singularity!!!
% T_b_tcp_desired(1, 4) = T_b_tcp_desired(1, 4) + 1e-3;

% this is the RCM pose
T_b_rcm_base = [
    1 0 0  0.5
    0 1 0  0
    0 0 1  0.0
    0 0 0  1
];

% Load the urdf model
robot = importrobot(['mops' num2str(adapter) '.urdf'], DataFormat='row');

% choose right IK function
if adapter == 1
    ik_mops = @(TCP, RCM) ik_mops1(TCP, RCM);
elseif adapter==2
    ik_mops = @(TCP, RCM) ik_mops2(TCP, RCM);
end

%% Testing
phi_sweep = [-90 : 15: 90];
theta_sweep = [-90 : 15: 90];
domains = [0:7];
F = 1;

wmin = 0.0;
wmax = 0;

data = {};

plotsDir = ['plots_xzA' num2str(adapter)];
if exist(plotsDir, 'dir') == 0
    mkdir(plotsDir);
end

% progress tracking
t0 = tic;
totalPhi = numel(phi_sweep);
totalTheta = numel(theta_sweep);
totalDomains = numel(domains); 
% Total tasks: (Part 1) one plot per theta per domain + (Part 2) one plot per phi per domain
totalTasks = (totalTheta + totalPhi) * totalDomains;
taskIdx = 0;

% Part 1: vary theta with fixed phi=0
for theta_deg = theta_sweep
    % avoid exact ±90 deg to prevent singularities
    if abs(theta_deg) == 90
        theta_deg = sign(theta_deg) * 89.99;
    end
    phi_deg = 0;
    phi = deg2rad(phi_deg);
    theta = deg2rad(theta_deg);
    R_ws = eul2rotm([phi, theta, 0], 'xyz'); % Rotation matrix for workspace
    T_b_rcm = T_b_rcm_base;
    T_b_rcm(1:3, 1:3) = R_ws * T_b_rcm(1:3, 1:3); % Rotate RCM orientation

for domain = domains

    iy = 1;
    yslice = 0.0;
        
        % Build grid in RCM frame (x_rcm, z_rcm) at y_rcm = 0
        [Xrcm, Zrcm] = meshgrid(-0.5:res:0.5, -0.5:res:0.0);
        [m, n] = size(Xrcm);
        W = zeros(m, n);
        q0 = NaN;
        
        
        % parallel computation
        parfor j = 1:n
            Wcol = zeros(m, 1);
            for i = 1:m
                % Map RCM point to world
                p_rcm = T_b_rcm(1:3, 4); % extracts the RCM position.
                p_world = T_b_rcm(1:3, 1:3) * [Xrcm(i, j); 0.0; Zrcm(i, j)] + p_rcm; % convert RCM to world frame.
                T = T_b_tcp_desired;
                %T(1:3, 1:3) = T_b_tcp_desired(1:3, 1:3); % keep TCP orientation fixed.
                T(1:3, 1:3) = R_ws * T(1:3, 1:3); % rotate TCP orientation by the same R_ws.
                T(1:3, 4) = p_world;
                solutions = ik_mops(T, T_b_rcm);
                if numel(solutions) < 1
                    Wcol(i) = 0;
                    continue
                end
                % use only selected postures
                D = classify_ur5_domains(solutions);
                I = all([D(:, 1) == (bitand(domain, 4)>0) D(:, 2) == (bitand(domain, 2)>0) D(:, 3) == (bitand(domain, 1)>0)], 2);
                solutions = solutions(I, :);
                sol = solutions;
                if numel(solutions) < 1
                    Wcol(i) = 0;
                    continue
                end
                q = q2urdf(sol);
                Wcol(i) = wconstr1(q, robot);
            end
            W(:, j) = Wcol;
        end
        
        % figure(F);
        % subplot(2, 4, domain+1);
        data{domain+1, iy, 1} = Xrcm;
        data{domain+1, iy, 2} = Zrcm;
        data{domain+1, iy, 3} = abs(W);

        % plot and save
        fig = figure;
        imagesc(Xrcm(1, :), Zrcm(:, 1), W);
        set(gca, 'YDir', 'normal');
        xlabel('x_{RCM}')
        ylabel('z_{RCM}')
        wmax_local = max(W, [], 'all');
        if isfinite(wmax_local) && wmax_local > 0
            clim([0 wmax_local]);
        end
        title(['w(q) on X''-Z'' (RCM) at Y=0, domain=' num2str(domain) ...
            ', angles: phi=' num2str(phi_deg) '°' ', theta=' num2str(theta_deg) '°'])
        base = ['xz_d' num2str(domain) '_y0_phi' num2str(phi_deg) '_theta' num2str(theta_deg)];
        domainDir = fullfile(plotsDir, ['domain_' num2str(domain)]);
        if exist(domainDir, 'dir') == 0
            mkdir(domainDir);
        end
        saveas(fig, fullfile(domainDir, [base '.fig']));
        saveas(fig, fullfile(domainDir, [base '.png']));
        % close;
        
        F = F+1;

        % update progress after each plot
        taskIdx = taskIdx + 1;
        elapsed = toc(t0);
        progress = taskIdx / totalTasks;
        if progress > 0
            estTotal = elapsed / progress;
            eta = estTotal - elapsed;
        else
            estTotal = NaN; eta = NaN;
        end
        disp(sprintf('Progress: %d/%d (%.1f%%) | phi=%d°, theta=%d°, domain=%d | elapsed=%.1fs, ETA=%.1fs', ...
            taskIdx, totalTasks, 100*progress, round(phi_deg), round(theta_deg), domain, elapsed, eta));

end

NF=F-1;

end

% Part 2: vary phi with fixed theta=0
for phi_deg = phi_sweep
    theta_deg = 0;
    if abs(phi_deg) == 90
        phi_deg = sign(phi_deg) * 89.99;
    end
    % compute workspace rotation
    phi = deg2rad(phi_deg);
    theta = deg2rad(theta_deg);
    R_ws = eul2rotm([phi, theta, 0], 'xyz');
    % update RCM orientation
    T_b_rcm = T_b_rcm_base;
    T_b_rcm(1:3, 1:3) = R_ws * T_b_rcm_base(1:3, 1:3);

for domain = domains

    iy = 1;
    yslice = 0.0;
        
        % Build grid in RCM frame
        [Xrcm, Zrcm] = meshgrid(-0.5:res:0.5, -0.5:res:0.0);
        [m, n] = size(Xrcm);
        W = zeros(m, n);
        q0 = NaN;
        
        
        
        parfor j = 1:n
            Wcol = zeros(m, 1);
            for i = 1:m
                p_rcm = T_b_rcm(1:3, 4);
                p_world = T_b_rcm(1:3, 1:3) * [Xrcm(i, j); 0.0; Zrcm(i, j)] + p_rcm;
                T = T_b_tcp_desired;
                %T(1:3, 1:3) = T_b_tcp_desired(1:3, 1:3);
                T(1:3, 1:3) = R_ws * T(1:3, 1:3); % rotate TCP orientation by the same R_ws.
                T(1:3, 4) = p_world;
                solutions = ik_mops(T, T_b_rcm);
                if numel(solutions) < 1
                    Wcol(i) = 0;
                    continue
                end
                D = classify_ur5_domains(solutions);
                I = all([D(:, 1) == (bitand(domain, 4)>0) D(:, 2) == (bitand(domain, 2)>0) D(:, 3) == (bitand(domain, 1)>0)], 2);
                solutions = solutions(I, :);
                sol = solutions;
                if numel(solutions) < 1
                    Wcol(i) = 0;
                    continue
                end
                q = q2urdf(sol);
                Wcol(i) = wconstr1(q, robot);
            end
            W(:, j) = Wcol;
        end
        
        data{domain+1, iy, 1} = Xrcm;
        data{domain+1, iy, 2} = Zrcm;
        data{domain+1, iy, 3} = abs(W);

        % plot and save immediately in domain subdir (RCM frame axes)
        fig = figure;
        imagesc(Xrcm(1, :), Zrcm(:, 1), W);
        set(gca, 'YDir', 'normal');
        xlabel('x_{RCM}')
        ylabel('z_{RCM}')
        wmax_local = max(W, [], 'all');
        if isfinite(wmax_local) && wmax_local > 0
            clim([0 wmax_local]);
        end
        title(['w(q) on X''-Z'' (RCM) at Y=0, domain=' num2str(domain) ...
            ', angles: phi=' num2str(phi_deg) '°' ', theta=' num2str(theta_deg) '°'])
        base = ['xz_d' num2str(domain) '_y0_phi' num2str(phi_deg) '_theta' num2str(theta_deg)];
        domainDir = fullfile(plotsDir, ['domain_' num2str(domain)]);
        if exist(domainDir, 'dir') == 0
            mkdir(domainDir);
        end
        saveas(fig, fullfile(domainDir, [base '.fig']));
        saveas(fig, fullfile(domainDir, [base '.png']));
        
        F = F+1;
        
        % update progress after each plot
        taskIdx = taskIdx + 1;
        elapsed = toc(t0);
        progress = taskIdx / totalTasks;
        if progress > 0
            estTotal = elapsed / progress;
            eta = estTotal - elapsed;
        else
            estTotal = NaN; eta = NaN;
        end
        disp(sprintf('Progress: %d/%d (%.1f%%) | phi=%d°, theta=%d°, domain=%d | elapsed=%.1fs, ETA=%.1fs', ...
            taskIdx, totalTasks, 100*progress, round(phi_deg), round(theta_deg), domain, elapsed, eta));

NF=F-1;

end
end

% for F = 1:NF
%     figure(F);
%     clim([wmin wmax]);
%     [i, d] = ind2sub([11, 8], F);
%     saveas(F,['xz_d' num2str(d-1) '_' num2str(i) '.png']);
%     close
% end


%% Functions
function J = Jfull(q, robot)
J1 = geometricJacobian(robot, q, 'tool_tcp0');
J = [J1(4:6, :); J1(1:3, :)];
end

function J = Jconstr(q, robot)
J1 = geometricJacobian(robot, q, 'rcm_link');
J2 = [J1(4:6, :); J1(1:3, :)];
J = J2(1:3, :);
end

function w = wconstr(q, P)
Jf = Jfull(q, P);
Jc = Jconstr(q, P);
N = eye(10) - pinv(Jc) * Jc;
Je = (N * Jf')';
w = sqrt(det(Je*Je'));
end

function w = wconstr1(q, P)
Jf = Jfull(q, P);
Jc = Jconstr(q, P);
N = eye(10) - pinv(Jc) * Jc;
Je = (N * Jf')';
Je1 = Je(:, [1:6 8:10]);
w = sqrt(det(Je1*Je1'));
end

function q = find_closest_solution(solutions, q0, threshold)
if nargin < 3
    threshold = pi/4;
end

n = size(solutions, 1);
D = zeros(n, 1);
for i = 1:n
    D(i) = norm(solutions(i, :) - q0);
end
[min_d, idx] = min(D);
if min_d < threshold
    q = solutions(idx, :);
else
    q = NaN;
end
end