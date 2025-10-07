%clc; clear all; close all;
function [XR, ZR, WR]=test9_rcm_xz_domains_fun_average_adapter1(domain)
%% Configuration
tcp_depth = -0.1;
res = 0.01;
cube_res = 0.025;
adapter = 1;

% this is the desired TCP pose
T_b_tcp_desired = [
    0   0 -1  0.5
    0  -1  0  0.0
    -1  0  0  tcp_depth
    0   0  0  1
];

T_b_tcp_desired(1:3, 1:3) = eul2rotm([0, -pi/2+1e-3, 0], 'zyx'); % !!!change the zero here to small value like 0.001 to avoid singularity!!!

% this is the RCM pose
T_b_rcm = [
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
F = 1;
wmax = 0.0;

%for domain = 0:7

[XR, ZR] = meshgrid(0.0:res:1.0, -0.5:res:0.5);
%XR=XR+0.000001;
%ZR=ZR+0.000001;
[M, N] = size(XR);
WR = zeros(M, N);

dt = 0;
t_elapsed = 0;
cnt = 0;
for i2 = 1 : M*N
    tic;
    [I, J] = ind2sub([M, N], i2);
    t_est = ((M*N) / (i2)) * t_elapsed - t_elapsed;

    % if ~mod(cnt, 100)
    %     disp([num2str(domain) ' ' num2str(100*((J-1)*N+I)/(M*N)) '% ' num2str(t_elapsed) 's / ' num2str(t_est) 's'])
    % end
    cnt = cnt+1;

    TR = T_b_rcm;
    TR(1, 4) = XR(I, J);
    TR(3, 4) = ZR(I, J);

    [X, Y, Z] = ndgrid(-0.05:cube_res:0.05, -0.05:cube_res:0.05, -0.05:cube_res:0.05);
    [m, n, o] = size(X);
    W = zeros(m, n, o);
    q0 = NaN;
    
    for i1 = 1 : m*n*o
        [i, j, k] = ind2sub([m, n, o], i1);
        % disp([num2str(i) ' ' num2str(j) '' num2str(k)])
        T = T_b_tcp_desired;
        T(1, 4) = TR(1, 4) + X(i, j, k) + 0.0;
        T(2, 4) = TR(2, 4) + Y(i, j, k) + 0.0;
        T(3, 4) = TR(3, 4) + Z(i, j, k) - 0.2;

        % T, TR
    
        solutions = ik_mops(T, TR);
        if numel(solutions) < 1
            W(i, j, k) = 0;
            continue
        end
        
        % use only selected postures
        D = classify_ur5_domains(solutions);
        ID = all([D(:, 1) == (bitand(domain, 4)>0) D(:, 2) == (bitand(domain, 2)>0) D(:, 3) == (bitand(domain, 1)>0)], 2);
        solutions = solutions(ID, :);
        sol = solutions;
    
        if numel(solutions) < 1
            W(i, j, k) = 0;
            continue
        end
    
        q = q2urdf(sol);
        W(i, j, k) = wconstr(q, robot);
    end

    % W
 %   WR(I, J) = min(min(min(W)));
 WR(I, J)=mean(W,"all");
  %  if WR(I, J) > wmax
   %     wmax = WR(I, J);
  %  end


    dt = toc;
    t_elapsed = t_elapsed + dt;
end

% figure(F);
% % subplot(2, 4, domain+1);
% h = pcolor(XR, ZR, abs(WR))
% set(h, 'EdgeColor', 'none')
% xlabel('x')
% ylabel('z')
% % zlabel('w(q)')
% title(['min w(q) for RCM on X-Z plane, domain=' num2str(domain)]) % num2str(tcp_depth)])
% % view([45 45])
% 
% F = F + 1;
end

% NF = F-1;
% 
% for F = 1:NF
%     figure(F);
%     clim([0 wmax]);
%     % saveas(F,['xz_y' num2str(F) '_d' num2str(domain) '.png']);
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
