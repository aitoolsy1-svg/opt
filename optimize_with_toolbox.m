function optimize_with_toolbox(population_size, crossover_fraction, migration_fraction, migration_interval)
    % Main script to optimize RCM position using the MATLAB Global Optimization Toolbox.
    %
    % ARGS:
    %   population_size (optional): The number of individuals in the population. Default is 100.
    %   crossover_fraction (optional): The fraction of the next generation created by crossover. Default is 0.8.
    %   migration_fraction (optional): Fraction of individuals that migrate between subpopulations. Default is 0.2.
    %   migration_interval (optional): Number of generations between migrations. Default is 20.

    % --- Handle optional arguments ---
    if nargin < 4
        migration_interval = 20;
    end
    if nargin < 3
        migration_fraction = 0.2;
    end
    if nargin < 2
        crossover_fraction = 0.8;
    end
    if nargin < 1
        population_size = 100;
    end

    % --- Dependencies ---
    % This script requires the following:
    % - MATLAB Global Optimization Toolbox (for the 'ga' function).
    % - MATLAB Robotics System Toolbox.
    % - External functions: `ik_mops1.m`, `classify_ur5_domains.m`, `q2urdf.m`.
    % - MATLAB Parallel Computing Toolbox (recommended for speed).

    % --- Optimization Bounds ---
    % [X, Z, theta]
    nvars = 3;
    lb = [0.0, -0.5, -pi/2]; % Lower bounds
    ub = [1.0,  0.5,  pi/2]; % Upper bounds

    % --- Main Loop ---
    domains = 0:7;
    results = struct('domain', {}, 'optimal_params', {}, 'max_wr', {});

    % Load robot model once to be passed to the objective function
    disp('Loading robot model...');
    robot = importrobot('mops1.urdf', DataFormat='row');

    for i = 1:numel(domains)
        domain = domains(i);
        fprintf('Optimizing for domain %d...\n', domain);

        % Define the objective function for the current domain
        objFun = @(params) -objectiveFunction(params, domain, robot);

        % Configure GA options
        % For flat fitness landscapes, increasing exploration is key.
        % - Use migration to move individuals between subpopulations to escape local optima.
        % - Increase mutation by lowering the CrossoverFraction (e.g., to 0.6 or 0.7).
        %   The remaining fraction (e.g., 0.4 or 0.3) will be created by mutation.
        options = optimoptions('ga', ...
            'PopulationSize', population_size, ...
            'CrossoverFraction', crossover_fraction, ...
            'MigrationDirection', 'both', ...
            'MigrationInterval', migration_interval, ...
            'MigrationFraction', migration_fraction, ...
            'MaxStallGenerations', 10, ...
            'Display', 'iter', ...
            'PlotFcn', @gaplotbestf, ...
            'UseParallel', true);

        % Run the GA
        [optimal_params, max_wr_neg] = ga(objFun, nvars, [], [], [], [], lb, ub, [], options);

        % Store results
        results(i).domain = domain;
        results(i).optimal_params = optimal_params;
        results(i).max_wr = -max_wr_neg; % Convert back to positive for maximization

        fprintf('Domain %d: Optimal [X, Z, theta] = [%.4f, %.4f, %.4f], Max WR = %.6f\n', ...
            domain, optimal_params(1), optimal_params(2), optimal_params(3), -max_wr_neg);
    end

    % --- Save and Display Final Results ---
    fprintf('\n--- Optimization Complete ---\n');
    disp(results);
    save('toolbox_optimization_results.mat', 'results');
    fprintf('Results saved to toolbox_optimization_results.mat\n');
end

function fitness = objectiveFunction(params, domain, robot)
    % Objective function to be maximized (returns WR)
    % params: [X, Z, theta]

    X_rcm = params(1);
    Z_rcm = params(2);
    theta = params(3);

    % --- Configuration from test9 ---
    tcp_depth = -0.1;
    cube_res = 0.025; % Coarser resolution for faster optimization

    % Desired TCP pose (base)
    T_b_tcp_desired = [
        0   0 -1  0.5;
        0  -1  0  0.0;
        -1  0  0  tcp_depth;
        0   0  0  1
    ];
    T_b_tcp_desired(1:3, 1:3) = eul2rotm([0, -pi/2+1e-3, 0], 'zyx');

    % RCM pose from GA parameters
    T_b_rcm = [
        1 0 0  X_rcm;
        0 1 0  0;
        0 0 1  Z_rcm;
        0 0 0  1
    ];
    R_rcm = eul2rotm([0, theta, 0], 'zyx');
    T_b_rcm(1:3, 1:3) = R_rcm;

    % --- IK function ---
    ik_mops = @(TCP, RCM) ik_mops1(TCP, RCM);

    % --- Evaluate manipulability over a small cube ---
    [X, Y, Z] = ndgrid(-0.05:cube_res:0.05, -0.05:cube_res:0.05, -0.05:cube_res:0.05);
    [m, n, o] = size(X);
    W = zeros(m, n, o);

    for i1 = 1:m*n*o
        [i, j, k] = ind2sub([m, n, o], i1);
        T = T_b_tcp_desired;
        offset_rcm = [X(i, j, k); Y(i, j, k); Z(i, j, k)];
        offset_world = R_rcm * offset_rcm;

        T(1:4, 4) = T_b_rcm(1:4, 4) + [offset_world; 0];
        T(3, 4) = T(3, 4) - 0.2; % From test9 logic

        solutions = ik_mops(T, T_b_rcm);
        if isempty(solutions)
            W(i, j, k) = 0;
            continue;
        end

        D = classify_ur5_domains(solutions);
        ID = all([D(:, 1) == (bitand(domain, 4)>0) D(:, 2) == (bitand(domain, 2)>0) D(:, 3) == (bitand(domain, 1)>0)], 2);
        sol = solutions(ID, :);

        if isempty(sol)
            W(i, j, k) = 0;
            continue;
        end

        q = q2urdf(sol);
        W(i, j, k) = wconstr(q, robot);
    end

    fitness = mean(W, 'all');
    if isnan(fitness)
        fitness = 0;
    end
end

% --- Helper Functions ---

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
    if isempty(q)
        w = 0;
        return;
    end

    num_solutions = size(q, 1);
    total_w = 0;

    for i = 1:num_solutions
        q_single = q(i, :);
        Jf = Jfull(q_single, P);
        Jc = Jconstr(q_single, P);
        N = eye(10) - pinv(Jc) * Jc;
        Je = (N * Jf')';

        if rcond(Je*Je') < 1e-12
            total_w = total_w + 0;
        else
            total_w = total_w + sqrt(det(Je*Je'));
        end
    end

    w = total_w / num_solutions;
end