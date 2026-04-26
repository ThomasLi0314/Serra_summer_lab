% Obtain the partical trajectory of a set of initial position X0 from teh velocity field u and v

function [X_traj, Y_traj] = partical_traj(X_grid, Y_grid, T_grid, u, v, X0, Y0, t_0, dt, n_steps)

    num_particles = length(X0);

    % Interpolate the velocity field, other options include:
    % Changing 'linear' to 'cubic' or 'spline' for higher-order interpolation
    x_vec = X_grid(1, :);    
    y_vec = Y_grid(:, 1).';      
    % Permute so first dim = x, second = y, third = t (matches F_u(x,y,t) calls)
    u_xyt = permute(u, [2, 1, 3]);
    v_xyt = permute(v, [2, 1, 3]);

    F_u = griddedInterpolant({x_vec, y_vec, T_grid}, u_xyt, 'linear', 'none');
    F_v = griddedInterpolant({x_vec, y_vec, T_grid}, v_xyt, 'linear', 'none');

    % Preallocate trajectory arrays for speed
    X_traj = zeros(num_particles, n_steps + 1);
    Y_traj = zeros(num_particles, n_steps + 1);

    % Initial position 
    X_traj(:, 1) = X0;
    Y_traj(:, 1) = Y0;

    current_t = t_0;

    % Implement a RK4 scheme here
    for step = 1 : n_steps 
        x_n = X_traj(:, step);
        y_n = Y_traj(:, step);
        t_n = repmat(current_t, size(x_n));

        % Step 1
        u1 = F_u(x_n, y_n, t_n);
        v1 = F_v(x_n, y_n, t_n);

        % March to step 2
        x2 = x_n + 0.5 * u1 * dt;
        y2 = y_n + 0.5 * v1 * dt;
        t2 = t_n + 0.5 * dt;
        u2 = F_u(x2, y2, t2);
        v2 = F_v(x2, y2, t2);

        % March to step 3
        x3 = x_n + 0.5 * dt * u2;
        y3 = y_n + 0.5 * dt * v2;
        t3 = t_n + 0.5 * dt;
        u3 = F_u(x3, y3, t3);
        v3 = F_v(x3, y3, t3);
        
        % March to step 4
        x4 = x_n + dt * u3;
        y4 = y_n + dt * v3;
        t4 = t_n + dt;
        u4 = F_u(x4, y4, t4);
        v4 = F_v(x4, y4, t4);

        % Update positions
        X_traj(:, step+1) = x_n + (dt / 6) * (u1 + 2*u2 + 2*u3 + u4);
        Y_traj(:, step+1) = y_n + (dt / 6) * (v1 + 2*v2 + 2*v3 + v4);

        % update time 
        current_t = current_t + dt;
    end
end

        

        
