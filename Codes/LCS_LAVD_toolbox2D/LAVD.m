% Compute the Lagrangian-Averaged Vorticity Deviation (LAVD) at a list of
% initial positions (X0s, Y0s) over the window [t0, t1].
% The routine:
%   (1) advects particles from the initial seeds via partical_traj,
%   (2) obtains the spatial-mean vorticity bar{omega}(s) via mean_omega,
%   (3) interpolates omega_3 along each trajectory,
%   (4) integrates the absolute deviation in time (trapezoidal rule).

function LAVD = LAVD(X_grid, Y_grid, T_grid, u, v, X0s, Y0s, t0, t1, dx, dy, dt, t_0, n_steps)


    input_shape   = size(X0s);
    X0s           = X0s(:);
    Y0s           = Y0s(:);
    num_particles = numel(X0s);

    % Get particle trajectories
    [X_traj, Y_traj] = partical_traj(X_grid, Y_grid, T_grid, u, v, X0s, Y0s, t_0, dt, n_steps);

    % Return mean omega as a function of time
    mean_omega_t = mean_omega(u, v, t0, t1, dx, dy);

    % Calculate the vorticity field
    Nt_field = size(u, 3);
    omega3   = zeros(size(u));
    for i = 1 : Nt_field
        [~,dudy_t] = gradient(u(:, :, i), dx, dy);
        [dvdx_t,~] = gradient(v(:, :, i), dx, dy);
        omega3(:, :, i)  = dvdx_t - dudy_t;
    end

    % Interpolate the vorticity field
    x_vec      = X_grid(1, :);
    y_vec      = Y_grid(:, 1).';
    omega3_xyt = permute(omega3, [2, 1, 3]);
    F_omega    = griddedInterpolant({x_vec, y_vec, T_grid}, omega3_xyt, 'linear', 'none');

    % Calculate the integrand of LAVD, |omega_3(X(s), s) - bar{omega}(s)| 
    traj_t    = t_0 + (0 : n_steps) * dt;
    integrand = zeros(num_particles, n_steps + 1);

    for k = 1 : (n_steps + 1)
        t_k           = traj_t(k);
        t_k_vec       = repmat(t_k, size(X_traj(:, k)));
        omega_at_traj = F_omega(X_traj(:, k), Y_traj(:, k), t_k_vec);
        mean_at_tk    = interp1(T_grid, mean_omega_t, t_k, 'linear', 'extrap');
        integrand(:, k) = abs(omega_at_traj - mean_at_tk);
    end

    % Integrate to get the LAVD
    LAVD_flat = trapz(traj_t, integrand, 2); 
    LAVD = reshape(LAVD_flat, input_shape);

end
