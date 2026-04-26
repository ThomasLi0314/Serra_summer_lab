% This function computes the mean vorticity at a given 

function mean_omega = mean_omega(u, v, t0, t1, dx, dy)
    % Input : 
    % u : x-component of velocity field on a spatial domain U(t)
    % v : y-component of velocity field on a spatial domain U(t)
    % t0, t1 : time interval
    % dx, dy : grid spacing
    % Output : 
    % mean_xi : mean vorticity on U(t) as a function of time

    % The velocity field has a boundary, outside the boundary the matrix entry is 0. 
    
    total_t_steps = size(u, 3);
    t = linspace(t0, t1, total_t_steps);

    % Initialize storage
    mean_omega = zeros(total_t_steps, 1);

    for i = 1 : total_t_steps
        u_t = u(:, :, i);
        v_t = v(:, :, i);

        is_outside_t = (abs(u_t) < 1e-10) & (abs(v_t) < 1e-10);
        is_interior  = ~is_outside_t;

        [~, dudy_t] = gradient(u_t, dx, dy);
        [dvdx_t, ~] = gradient(v_t, dx, dy);
        omega_t = dvdx_t - dudy_t;

        % Exclude the NaN-free interior only (one-cell erosion of the mask)
        valid = is_interior & ~isnan(omega_t);
        mean_omega(i) = sum(omega_t(valid)) / nnz(valid);
    end
end
    