% Initialize the velocity field and build x, y, t grids.
%
% Sets (in the caller workspace):
%   X_grid, Y_grid, T_grid, u, v         : field and grids
%   Lx, Ly, Nx, Ny, dx, dy               : spatial geometry
%   t_0, t_1, dt, Nt                     : time geometry
%   case_id                              : which example was used


%% Choose which example to use
%   1 : Point vortex         (irrotational everywhere except origin)
%   2 : Lamb-Oseen vortex    (rotational, Gaussian vorticity core)
%   3 : Multi Lamb-Oseen     (four well-separated Gaussian-core vortices)
if ~exist('case_id', 'var')
    case_id = 2;
end


%% Shared spatial setup
Lx = 5;
Ly = 5;
Nx = 201;
Ny = 201;
dx = 2*Lx / (Nx - 1);       % grid covers [-Lx, Lx]
dy = 2*Ly / (Ny - 1);
[X_grid, Y_grid] = meshgrid(-Lx:dx:Lx, -Ly:dy:Ly);


%% Shared time setup
t_0 = 0;
t_1 = 10;
dt  = 0.01;
T_grid = t_0 : dt : t_1;
Nt     = numel(T_grid);

R2 = X_grid.^2 + Y_grid.^2;


%% Case-specific velocity field (steady)
switch case_id

    case 1
        % ---- Case 1: Point vortex (irrotational except at origin) ----
        %   v = [ -alpha*y / r^2 ,  alpha*x / r^2 ]
        %   Flow defined only on the disk { x^2 + y^2 < 2*Lx/3 }.
        alpha  = 1;
        mask   = R2 < (2*Lx/3);
        R2_reg = max(R2, eps);

        u_steady = -alpha * Y_grid ./ R2_reg;
        v_steady =  alpha * X_grid ./ R2_reg;
        u_steady(~mask) = 0;
        v_steady(~mask) = 0;

    case 2
        % ---- Case 2: Rotational flow
        % v = [ -alpha * y, alpha * x
        % Similarly defineed on the disk

        alpha = 1;
        mask = R2 < (2 * Lx / 3);

        u_steady = -alpha * Y_grid;
        v_steady = alpha * X_grid;
        u_steady(~mask) = 0;
        v_steady(~mask) = 0;

    case 3
        % ---- Case 3: Superposition of Lamb-Oseen vortices, moving ----
        % Each vortex k contributes a tangential velocity
        %   u_theta(r) = (Gamma_k / (2*pi*r)) * (1 - exp(-r^2 / r_c^2))
        % about a *time-dependent* center (xc_k(t), yc_k(t)).
        % Centers orbit a slow common circle plus small individual wobbles
        % so the four cores remain well separated for all t in [t_0, t_1].
        base_centers = [ -2.0, -2.0;
                          2.0, -2.0;
                         -2.0,  2.0;
                          2.0,  2.0 ];
        Gamma  = [  6.0,  -6.0,  -6.0,   6.0 ];   % alternating signs
        r_c    = 0.6;

        Omega    = 2*pi / (t_1 - t_0);   % one full orbit over [t_0, t_1]
        A_orbit  = 0.6;                  % radius of the common slow rotation
        A_wobble = 0.3;                  % per-vortex wobble amplitude
        Omega_w  = 2 * Omega;            % wobble frequency
        Nv       = size(base_centers, 1);

        u = zeros(Ny, Nx, Nt);
        v = zeros(Ny, Nx, Nt);
        for ti = 1 : Nt
            t = T_grid(ti);
            U_t = zeros(size(X_grid));
            V_t = zeros(size(X_grid));
            for k = 1 : Nv
                phase_k = 2*pi*(k-1)/Nv;
                xc = base_centers(k,1) ...
                     + A_orbit  * cos(Omega*(t - t_0) + phase_k) ...
                     + A_wobble * sin(Omega_w*(t - t_0) + phase_k);
                yc = base_centers(k,2) ...
                     + A_orbit  * sin(Omega*(t - t_0) + phase_k) ...
                     + A_wobble * cos(Omega_w*(t - t_0) + phase_k);

                dxk  = X_grid - xc;
                dyk  = Y_grid - yc;
                r2k  = dxk.^2 + dyk.^2;
                r2r  = max(r2k, eps);
                core = (Gamma(k) / (2*pi)) .* (1 - exp(-r2k / r_c^2)) ./ r2r;
                U_t  = U_t - core .* dyk;
                V_t  = V_t + core .* dxk;
            end
            u(:,:,ti) = U_t;
            v(:,:,ti) = V_t;
        end
end


%% Broadcast to time for steady cases
if case_id ~= 3
    u = repmat(u_steady, [1, 1, Nt]);
    v = repmat(v_steady, [1, 1, Nt]);
end


%% Load Field Data
