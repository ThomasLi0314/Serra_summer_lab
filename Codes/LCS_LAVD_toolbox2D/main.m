% LCS_LAVD_toolbox2D — driver script for the LAVD vortex-detection pipeline.
%
% Steps:
%   (1) initialize.m  : build the velocity field and grids (pick case_id).
%   (2) LAVD.m        : compute LAVD on a seed grid.
%   (3) local maxima  : identify vortex centers.
%   (4) extract       : outermost closed-convex iso-contour around each center.
%   (5) plot          : LAVD heatmap with centers and boundaries overlaid.

clear;
clc;
close all;


%% (1) Velocity field and grids
case_id = 3;
initialize;


%% (1b) Velocity-field video over [t_0, t_1]
make_video    = true;        % set false to skip
video_path    = sprintf('velocity_case%d.mp4', case_id);
frame_stride  = 5;           % render every Nth time step
quiver_stride = 10;          % subsample grid for arrow readability
fps           = 20;

if make_video
    fprintf('Writing velocity-field video to %s ...\n', video_path);

    speed_all = sqrt(u.^2 + v.^2);
    smax      = max(speed_all(:));
    if smax == 0, smax = 1; end

    Xq = X_grid(1:quiver_stride:end, 1:quiver_stride:end);
    Yq = Y_grid(1:quiver_stride:end, 1:quiver_stride:end);

    try
        vw = VideoWriter(video_path, 'MPEG-4');
    catch
        % Fallback for platforms without MPEG-4 support
        video_path = strrep(video_path, '.mp4', '.avi');
        vw = VideoWriter(video_path, 'Motion JPEG AVI');
    end
    vw.FrameRate = fps;
    open(vw);

    fig_v = figure('Color', 'w', 'Position', [100 100 720 620], 'Visible', 'off');
    ax_v  = axes(fig_v); %#ok<LAXES>

    frame_idx = 1 : frame_stride : Nt;
    for ii = 1 : numel(frame_idx)
        ti = frame_idx(ii);
        cla(ax_v);
        hold(ax_v, 'on');
        imagesc(ax_v, X_grid(1,:), Y_grid(:,1), speed_all(:,:,ti));
        set(ax_v, 'YDir', 'normal');
        caxis(ax_v, [0, smax]);
        colormap(ax_v, parula);
        cb = colorbar(ax_v); cb.Label.String = '|v|';

        Uq = u(1:quiver_stride:end, 1:quiver_stride:end, ti);
        Vq = v(1:quiver_stride:end, 1:quiver_stride:end, ti);
        quiver(ax_v, Xq, Yq, Uq, Vq, 1.2, 'Color', 'k', 'LineWidth', 0.5);

        axis(ax_v, 'equal');
        axis(ax_v, [-Lx Lx -Ly Ly]);
        xlabel(ax_v, 'x'); ylabel(ax_v, 'y');
        title(ax_v, sprintf('Velocity field — case %d, t = %.2f', case_id, T_grid(ti)));
        hold(ax_v, 'off');

        drawnow;
        writeVideo(vw, getframe(fig_v));
    end
    close(vw);
    close(fig_v);
    fprintf('   wrote %d frames.\n', numel(frame_idx));
end


%% Seed grid and LAVD computation, initial positions of tracer particle that gets advected through the flow. LAVD is only calculated at the Xseed Yseed points.
seed_skip = 1;                           
Xseed = X_grid(1:seed_skip:end, 1:seed_skip:end);
Yseed = Y_grid(1:seed_skip:end, 1:seed_skip:end);

t0_LAVD = t_0;
t1_LAVD = t_1;
n_steps = round((t1_LAVD - t0_LAVD) / dt);

fprintf('Computing LAVD on %dx%d seed grid over t in [%.2f, %.2f] ...\n', ...
        size(Xseed,1), size(Xseed,2), t0_LAVD, t1_LAVD);
tic;
LAVD_field = LAVD(X_grid, Y_grid, T_grid, u, v, Xseed, Yseed, ...
                  t0_LAVD, t1_LAVD, dx, dy, dt, t_0, n_steps);
fprintf('   done in %.2f s\n', toc);


%% (3) Vortex centers: local maxima above a fraction of the global max
max_frac_threshold = 0.3;

is_max = local_maxima_2d(LAVD_field);
L_max  = max(LAVD_field(:), [], 'omitnan');
is_max = is_max & (LAVD_field > max_frac_threshold * L_max) & ~isnan(LAVD_field);

[ci, cj] = find(is_max);
center_x = Xseed(sub2ind(size(Xseed), ci, cj));
center_y = Yseed(sub2ind(size(Yseed), ci, cj));
fprintf('Found %d candidate vortex center(s).\n', numel(center_x));


%% (4) Vortex boundaries: outermost closed convex iso-contours
d_max    = 1e-1;      % convexity-deficiency tolerance
n_levels = 80;        % number of iso-levels to sweep

xseed_vec = Xseed(1, :);
yseed_vec = Yseed(:, 1).';
boundaries = extract_vortex_boundaries(xseed_vec, yseed_vec, LAVD_field, ...
                                        center_x, center_y, d_max, n_levels);


%% (5) Plot LAVD field with vortex centers and boundaries
figure('Color', 'w', 'Position', [100 100 720 620]);
hold on;

contourf(Xseed, Yseed, LAVD_field, 30, 'LineStyle', 'none');
colormap(parula);
cb = colorbar;  cb.Label.String = 'LAVD';

% Faint reference contours
contour(Xseed, Yseed, LAVD_field, 12, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.3);

% Vortex centers (red stars)
if ~isempty(center_x)
    plot(center_x, center_y, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
end

% Vortex boundaries (red solid lines)
for k = 1 : numel(boundaries)
    if isempty(boundaries{k}), continue; end
    poly = boundaries{k};
    plot(poly(:,1), poly(:,2), 'r-', 'LineWidth', 2);
end

axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('LAVD — case %d (centers *, boundaries --)', case_id));
hold off;


% =======================================================================
% Local function: strict 3x3 local maximum (no Image Processing Toolbox).
% =======================================================================
function is_max = local_maxima_2d(F)
    pad                       = -inf(size(F) + 2);
    pad(2:end-1, 2:end-1)     = F;
    is_max                    = true(size(F));
    for di = -1:1
        for dj = -1:1
            if di == 0 && dj == 0, continue; end
            neighbor = pad(2+di : end-1+di, 2+dj : end-1+dj);
            is_max   = is_max & (F > neighbor);
        end
    end
end
