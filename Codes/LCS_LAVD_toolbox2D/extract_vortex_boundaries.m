function boundaries = extract_vortex_boundaries(xvec, yvec, F, cx, cy, d_max, n_levels)
% Extract the vortex boundary around each (cx(i), cy(i)) as the outermost
% closed, convex iso-contour of F that encloses that centre and no other.
%
% Inputs:
%   xvec, yvec : 1-D coordinate vectors; F is defined on meshgrid(xvec, yvec).
%   F          : scalar field (e.g. LAVD) of size [numel(yvec), numel(xvec)].
%   cx, cy     : vortex-center coordinates (column vectors, same length).
%   d_max      : convexity-deficiency tolerance (typical 1e-2 ... 1e-1).
%   n_levels   : number of iso-levels to sweep.
%
% Output:
%   boundaries : cell array of length numel(cx). boundaries{i} is an [N x 2]
%                matrix of (x,y) points defining the vortex boundary around
%                (cx(i), cy(i)), or [] if no admissible contour was found.

    num_c      = numel(cx);
    boundaries = cell(num_c, 1);

    F_hi = max(F(:), [], 'omitnan');
    F_lo = min(F(:), [], 'omitnan');
    if ~isfinite(F_hi) || F_hi <= F_lo
        return;
    end

    % Sweep from high (near the peak) to low (outward). Skip the degenerate
    % topmost level because it yields empty or single-point contours.
    levels = linspace(F_hi, F_lo, n_levels);
    levels(1) = [];

    for i = 1 : num_c
        best_poly  = [];
        best_level = Inf;
        other_ix   = setdiff(1:num_c, i);

        for L = levels
            C   = contourc(xvec, yvec, F, [L L]);
            cts = parse_contourc(C);

            for c = 1 : numel(cts)
                poly = cts{c};

                if ~is_closed(poly),                                     continue; end
                if ~inpolygon(cx(i), cy(i), poly(:,1), poly(:,2)),       continue; end
                if ~isempty(other_ix) && ...
                   any(inpolygon(cx(other_ix), cy(other_ix), ...
                                 poly(:,1), poly(:,2))),                 continue; end
                if convexity_deficiency(poly) > d_max,                   continue; end

                % Admissible. Keep the outermost (lowest-level) one.
                if L < best_level
                    best_poly  = poly;
                    best_level = L;
                end
            end
        end

        boundaries{i} = best_poly;
    end
end


% -----------------------------------------------------------------------
function contours = parse_contourc(C)
% Parse the packed output of contourc into a cell array of [N x 2] polylines.
    contours = {};
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        contours{end+1} = C(:, k+1 : k+n).';     %#ok<AGROW>
        k = k + n + 1;
    end
end


function tf = is_closed(poly)
% A contour is "closed" if its endpoints coincide within 1% of its bounding
% box scale. Also require a minimum point count to avoid degenerate shapes.
    if size(poly, 1) < 4
        tf = false; return;
    end
    scale = max( max(poly,[],1) - min(poly,[],1) );
    tf    = norm(poly(1,:) - poly(end,:)) < 0.01 * scale;
end


function d = convexity_deficiency(poly)
% Relative convexity deficiency:
%   d = (perimeter of contour - perimeter of its convex hull) / hull perimeter.
% 0 for a convex shape; grows with indentations.
    k      = convhull(poly(:,1), poly(:,2));
    hull   = poly(k, :);
    arc    = sum( sqrt( sum( diff(poly).^2, 2 ) ) );
    hull_a = sum( sqrt( sum( diff(hull).^2, 2 ) ) );
    d      = (arc - hull_a) / max(hull_a, eps);
end
