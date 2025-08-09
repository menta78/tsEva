function u = tsCopulaRnd(family, copulaPar, N, uProb)

% generates the random vector from a copula distribution.
% if multivariate t or gaussian or bivariate archimedean, it uses copularnd
% else (this would be the case of multivariate archimedean) it implements a
% simple C-Vine copula
%
% parameters:
%   family: copula family
%   copulaPar: rho for gaussian, alpha for archimedean
%   N: sample size
%   uProb: data (in uniform probability space) from which the copula was built. Necessary in case the
%   c-vine copula must be constructed

if strcmpi(family, 't')
    error("copula family unsupported: t");
end

cndCopulaRnd0 = strcmpi(family, 'gaussian');
cndCopulaRnd1 = isscalar(copulaPar);
cndCopulaRnd = cndCopulaRnd0 | cndCopulaRnd1;

if cndCopulaRnd
    u = copularnd(family, copulaPar, N);
else
    if ~strcmpi(family, 'gumbel')
        error("multivariate archimedean only supported for gumbel");
    end
    alpha = copulaPar;
    order = cvineOrderFromCopulaPar(alpha);          % root-first ordering
    theta = fitGumbelCvine(uProb, alpha, order);   % estimate θ on all trees (pseudo-obs)
    u = simGumbelCvine(N, order, theta);      % sample (the sim I gave you)
end
end

function order = cvineOrderFromCopulaPar(alpha)
% Root-first C-vine order from a Gumbel-θ matrix (alpha).
% Heuristic: pick node with max total |tau| as root, then greedily add the next.
    d = size(alpha,1);
    if size(alpha,2) ~= d, error('alpha must be square'); end
    if any(diag(alpha) ~= 1), warning('alpha diagonal should be 1 for Gumbel'); end

    % convert to Kendall tau (τ = 1 - 1/θ)
    tau = 1 - 1./max(alpha, 1);  % guard θ>=1
    tau(1:d+1:end) = 0;
    absTau = abs(tau);

    scores = sum(absTau,2);
    picked = false(d,1);
    order  = zeros(1,d);

    [~, root] = max(scores);
    order(1) = root;
    picked(root) = true;

    for k = 2:d
        remaining = find(~picked);
        sc = absTau(remaining, order(1:k-1)) * ones(k-1,1);
        [~, ix] = max(sc);
        order(k) = remaining(ix);
        picked(order(k)) = true;
    end
end

function theta = fitGumbelCvine(uProb, alpha, order)
% Fit all Gumbel pair-copulas in a C-vine using pseudo-observations.
% Tree 1: take θ directly from alpha.
% Trees >=2: estimate τ from pseudo-obs, then θ = 1/(1-τ).

    [n,d] = size(uProb);
    if d ~= numel(order), error('uProb/order size mismatch'); end

    % reorder data by vine order
    U = uProb(:, order);
    Theta = cell(d-1,d);

    % running conditioned variables
    Ucond = U;

    % Tree 1: direct from alpha
    for j = 2:d
        Theta{1,j} = max(alpha(order(1), order(j)), 1); % θ>=1
    end

    % Trees 2..d-1
    for t = 2:d-1
        % build pseudo-obs for edges at tree t-1 so Ucond(:,j) = U_{v_j | v_1..v_{t-1}}
        for j = t:d
            % already have Ucond from previous loop when t>2
            % (nothing to do here before estimating; h() applied after using t)
        end

        % estimate θ on edges (t,j), j>t
        for j = t+1:d
            % Kendall tau on current pseudo-obs
            tau_tj = corr(Ucond(:,t), Ucond(:,j), 'type','Kendall', 'rows','pairwise');
            tau_tj = max(min(tau_tj, 0.999), 1e-6);
            Theta{t,j} = 1.0 / (1.0 - tau_tj); % θ = 1/(1-τ)
        end

        % update pseudo-obs for next tree (condition on pivot v_t)
        for j = t+1:d
            th = Theta{t,j};
            Ucond(:,j) = h_gumbel(Ucond(:,j), Ucond(:,t), th);
        end
    end

    theta = Theta;
end

function u = simGumbelCvine(N, order, Theta)
% Simulate from a Gumbel C-vine with parameters Theta in vine order.
    d = numel(order);
    U = zeros(N,d);
    W = rand(N,d);

    U(:,1) = W(:,1);
    for k = 2:d
        uk = W(:,k);
        for j = k-1:-1:1
            th = Theta{j,k};
            uk = hinv_gumbel(uk, U(:,j), th);
        end
        U(:,k) = uk;
    end

    % map back to original column order
    u = zeros(N,d);
    u(:, order) = U;
end

function h = h_gumbel(u, v, theta)
    epsv = 1e-12;
    u = min(max(u, epsv), 1-epsv);
    v = min(max(v, epsv), 1-epsv);
    theta = max(theta, 1);

    a = (-log(u)).^theta;
    b = (-log(v)).^theta;
    s = (a + b).^(1/theta);
    C = exp(-s);

    h = C .* (a + b).^(1/theta - 1) .* (-log(v)).^(theta - 1) ./ v;
    h = min(max(h, epsv), 1 - epsv);
end

function u = hinv_gumbel(z, v, theta)
    epsv = 1e-12;
    z = min(max(z, epsv), 1-epsv);
    v = min(max(v, epsv), 1-epsv);
    theta = max(theta, 1);

    n = numel(z);
    u = zeros(size(z));
    lo = epsv; hi = 1 - epsv;

    for i = 1:n
        zi = z(i); vi = v(i);
        f = @(uu) h_gumbel(uu, vi, theta) - zi;

        ui = min(max(zi, 1e-3), 1-1e-3); % starter
        ok = true;
        for it = 1:6
            du = 1e-6 * max(ui, 1e-3);
            g1 = f(min(max(ui - du, lo), hi));
            g2 = f(min(max(ui + du, lo), hi));
            der = (g2 - g1) / (2*du);
            val = f(ui);
            ui = min(max(ui - val / max(abs(der),1e-8), lo), hi);
            if abs(val) < 1e-8, break; end
            if ~isfinite(ui), ok = false; break; end
        end
        if ~ok || abs(f(ui)) > 1e-6
            try
                ui = fzero(f, [lo, hi]);
            catch
                ui = min(max(ui, lo), hi);
            end
        end
        u(i) = ui;
    end
end






