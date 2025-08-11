classdef tsGumbelCVine
% GUMBELCVINE  Utilities for a Gumbel C-vine built from pair-copulas.
% Static helpers:
%   order = tsGumbelCVine.cvineOrder(alpha)
%   Theta = tsGumbelCVine.fit(uProb, alpha, order)
%   U     = tsGumbelCVine.simulate(N, order, Theta)
%   h     = tsGumbelCVine.h(u, v, theta)
%   u     = tsGumbelCVine.hinv(z, v, theta)

methods(Static)
    
    function order = cvineOrder(alpha)
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

    function Theta = fit(uProb, alpha, order)
    % Fit all Gumbel pair-copulas in a C-vine using pseudo-observations.
    % Tree 1: take θ directly from alpha (matrix or scalar).
    % Trees >=2: estimate τ from pseudo-obs, then θ = 1/(1-τ).
        [n,d] = size(uProb);
        if d ~= numel(order), error('uProb/order size mismatch'); end
        if isscalar(alpha)
            A = alpha * ones(d);
        else
            A = alpha;
            if ~isequal(size(A), [d,d])
                error('alpha must be scalar or d-by-d matrix matching uProb cols (in original order).');
            end
            % reindex alpha to vine order
            A = A(order, order);
        end

        % reorder data by vine order
        U = uProb(:, order);

        Theta = cell(d-1, d);
        Ucond = U; % will be replaced by conditional pseudo-obs as we go

        % Tree 1: parameters from alpha (θ>=1)
        for j = 2:d
            Theta{1,j} = max(A(1,j), 1);
        end

        % Trees 2..d-1
        for t = 2:d-1
            % estimate θ on edges (t,j), j>t using current (t-1)-conditioned pseudo-obs
            for j = t+1:d
                tau_tj = corr(Ucond(:,t), Ucond(:,j), 'type','Kendall', 'rows','pairwise');
                tau_tj = max(min(tau_tj, 0.999), 1e-6);
                Theta{t,j} = 1.0 / (1.0 - tau_tj); % θ = 1/(1-τ), clamp implicit via tau
            end
            % update pseudo-obs for next tree by conditioning on pivot v_t
            for j = t+1:d
                th = Theta{t,j};
                Ucond(:,j) = tsGumbelCVine.h(Ucond(:,j), Ucond(:,t), th);
            end
        end
    end

    function U = simulate(N, order, Theta)
    % Simulate from a Gumbel C-vine with parameters Theta in vine order.
        d = numel(order);
        Uvine = zeros(N,d);
        W = rand(N,d);

        Uvine(:,1) = W(:,1);
        for k = 2:d
            uk = W(:,k);
            for j = k-1:-1:1
                th = Theta{j,k};
                uk = tsGumbelCVine.hinv(uk, Uvine(:,j), th);
            end
            Uvine(:,k) = uk;
        end

        % map back to original column order
        U = zeros(N,d);
        U(:, order) = Uvine;
    end

    function h = h(u, v, theta)
    % Gumbel h-function: C_{U|V}(u|v; theta) for pair-copula (U,V).
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

    function u = hinv(z, v, theta)
    % Inverse h-function: solve h(u,v,theta)=z for u in (0,1)
        epsv = 1e-12;
        z = min(max(z, epsv), 1-epsv);
        v = min(max(v, epsv), 1-epsv);
        theta = max(theta, 1);

        n = numel(z);
        u = zeros(size(z));
        lo = epsv; hi = 1 - epsv;

        for i = 1:n
            zi = z(i); vi = v(i);
            f = @(uu) tsGumbelCVine.h(uu, vi, theta) - zi;

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

end
end
