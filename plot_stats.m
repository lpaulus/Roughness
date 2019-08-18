function [] = plot_stats(fig, x, time_idx)
% We know that (mean - mu) / sqrt(S^2 / n) follows a student distribution
% with degree of freedom n - 1.
% We want a confidence interval of 95% so we take t so that 2.5% is higher
% than t so that [-t, t] has 95% of the mass of the distribution.
% The confidence interval is therefore
% [mean - t * sqrt(S^2 / n), mean + t * sqrt(S^2 / n)]
m = mean(x);
s2 = var(x);
n = length(x);
t = tinv(0.975, n - 1);
err = t * sqrt(s2 / n);
hold(fig, 'on');
bar = errorbar(fig, time_idx, m, err, 'LineWidth', 2, 'CapSize', 30, 'Marker', '+', 'MarkerSize', 15);
%bar.XNegativeDelta = min(m - 1e-6, bar.XNegativeDelta);
scatter(fig, time_idx * ones(n, 1), x, 64, 'Marker', 'x', 'LineWidth', 2);
end