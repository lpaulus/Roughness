function [] = roughness_stats(prefix, times, time_labels, ids, SmoothRadius, CurvatureRadius, AverageRadius)
figure('Name', 'Ra');
Ra_fig = axes();
xlim(Ra_fig, [0, length(times) + 1]);
set(Ra_fig,'xtick',1:length(time_labels),'xticklabel',time_labels);

figure('Name', 'Rq');
Rq_fig = axes();
xlim(Rq_fig, [0, length(times) + 1]);
set(Rq_fig,'xtick',1:length(time_labels),'xticklabel',time_labels);

figure('Name', 'Rsk');
Rsk_fig = axes();
xlim(Rsk_fig, [0, length(times) + 1]);
set(Rsk_fig,'xtick',1:length(time_labels),'xticklabel',time_labels);

figure('Name', 'Rku');
Rku_fig = axes();
xlim(Rku_fig, [0, length(times) + 1]);
set(Rku_fig,'xtick',1:length(time_labels),'xticklabel',time_labels);

for time_idx = 1:length(times)
    n = length(ids{time_idx});
    Ra = zeros(n, 1);
    Rq = zeros(n, 1);
    Rsk = zeros(n, 1);
    Rku = zeros(n, 1);
    for id_idx = 1:n
        id = ids{time_idx}(id_idx);
        if isempty(times{time_idx})
            name = sprintf('%s%d', prefix, id);
        else
            name = sprintf('%s%d_%s', prefix, id, times{time_idx});
        end
        [Ra(id_idx), Rq(id_idx), Rsk(id_idx), Rku(id_idx)] = roughness_params(name, SmoothRadius, CurvatureRadius, AverageRadius);
    end
    plot_stats(Ra_fig, Ra, time_idx);
    plot_stats(Rq_fig, Rq, time_idx);
    plot_stats(Rsk_fig, Rsk, time_idx);
    plot_stats(Rku_fig, Rku, time_idx);
end

suffix = sprintf('%s_%s_%s_%s', prefix, SmoothRadius, CurvatureRadius, AverageRadius);
saveas(Ra_fig, sprintf('results/Ra_%s.eps', suffix), 'epsc');
saveas(Rq_fig, sprintf('results/Rq_%s.eps', suffix), 'epsc');
saveas(Rsk_fig, sprintf('results/Rsk_%s.eps', suffix), 'epsc');
saveas(Rku_fig, sprintf('results/Rku_%s.eps', suffix), 'epsc');

end