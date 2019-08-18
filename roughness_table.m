function [] = roughness_table(prefix, times, hours, ids, SmoothRadii, CurvatureRadii, AverageRadii)

N = 0;
for time_idx = 1:length(times)
    N = N + length(ids{time_idx});
end

M = 2 + length(SmoothRadii) * length(CurvatureRadii) * length(AverageRadii);

Ra = zeros(N, M);
Rq = zeros(N, M);
Rsk = zeros(N, M);
Rku = zeros(N, M);
row = 0;
for time_idx = 1:length(times)
    n = length(ids{time_idx});
    for id_idx = 1:n
        id = ids{time_idx}(id_idx);
        row = row + 1;
        Ra(row, 1) = id;
        Ra(row, 2) = hours(time_idx);
        Rq(row, 1) = id;
        Rq(row, 2) = hours(time_idx);
        Rsk(row, 1) = id;
        Rsk(row, 2) = hours(time_idx);
        Rku(row, 1) = id;
        Rku(row, 2) = hours(time_idx);
        if isempty(times{time_idx})
            name = sprintf('%s%d', prefix, id);
        else
            name = sprintf('%s%d_%s', prefix, id, times{time_idx});
        end
        col = 2;
        for SmoothRadius = SmoothRadii
            for CurvatureRadius = CurvatureRadii
                for AverageRadius = AverageRadii
                    col = col + 1;
                    [Ra(row, col), Rq(row, col), Rsk(row, col), Rku(row, col)] = roughness_params(name, SmoothRadius{1}, CurvatureRadius{1}, AverageRadius{1});
                end
            end
        end
    end
end

print_table(sprintf('results/Ra_%s.tex', prefix), Ra);
print_table(sprintf('results/Rq_%s.tex', prefix), Rq);
print_table(sprintf('results/Rsk_%s.tex', prefix), Rsk);
print_table(sprintf('results/Rku_%s.tex', prefix), Rku);

end