function [] = print_table(name, T)
file = fopen(name, 'w');
T = sortrows(T);
for row = 1:size(T, 1)
    fprintf(file, '%d & %d', T(row, 1), T(row, 2));
    for col = 3:size(T, 2)
        fprintf(file, ' & %f', T(row, col));
    end
    fprintf(file, '\\\\\n');
    if mod(row, 2) == 0
        fprintf(file, '\\hline');
    end
end
fclose(file);

end