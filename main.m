function [] = main(SmoothRadius, CurvatureRadius, AverageRadius, log_scale)
twip_times = {'', 't5'};
twip_time_labels = {'0 h', '5 h'};
twip_hours = [0, 5];
twip_indices = {1:5, 1:5};
roughness_stats('twip', twip_times, twip_time_labels, twip_indices, SmoothRadius, CurvatureRadius, AverageRadius, log_scale);
purefe_times = {'', 't5', 't24', 't7'};
purefe_time_labels = {'0 h', '5 h', '1 day', '7 days'};
purefe_hours = [0, 5, 24, 168];
purefe_indices = {1:15, 11:15, 6:10, 1:5};
roughness_stats('purefe', purefe_times, purefe_time_labels, purefe_indices, SmoothRadius, CurvatureRadius, AverageRadius, log_scale);
end
