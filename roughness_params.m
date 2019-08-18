function [Ra,Rq,Rsk,Rku] = roughness_params(name, SmoothRadius, CurvatureRadius, AverageRadius)
file = sprintf('results/%s/%s/%s/%s/roughness.txt', name, SmoothRadius, CurvatureRadius, AverageRadius);
y = load(file);
Ra = mean(abs(y));
var_y = var(y, 1);
Rq = sqrt(var_y);
Rsk = skewness(y, 1);
Rku = kurtosis(y, 1);
end