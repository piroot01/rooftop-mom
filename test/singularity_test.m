clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

a = 0:0.01:10;
b = 0:0.02:20;

sings1 = nan(1, numel(a));
sings2 = nan(1, numel(b));
for n = 1:numel(a)
    sings1(n) = mom.singularity.singularity_1(a(n), b(n));
    sings2(n) = mom.singularity.singularity_5(a(n), b(n));
end

figure;
loglog(a, sings1);
hold on;
loglog(a, sings2);