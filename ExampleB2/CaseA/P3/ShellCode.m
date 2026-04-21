clc; clear; close all

no_sim = 1000;
n = [200, 300, 400];

for i = 1:length(n)
    RunSim(n(i), no_sim)
end