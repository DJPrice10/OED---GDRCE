function [data_SEKI_inf]=getdataSEKI(transmission_rate, exposure_rate, N, e1, sample_times, nclasses)
% Simulate the SEKI model and get infected no. at specified times
% transmission_rate, effective rate of transmission
% expsure_rate, rate of moving through exposed classes
% N, population size
% e1, initial state (no. exposed)
% Tmax, maximum time, end of observing epidemic
% sample_times, observation times of process

% transmission_rate=2;
% exposure_rate=1;
% N=20;
% e1=10;
% sample_times=[1:10];

Tmax=max(sample_times);

[n_infectious, time]=SEKI_simulate(transmission_rate, exposure_rate, N, e1, Tmax, nclasses);
data_SEKI_inf = infected_pop_at_time(sample_times, n_infectious, time);
