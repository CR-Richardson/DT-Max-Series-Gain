
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% CR Richardson and MC Turner 
% ECS
% University of Southampton
% UK
%
% Date: 26/11/23
%
% Purpose:
% This script builds various discrete time linear systems and computes the
% maximum series gain (alpha) according to various criteria for which the
% Lurie system is stable when the repeated ReLU is placed in the feebdack path.  
%
% Scripts
% DT_Examples:    Contains example discrete-time linear systems
%
% Functions
% DT_Circle:      Discrete time Circle Criterion - See (Haddad and Bernstein, 1994)
% DT_Circle_Like: Discrete time Circle-Like Criterion - See Theorem 12
% DT_Popov:       Discrete time Popov Criterion - See (Haddad and Bernstein, 1994)
% DT_Popov_Like:  Discrete time Relaxed Popov-Like Criterion - See Remark 15
% DT_Park:        Discrete time Park Criterion - See (Park, 2019)

% Variables
% Total_Ex: Total number of examples
% Ex_array: The selected set of example systems to test
% eps:      Loop termination accuracy - try 1e-6 for Examples 1-8 and 1e-2 for Examples 9-12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex = 12;
Ex_array = 1:8;
eps      = 1e-6;

%% Makes example systems accessible to script
DT_Examples;

%% Nyquist gains of each example
alpha_up = [6666.6651, 89.9, 0.6983, 0.0020, 0.0869, 0.8212, ...
            0.2008, 2.0221, 2.0600, 2.7000, 2.2200, 2.7500];

alpha_up = alpha_up + 0.1; % ensures alpha_low != alpha_up for each example.

%% Arrays for storing the maximum series gain (alpha) and # of decision variables.

All_ex = 1:Total_Ex; % All example systems

% Arrays for storing the maximum series gain (alpha) for each example.
alpha_Circle        = zeros(All_ex);
alpha_Circle_Like   = zeros(All_ex);
alpha_Popov         = zeros(All_ex);
alpha_Popov_Like_H  = zeros(All_ex);
alpha_Park          = zeros(All_ex);

% Arrays for storing the # of decision variables for each example.
decs_Circle        = zeros(All_ex);
decs_Circle_Like   = zeros(All_ex);
decs_Popov         = zeros(All_ex);
decs_Popov_Like_H  = zeros(All_ex);
decs_Park          = zeros(All_ex);

%% Calculate maximum series gain using various criteria

for i=Ex_array
    disp(['Example ',num2str(i),' ']);
    
    disp('DT Circle calculations...');
    tic;
    [alpha_Circle(i), data_Circle(i), decs_Circle(i)] = DT_Circle(Syst{i}, eps, alpha_up(i));
    toc;

    disp('DT Circle-like calculations...'); 
    tic;
    [alpha_Circle_Like(i), data_Circle_like(i), decs_Circle_Like(i)] = DT_Circle_Like(Syst{i}, eps, alpha_up(i), alpha_Circle(i));
    toc;

    disp('DT Popov calculations...');
    tic;
    [alpha_Popov(i), data_Popov(i), decs_Popov(i)] = DT_Popov(Syst{i}, eps, alpha_up(i), alpha_Circle(i));
    toc;

    disp('DT Popov-like (H relax) calculations...');
    tic;
    [alpha_Popov_Like_H(i), data_Popov_Like_H(i), decs_Popov_Like_H(i)] = DT_Popov_Like(Syst{i}, eps, alpha_up(i), alpha_Popov(i));
    toc;

    disp('DT Park calculations...');
    disp('Very slow for high-dimensional systems.');
    tic;
    [alpha_Park(i), data_Park(i), decs_Park(i)] = DT_Park(Syst{i}, eps, alpha_up(i), alpha_Circle(i));
    toc;
    
end

%% Display max series gain
disp(' ');
disp('Max. series gain');
title_str=['        Example', '         Circle', '    Circle-Like', '          Popov', ...
           '   Popov-Like (H)', '           Park'];
mat_data =[Ex_array' alpha_Circle(Ex_array)' alpha_Circle_Like(Ex_array)' alpha_Popov(Ex_array)' ...
           alpha_Popov_Like_H(Ex_array)' alpha_Park(Ex_array)']; 
           
fprintf('%15s %15s %15s %15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.4f %14.4f %14.4f %16.4f %14.4f\n',mat_data');

%% Display # of decision variables
disp(' ');
disp('# of decision variables');
title_str=['        Example', '         Circle', '    Circle-Like', '          Popov', ...
           '   Popov-Like (H)', '           Park'];
mat_data =[Ex_array' decs_Circle(Ex_array)' decs_Circle_Like(Ex_array)' decs_Popov(Ex_array)' ...
           decs_Popov_Like_H(Ex_array)' decs_Park(Ex_array)'];

fprintf('%15s %15s %15s %15s %15s %15s \n',title_str);
disp(' ');
fprintf('%15d %14.0f %14.0f %14.0f %16.0f %14.0f\n',mat_data');
