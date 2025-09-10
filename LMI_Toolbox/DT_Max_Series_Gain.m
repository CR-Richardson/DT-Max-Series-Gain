
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
% Note: Only Circle and Circle-like Criteria can handle D~=0. All other 
% criteria will return nan. 
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
% Total_Ex:       Total number of examples
% Ex_array:       The selected set of example systems to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex   = 8;
Ex_array   = 1:8;

%% Makes example systems accessible to script
DT_Examples;

%% Calculate maximum series gain using various criteria

All_ex = 1:Total_Ex; % All example systems

% Arrays for storing the maximum series gain (alpha) for each example.
Alpha_Circle        = zeros(All_ex);
Alpha_Circle_Like   = zeros(All_ex);
Alpha_Popov         = zeros(All_ex);
Alpha_Popov_Like_H  = zeros(All_ex);
Alpha_Park          = zeros(All_ex);

% Arrays for storing the # of decision variables for each example.
decs_Circle        = zeros(All_ex);
decs_Circle_Like   = zeros(All_ex);
decs_Popov         = zeros(All_ex);
decs_Popov_Like_H  = zeros(All_ex);
decs_Park          = zeros(All_ex);

for i=Ex_array
    disp(['Example ',num2str(i),' ']);
    
    disp('DT Circle calculations...'); 
    [Alpha_Circle(i), data1(i), decs_Circle(i)] = DT_Circle(Syst{i});
    
    disp('DT Circle-like calculations...'); 
    [Alpha_Circle_Like(i), data2(i), decs_Circle_Like(i)] = DT_Circle_Like(Syst{i});

    disp('DT Popov calculations...'); 
    [Alpha_Popov(i), data3(i), decs_Popov(i)] = DT_Popov(Syst{i});
    
    disp('DT Popov-like (H relax) calculations...'); 
    [Alpha_Popov_Like_H(i), data4(i), decs_Popov_Like_H(i)] = DT_Popov_Like(Syst{i});

    disp('DT Park calculations...'); 
    [Alpha_Park(i), data5(i), decs_Park(i)] = DT_Park(Syst{i});
    
end

%% Display max series gain
disp(' ');
disp('Max. series gain');
title_str=['        Example', '         Circle', '    Circle-Like', '          Popov', ...
           '   Popov-Like (H)', '           Park'];
mat_data =[Ex_array' Alpha_Circle(Ex_array)' Alpha_Circle_Like(Ex_array)' Alpha_Popov(Ex_array)' ...
           Alpha_Popov_Like_H(Ex_array)' Alpha_Park(Ex_array)']; 
           
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
