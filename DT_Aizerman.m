
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% MC Turner
% ECS
% University of Southampton
% UK
%
% Date: 15/05/24
%
% Purpose: 
% Computes the Nyquist gain (alpha) for the example discrete-time systems
% based on the Aizerman conjecture. The algorithm assumes the set of
% stabilising gains is connected, which is often, but not always the case.
%
% Scripts
% DT_Examples: Contains example discrete-time linear systems (Syst{i})
%
% Variables:
% alpha: Nyquist gain (float)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex   = 8;
Ex_array   = 1:8;
eps        = 1e-6;

%% Makes example systems accessible to script
DT_Examples;

%% Find largest slope size for each example

alpha = zeros(1,Total_Ex);

for i=Ex_array
    % calculate initial upper and lower slope size
    curlyA_u=10000;         % Upper value of slope size
    curlyA_l=0;             % Lower value of slope size (we know 0 is always feasible as system is stable)
    curlyA=curlyA_u;

    % Load system
    A     = Syst{i}.a;
    B     = Syst{i}.b;
    C     = Syst{i}.c;
    D     = Syst{i}.d;
    [~,m] = size(B);
    
    % Determine largest slope size (curlyA)
    while ((curlyA_u-curlyA_l)/curlyA_u) > eps

        Acl = A+B*curlyA*(inv(eye(m)-D*curlyA))*C;
        E   = eig(Acl);
        rho = max(abs(E));

        if rho < 1                  % System is stable
            curlyA_l=curlyA;
        else 
            curlyA_u=curlyA;        % System is unstable
        end
  
        curlyA=(curlyA_u+curlyA_l)/2;

    end

    alpha(i) = curlyA_l;

end

