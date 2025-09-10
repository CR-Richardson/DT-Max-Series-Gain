% ECS
% University of Southampton
% UK
%
% Date: 5/09/25
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the DT Popov like criterion
% as defined by Remark 15.
%
% Parameters:
% syst:      Structure containing the system matrices of an example.
% eps:       Loop termination accuracy.
% alpha_up:  Largest series gain to be tested.
% alpha_low: Lowest series gain to be tested.
%
% Returns:
% alpha_low: Maximum series gain.
% data:      Structure containing solutions of the LMI parametrised by alpha.
% dec:       # of decision variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha_low,data,dec]=DT_Popov_Like(syst,eps,alpha_up,alpha_low)

%% Parameters
A     = syst.a;
B     = syst.b;
C     = syst.c;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

% Raise error if D != 0.
if ~isequal(D, zeros(m))
    error('D terms must be added to LMIs.');
end

%% Initialising alpha
alpha = alpha_up;

%% Determine max alpha by repeatedly solving LMI for largest feasible alpha

while ((alpha_up - alpha_low)/alpha_up) > eps

% define variables
P    = sdpvar(n,n,'symmetric');
L    = sdpvar(m,m,'diagonal');
V    = sdpvar(m,m,'full');
Q11  = sdpvar(m,m,'full');

% define constraints
epsilon = 1e-3;
In = eye(n);
Im = eye(m);
c1 = P >= epsilon*In; % PD
c2 = L >= epsilon*Im; % PD
c3 = V(~Im) <= 0; % Z-matrix
c4 = Q11(:) >= 0; % Non-negative matrix
% c5 = Q11_(:) >= 0; % Non-negative matrix

X11 = A'*P*A - P + ( (A - In)'*C'*L*C*(A - In) ) + ( (A - In)'*C'*L*C*(A - In) )';
X12 = alpha*A'*P*B + C'*V' + (A - In)'*C'*L + 2*alpha*(A - In)'*C'*L*C*B;
X22 = alpha*alpha*B'*P*B + ( Q11 - V + alpha*L*C*B + alpha*alpha*B'*C'*L*C*B ) + ( Q11 - V + alpha*L*C*B + alpha*alpha*B'*C'*L*C*B )';

c5 = [X11,  X12;
      X12', X22] <= -epsilon*eye(n+m);

Constraints = [c1, c2, c3, c4, c5];

% solve LMI
Objective = [];
options = sdpsettings('solver', 'mosek', 'verbose', 0);
result = optimize(Constraints, Objective, options);

% Update alpha upper/lower bound plus new test value
if result.problem == 0 % if LMIs are feasible
    alpha_low = alpha;
else 
    alpha_up = alpha; % if LMIs are infeasible
end
  
alpha = (alpha_up + alpha_low)/2;
  
end

%% Output

% Get variable indices
vars = [P(:); L(:); V(:); Q11(:)];
dec  = length(unique(getvariables(vars)));

% solutions
data.P    = value(P);
data.L    = value(L);
data.V    = value(V);
data.Q11  = value(Q11);

end

