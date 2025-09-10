% ECS
% University of Southampton
% UK
%
% Date: 5/09/25
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the DT Park criterion
% as defined by Lemma 2 & Theorem 1 (Park, 2019).
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

function [alpha_low,data,dec]=DT_Park(syst,eps,alpha_up,alpha_low)

%% Parameters

A     =  syst.a;
B     = -syst.b; % adopting negative B matrix convention as in Park 2019
C     =  syst.c;
D     =  syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

% Raise error if D != 0.
if ~isequal(D, zeros(m))
    error('D terms must be added to LMIs.');
end

%% Initialising alpha
alpha = alpha_up;

%% Determine max alpha by repeatedly solving LMI for largest feasible alpha

while ((alpha_up-alpha_low)/alpha_up) > eps

% define variables
P11 = sdpvar(n,n,'symmetric');
P12 = sdpvar(n,n,'full');
P13 = sdpvar(n,m,'full');
P14 = sdpvar(n,m,'full');
P22 = sdpvar(n,n,'symmetric');
P23 = sdpvar(n,m,'full');
P24 = sdpvar(n,m,'full');
P33 = sdpvar(m,m,'symmetric');
P34 = sdpvar(m,m,'full');
P44 = sdpvar(m,m,'symmetric');

M1 = sdpvar(m,m,'diagonal');
M2 = sdpvar(m,m,'diagonal');

N1 = sdpvar(m,m,'diagonal');
N2 = sdpvar(m,m,'diagonal');
N3 = sdpvar(m,m,'diagonal');
N4 = sdpvar(m,m,'diagonal');

Q1 = sdpvar(m,m,'diagonal');
Q2 = sdpvar(m,m,'diagonal');
Q3 = sdpvar(m,m,'diagonal');

L1 = sdpvar(m,m,'diagonal');
L2 = sdpvar(m,m,'diagonal');
L3 = sdpvar(m,m,'diagonal');

G1 = sdpvar(n,n,'full');
G2 = sdpvar(n,n,'full');
G3 = sdpvar(m,n,'full');
G4 = sdpvar(m,n,'full');
G5 = sdpvar(m,n,'full');

% define constraints
epsilon = 1e-3;
In  = eye(n);
Im  = eye(m);
c1  = M1 >= epsilon*Im; % PD
c2  = M2 >= epsilon*Im; % PD
c3  = N1 >= epsilon*Im; % PD
c4  = N2 >= epsilon*Im; % PD
c5  = N3 >= epsilon*Im; % PD
c6  = N4 >= epsilon*Im; % PD
c7  = Q1 >= epsilon*Im; % PD
c8  = Q2 >= epsilon*Im; % PD
c9  = Q3 >= epsilon*Im; % PD
c10 = L1 >= epsilon*Im; % PD
c11 = L2 >= epsilon*Im; % PD
c12 = L3 >= epsilon*Im; % PD

% Constraint: preventing G values getting too large.
c13 = [-epsilon*In, G1',         G2',         G3',         G4',         G5';
       G1,          -epsilon*In, zeros(n),    zeros(n,m),  zeros(n,m),  zeros(n,m);
       G2,          zeros(n),    -epsilon*In, zeros(n,m),  zeros(n,m),  zeros(n,m);
       G3,          zeros(m,n),  zeros(m,n),  -epsilon*Im, zeros(m),    zeros(m);
       G4,          zeros(m,n),  zeros(m,n),  zeros(m),    -epsilon*Im, zeros(m);
       G5,          zeros(m,n),  zeros(m,n),  zeros(m),    zeros(m),    -epsilon*Im] <= 0;

% Constraint: Negative def. Lyapunov difference condition from Th. 1 where
% Gamma and Psi matrices are the identity for repeated ReLU.
X11 = -C'*(M2 + N2)*C - P11 + G1*A + A'*G1';
X12 = C'*M2*C - G1 + A'*G2' - P12;
X13 = C'*(L1 + Q1 + Q2) - alpha*G1*B + A'*G3' - P13 + C'*(M2 + N2);
X14 = A'*G4' - C'*(Q1 + M2 + N1) - P14;
X15 = A'*G5' - C'*Q2;
X22 = A'*C'*(M2 + N4)*C*A + C'*(N2 - N4)*C + P11 - P22 + A'*P22*A - G2 - G2' + P12*A + A'*P12' - A'*C'*M2*C - C'*M2'*C*A;
X23 = -C'*(N2 + M2) - alpha*G2*B - G3' - A'*C'*Q2 - C'*Q1 - P23;
X24 = C'*(M2 + N1 + L2 + M1 + N4 + Q1 + Q3) - alpha*A'*C'*(M2 + N4)*C*B - G4' + alpha*C'*M2*C*B + P13 - P24 - ...
      A'*C'*(M1 + N4 + Q3) - alpha*A'*P22*B + A'*P23 - alpha*P12*B;
X25 = A'*C'*(M1 + N3 + Q3 + L3 + Q2) - C'*(M1 + N3 + Q3) - G5' + P14 + A'*P24;
X33 = -2*(L1 + Q1 + Q2) - M1 - M2 - N1 - N2 - P33 - alpha*G3*B - alpha*B'*G3';
X34 = M1 + M2 + N1 + N2 + 2*Q1 + alpha*Q2*C*B - alpha*B'*G4' - P34;
X35 = 2*Q2 - alpha*B'*G5';
X44 = alpha*alpha*B'*C'*(M2 + N4)*C*B + P33 - P44 - 2*(M1 + M2 + L2 + Q1 + Q3) - N1 - N2 - N3 - N4 + ...
      alpha*alpha*B'*P22*B - alpha*B'*P23 - alpha*P23'*B + alpha*B'*C'*(M1 + N4 + Q3) + ...
      alpha*(M1 + N4 + Q3)'*C*B;
X45 = -alpha*B'*C'*(M1 + N3 + L3 + Q2 + Q3) + M1 + M2 + N3 + N4 + 2*Q3 + P34 - alpha*B'*P24;
X55 = -2*(Q2 + Q3 + L3) + P44 - M1 - M2 - N3 - N4; 

c14 = [X11,  X12,  X13,  X14,  X15;
       X12', X22,  X23,  X24,  X25;
       X13', X23', X33,  X34,  X35;
       X14', X24', X34', X44,  X45;
       X15', X25', X35', X45', X55] <= -epsilon*eye(2*n + 3*m);

% Constraint: Positive def. Lyapunov condition from Lemma 2 and Th. 1 where
% Gamma and Psi matrices are the identity for repeated ReLU.
Y11 = -P11 - C'*(M2 + N2)*C;
Y12 = -P12 + C'*M2*C;
Y13 = -P13 + C'*(M2 + N2);
Y14 = -P14 -C'*M2;
Y22 = -P22 -C'*(M2 + N4)*C;
Y23 = -P23 - C'*M2;
Y24 = -P24 + C'*(M2 + N4);
Y33 = -P33 - M1 - M2 - N1 - N2;
Y34 = -P34 + M1 + M2;
Y44 = -P44 - M1 - M2 - N3 - N4;

c15 = [Y11,  Y12,  Y13,  Y14;
       Y12', Y22,  Y23,  Y24;
       Y13', Y23', Y33,  Y34;
       Y14', Y24', Y34', Y44] <= -epsilon*eye(2*n + 2*m);


Constraints = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15];

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
vars = [P11(:); P12(:); P13(:); P14(:); P22(:); P23(:); P24(:); P33(:); P34(:); P44(:); ...
        M1(:); M2(:); N1(:); N2(:); N3(:); N4(:); Q1(:); Q2(:); Q3(:); L1(:); L2(:); L3(:); ...
        G1(:); G2(:); G3(:); G4(:); G5(:)];

dec = length(unique(getvariables(vars)));

% solutions
data.P11 = value(P11);
data.P12 = value(P12);
data.P13 = value(P13);
data.P14 = value(P14);
data.P22 = value(P22);
data.P23 = value(P23);
data.P24 = value(P24);
data.P33 = value(P33);
data.P34 = value(P34);
data.P44 = value(P44);
data.M1 = value(M1);
data.M2 = value(M2);
data.N1 = value(N1);
data.N2 = value(N2);
data.N3 = value(N3);
data.N4 = value(N4);
data.Q1 = value(Q1);
data.Q2 = value(Q2);
data.Q3 = value(Q3);
data.L1 = value(L1);
data.L2 = value(L2);
data.L3 = value(L3);
data.G1 = value(G1);
data.G2 = value(G2);
data.G3 = value(G3);
data.G4 = value(G4);
data.G5 = value(G5);

end
