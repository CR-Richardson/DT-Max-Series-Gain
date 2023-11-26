function [alpha,data,dec]=DT_Popov(syst)

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
% Compute the maximum series gain (alpha) when using the discrete-time
% Popov Criterion as defined by Section 6 (Haddad & Bernstein, 1994).
%
% Inputs:
% syst: Structure containing the system matrices of an example.
%
% Returns:
% alpha: Maximum series gain (float)
% data:  Structure containing solutions of the LMI parametrised by alpha
% dec:   # of decision variables
%
% Parameters:
% eps: a small tolerance for better conditioning of pos. def. matrix vars.
% eps2: while loop tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
A     = syst.a;
B     = syst.b;
C     = syst.c;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

eps  = 10^(-4);
eps2 = 10^(-6);

%% Initialising alpha

if m == 1
   Gm = margin(ss(A,B,-C,-D));
   if Gm > 10000
      Gm = 10000;
   end
else
    Gm = 1000;
end

% Determine initial upper/lower bound
alpha_up  = Gm*0.999;
alpha_low = 0; % alpha=0 is always feasible as system's are stable

%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible

while ((alpha_up - alpha_low)/alpha_up) > eps2

alpha = (alpha_up + alpha_low)/2;

setlmis([]);

P = lmivar(1,[n,1]);
V = lmivar(1,kron([1,0],ones(m,1)));
L = lmivar(1,kron([1,0],ones(m,1)));

% LMI
lmiterm([1,1,1,P],A',A);
lmiterm([1,1,1,P],-1,1);
lmiterm([1,1,1,L],(A-eye(n))'*C', C*(A-eye(n)),'s');

lmiterm([1,1,2,P],A',alpha*B);
lmiterm([1,1,2,L],2*(A-eye(n))'*C',alpha*C*B);
lmiterm([1,1,2,L],(A-eye(n))'*C',1);
lmiterm([1,1,2,V],C',1);

lmiterm([1,2,2,P],alpha*B',alpha*B);
lmiterm([1,2,2,L],alpha*B'*C',alpha*C*B,'s');
lmiterm([1,2,2,L],1,alpha*C*B,'s');
lmiterm([1,2,2,V],-1,1,'s');

% P > 0
lmiterm([2,1,1,P],-1,1);
lmiterm([2,1,1,0],eps*eye(n));

% V > 0
lmiterm([3,1,1,V],-1,1);
lmiterm([3,1,1,0],eps*eye(m));

% L > 0
lmiterm([4,1,1,L],-1,1);
lmiterm([4,1,1,0],eps*eye(m));

LMISYS = getlmis;
[tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);

% Update alpha upper/lower bound plus new test value
 if tmin < 0  % if LMIs are feasible
    alpha_low = alpha;
 else 
    alpha_up = alpha; % if LMIs are infeasible
 end
  
end

%% Return solutions
dec     =  decnbr(LMISYS);
data.P  =  dec2mat(LMISYS,xfeas,P);
data.V  =  dec2mat(LMISYS,xfeas,V);
data.L  =  dec2mat(LMISYS,xfeas,L);

if D ~= zeros(m)
    disp('D not equal to zero. Popov criterion may not be applied!');
    alpha   = nan;
    dec     = nan;
    data.P  = nan;
    data.V  = nan;
    data.L  = nan;
end

end

