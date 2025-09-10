function [alpha,data,dec]=DT_Circle_Like(syst)

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
% Circle-like Criterion as defined by Theorem 12.
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

P   = lmivar(1,[n,1]);
V   = lmivar(2,[m,m]);
Q11 = lmivar(2,[m,m]);

% LMI
lmiterm([1,1,1,P],A',A);
lmiterm([1,1,1,P],-1,1);

lmiterm([1,1,2,P],A',alpha*B);
lmiterm([1,1,2,-V],C',1);

lmiterm([1,2,2,P],alpha*B',alpha*B);
lmiterm([1,2,2,Q11],1,1,'s');
lmiterm([1,2,2,V],-1,eye(m)-alpha*D,'s');

% P > 0
lmiterm([2,1,1,P],-1,1);
lmiterm([2,1,1,0],eps*eye(n));

% V: Z-matrix conditions
JJ = eye(m);
count = 3;     
for i = 1:m
    for j = 1:m
        e1 = JJ(i,:); e2 = JJ(:,j);
        if i ~= j
           lmiterm([count,1,1,V],e1,e2,'s');
           count = count+1;
        end
    end
end

% Q11: Positivity matrix conditions
for i = 1:m
    for j = 1:m
        e1 = JJ(i,:); e2 = JJ(:,j);
        lmiterm([count,1,1,Q11],-e1,e2,'s');
        count = count+1;
     end
end      

LMISYS=getlmis;
[tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);

% Update alpha upper/lower bound plus new test value 
 if tmin < 0 % if LMIs are feasible
    alpha_low = alpha;
 else 
    alpha_up = alpha; % if LMIs are infeasible
 end
  
end

%% Return solutions
dec       =  decnbr(LMISYS);
data.P    =  dec2mat(LMISYS,xfeas,P);
data.V    =  dec2mat(LMISYS,xfeas,V);
data.Q11  =  dec2mat(LMISYS,xfeas,Q11);

end

