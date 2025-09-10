function [alpha,data,dec]=DT_Park(syst)

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
% Compute the maximum series gain (alpha) when using the Discrete-time 
% Park Criterion as defined by Lemma 2 & Theorem 1 (Park, 2019).
%
% Note: Due to the large values observed in the G matrices, another
% LMI was added to limit the max. singular values of these matrices.
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
% eps2: tolerance for limiting the max. singular values of the G matrices
% eps3: while loop tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

A     =  syst.a;
B     = -syst.b; % adopting negative B matrix convention as in Park 2019
C     =  syst.c;
D     =  syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

eps  = 10^(-4);
eps2 = 10^(3);
eps3 = 10^(-4); % used instead of 10^(-6) due to computation time

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

while ((alpha_up-alpha_low)/alpha_up) > eps3

alpha = (alpha_up + alpha_low)/2;

setlmis([]);

% Defining matrix variables
P11 = lmivar(1,[n,1]);
P12 = lmivar(2,[n,n]);
P13 = lmivar(2,[n,m]);
P14 = lmivar(2,[n,m]);
P22 = lmivar(1,[n,1]);
P23 = lmivar(2,[n,m]);
P24 = lmivar(2,[n,m]);
P33 = lmivar(1,[m,1]);
P34 = lmivar(2,[m,m]);
P44 = lmivar(1,[m,1]);

M1  = lmivar(1,kron([1,0],ones(m,1)));
M2  = lmivar(1,kron([1,0],ones(m,1)));

N1  = lmivar(1,kron([1,0],ones(m,1)));
N2  = lmivar(1,kron([1,0],ones(m,1)));
N3  = lmivar(1,kron([1,0],ones(m,1)));
N4  = lmivar(1,kron([1,0],ones(m,1)));

% Pi matrices
Q1  = lmivar(1,kron([1,0],ones(m,1)));
Q2  = lmivar(1,kron([1,0],ones(m,1)));
Q3  = lmivar(1,kron([1,0],ones(m,1)));

% Lambda matrices
L1  = lmivar(1,kron([1,0],ones(m,1)));
L2  = lmivar(1,kron([1,0],ones(m,1)));
L3  = lmivar(1,kron([1,0],ones(m,1)));

G1 = lmivar(2,[n,n]);
G2 = lmivar(2,[n,n]);
G3 = lmivar(2,[m,n]);
G4 = lmivar(2,[m,n]);
G5 = lmivar(2,[m,n]);

% First LMI - Positive def. Lyapunov condition from Lemma 2 and Th. 1
% Note: Gamma and Psi matrices are the identity for repeated ReLU

lmiterm([1,1,1,P11],-1,1);
lmiterm([1,1,1,M2],-C',C);
lmiterm([1,1,1,N2],-C',C);

lmiterm([1,1,2,P12],-1,1);
lmiterm([1,1,2,M2],C',C);

lmiterm([1,1,3,P13],-1,1);
lmiterm([1,1,3,M2],C',1);
lmiterm([1,1,3,N2],C',1);

lmiterm([1,1,4,P14],-1,1);
lmiterm([1,1,4,M2],-C',1);

lmiterm([1,2,2,P22],-1,1);
lmiterm([1,2,2,M2],-C',C);
lmiterm([1,2,2,N4],-C',C);

lmiterm([1,2,3,P23],-1,1);
lmiterm([1,2,3,M2],-C',1);

lmiterm([1,2,4,P24],-1,1);
lmiterm([1,2,4,M2],C',1);
lmiterm([1,2,4,N4],C',1);

lmiterm([1,3,3,P33],-1,1);
lmiterm([1,3,3,M1],-1,1);
lmiterm([1,3,3,M2],-1,1);
lmiterm([1,3,3,N1],-1,1);
lmiterm([1,3,3,N2],-1,1);

lmiterm([1,3,4,P34],-1,1);
lmiterm([1,3,4,M1],1,1);
lmiterm([1,3,4,M2],1,1);

lmiterm([1,4,4,P44],-1,1);
lmiterm([1,4,4,M1],-1,1);
lmiterm([1,4,4,M2],-1,1);
lmiterm([1,4,4,N3],-1,1);
lmiterm([1,4,4,N4],-1,1);

% Second LMI - Negative def. Lyapunov difference condition from Th. 1
% Note: Gamma and Psi matrices are the identity for repeated ReLU

lmiterm([2,1,1,M2],-C',C);
lmiterm([2,1,1,N2],-C',C);
lmiterm([2,1,1,P11],-1,1);
lmiterm([2,1,1,G1],1,A,'s');

lmiterm([2,1,2,M2],C',C);
lmiterm([2,1,2,G1],-1,1);
lmiterm([2,1,2,-G2],A',1);
lmiterm([2,1,2,P12],-1,1);

lmiterm([2,1,3,L1],C',1);
lmiterm([2,1,3,Q1],C',1);
lmiterm([2,1,3,Q2],C',1);
lmiterm([2,1,3,G1],1,-alpha*B);
lmiterm([2,1,3,-G3],A',1);
lmiterm([2,1,3,P13],-1,1);
lmiterm([2,1,3,M2],C',1);
lmiterm([2,1,3,N2],C',1);

lmiterm([2,1,4,-G4],A',1);
lmiterm([2,1,4,Q1],-C',1);
lmiterm([2,1,4,M2],-C',1);
lmiterm([2,1,4,N1],-C',1);
lmiterm([2,1,4,P14],-1,1);

lmiterm([2,1,5,-G5],A',1);
lmiterm([2,1,5,Q2],-C',1);

lmiterm([2,2,2,M2],A'*C',C*A);
lmiterm([2,2,2,N4],A'*C',C*A);
lmiterm([2,2,2,N2],C',C);
lmiterm([2,2,2,P11],1,1);
lmiterm([2,2,2,N4],-C',C);
lmiterm([2,2,2,P22],-1,1);
lmiterm([2,2,2,P22],A',A);
lmiterm([2,2,2,G2],-1,1,'s');
lmiterm([2,2,2,P12],1,A,'s');
lmiterm([2,2,2,M2],-A'*C',C,'s');

lmiterm([2,2,3,M2],-C',1);
lmiterm([2,2,3,N2],-C',1);
lmiterm([2,2,3,G2],-1,alpha*B);
lmiterm([2,2,3,-G3],-1,1);
lmiterm([2,2,3,Q2],-A'*C',1);
lmiterm([2,2,3,Q1],-C',1);
lmiterm([2,2,3,P23],-1,1);

lmiterm([2,2,4,M2],C',1);
lmiterm([2,2,4,N1],C',1);
lmiterm([2,2,4,L2],C',1);
lmiterm([2,2,4,M2],-A'*C',alpha*C*B);
lmiterm([2,2,4,-G4],-1,1);
lmiterm([2,2,4,M2],C',alpha*C*B);
lmiterm([2,2,4,P13],1,1);
lmiterm([2,2,4,P24],-1,1);
lmiterm([2,2,4,M1],-A'*C',1);
lmiterm([2,2,4,N4],-A'*C',1);
lmiterm([2,2,4,Q3],-A'*C',1);
lmiterm([2,2,4,M1],C',1);
lmiterm([2,2,4,N4],C',1);
lmiterm([2,2,4,Q1],C',1);
lmiterm([2,2,4,Q3],C',1);
lmiterm([2,2,4,N4],-A'*C',alpha*C*B);
lmiterm([2,2,4,P22],-A',alpha*B);
lmiterm([2,2,4,P23],A',1);
lmiterm([2,2,4,P12],-1,alpha*B);

lmiterm([2,2,5,M1],A'*C',1);
lmiterm([2,2,5,N3],A'*C',1);
lmiterm([2,2,5,Q3],A'*C',1);
lmiterm([2,2,5,M1],-C',1);
lmiterm([2,2,5,N3],-C',1);
lmiterm([2,2,5,Q3],-C',1);
lmiterm([2,2,5,L3],A'*C',1);
lmiterm([2,2,5,Q2],A'*C',1);
lmiterm([2,2,5,-G5],-1,1);
lmiterm([2,2,5,P14],1,1);
lmiterm([2,2,5,P24],A',1);

lmiterm([2,3,3,L1],-2,1);
lmiterm([2,3,3,M1],-1,1);
lmiterm([2,3,3,M2],-1,1);
lmiterm([2,3,3,N1],-1,1);
lmiterm([2,3,3,N2],-1,1);
lmiterm([2,3,3,Q1],-2,1);
lmiterm([2,3,3,Q2],-2,1);
lmiterm([2,3,3,P33],-1,1);
lmiterm([2,3,3,G3],-1,alpha*B,'s');

lmiterm([2,3,4,M1],1,1);
lmiterm([2,3,4,M2],1,1);
lmiterm([2,3,4,N1],1,1);
lmiterm([2,3,4,N2],1,1);
lmiterm([2,3,4,Q1],2,1);
lmiterm([2,3,4,Q2],1,alpha*C*B);
lmiterm([2,3,4,-G4],-alpha*B',1);
lmiterm([2,3,4,P34],-1,1);

lmiterm([2,3,5,Q2],2,1);
lmiterm([2,3,5,-G5],-alpha*B',1);

lmiterm([2,4,4,M2],alpha*B'*C',alpha*C*B);
lmiterm([2,4,4,P33],1,1);
lmiterm([2,4,4,N4],alpha*B'*C',alpha*C*B);
lmiterm([2,4,4,P44],-1,1);
lmiterm([2,4,4,M1],-2,1);
lmiterm([2,4,4,M2],-2,1);
lmiterm([2,4,4,N1],-1,1);
lmiterm([2,4,4,N2],-1,1);
lmiterm([2,4,4,N3],-1,1);
lmiterm([2,4,4,N4],-1,1);
lmiterm([2,4,4,L2],-2,1);
lmiterm([2,4,4,Q1],-2,1);
lmiterm([2,4,4,Q3],-2,1);
lmiterm([2,4,4,P22],alpha*B',alpha*B);
lmiterm([2,4,4,P23],-alpha*B',1,'s');
lmiterm([2,4,4,M1],alpha*B'*C',1,'s');
lmiterm([2,4,4,N4],alpha*B'*C',1,'s');
lmiterm([2,4,4,Q3],alpha*B'*C',1,'s');

lmiterm([2,4,5,M1],-alpha*B'*C',1);
lmiterm([2,4,5,N3],-alpha*B'*C',1);
lmiterm([2,4,5,L3],-alpha*B'*C',1);
lmiterm([2,4,5,Q2],-alpha*B'*C',1);
lmiterm([2,4,5,Q3],-alpha*B'*C',1);
lmiterm([2,4,5,M1],1,1);
lmiterm([2,4,5,M2],1,1);
lmiterm([2,4,5,N3],1,1);
lmiterm([2,4,5,N4],1,1);
lmiterm([2,4,5,Q3],2,1);
lmiterm([2,4,5,P34],1,1);
lmiterm([2,4,5,P24],-alpha*B',1);

lmiterm([2,5,5,Q2],-2,1);
lmiterm([2,5,5,Q3],-2,1);
lmiterm([2,5,5,P44],1,1);
lmiterm([2,5,5,L3],-2,1);
lmiterm([2,5,5,M1],-1,1);
lmiterm([2,5,5,M2],-1,1);
lmiterm([2,5,5,N3],-1,1);
lmiterm([2,5,5,N4],-1,1);

% LMIs constraining individual matrix variables
lmiterm([3,1,1,M1],-1,1);  % M1 > 0
lmiterm([3,1,1,0],eps*eye(m));
lmiterm([4,1,1,M2],-1,1);  % M2 > 0
lmiterm([4,1,1,0],eps*eye(m));
lmiterm([5,1,1,N1],-1,1);  % N1 > 0
lmiterm([5,1,1,0],eps*eye(m));
lmiterm([6,1,1,N2],-1,1);  % N2 > 0
lmiterm([6,1,1,0],eps*eye(m));
lmiterm([7,1,1,N3],-1,1);  % N3 > 0
lmiterm([7,1,1,0],eps*eye(m));
lmiterm([8,1,1,N4],-1,1);  % N4 > 0
lmiterm([8,1,1,0],eps*eye(m));
lmiterm([9,1,1,Q1],-1,1);  % Q1 > 0
lmiterm([9,1,1,0],eps*eye(m));
lmiterm([10,1,1,Q2],-1,1); % Q2 > 0
lmiterm([10,1,1,0],eps*eye(m));
lmiterm([11,1,1,Q3],-1,1); % Q3 > 0
lmiterm([11,1,1,0],eps*eye(m));
lmiterm([12,1,1,L1],-1,1); % L1 > 0
lmiterm([12,1,1,0],eps*eye(m));
lmiterm([13,1,1,L2],-1,1); % L2 > 0
lmiterm([13,1,1,0],eps*eye(m));
lmiterm([14,1,1,L3],-1,1); % L3 > 0
lmiterm([14,1,1,0],eps*eye(m));

% Prevent G values getting so large
lmiterm([15,1,1,0],-eps2*eye(n));
lmiterm([15,1,2,-G1],1,1);
lmiterm([15,1,3,-G2],1,1);
lmiterm([15,1,4,-G3],1,1);
lmiterm([15,1,5,-G4],1,1);
lmiterm([15,1,6,-G5],1,1);
lmiterm([15,2,2,0],-eps2*eye(n));
lmiterm([15,3,3,0],-eps2*eye(n));
lmiterm([15,4,4,0],-eps2*eye(m));
lmiterm([15,5,5,0],-eps2*eye(m));
lmiterm([15,6,6,0],-eps2*eye(m));

LMISYS       = getlmis;
[tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);

% Update alpha upper/lower bound plus new test value 
 if tmin < 0 % if LMIs are feasible
    alpha_low = alpha;
  else 
    alpha_up = alpha; % if LMIs are infeasible
 end

end

%% Return solution
dec       =  decnbr(LMISYS);
data.P11  =  dec2mat(LMISYS,xfeas,P11);
data.P12  =  dec2mat(LMISYS,xfeas,P12);
data.P13  =  dec2mat(LMISYS,xfeas,P13);
data.P14  =  dec2mat(LMISYS,xfeas,P14);
data.P22  =  dec2mat(LMISYS,xfeas,P22);
data.P23  =  dec2mat(LMISYS,xfeas,P23);
data.P24  =  dec2mat(LMISYS,xfeas,P24);
data.P33  =  dec2mat(LMISYS,xfeas,P33);
data.P34  =  dec2mat(LMISYS,xfeas,P34);
data.P44  =  dec2mat(LMISYS,xfeas,P44);
data.M1   =  dec2mat(LMISYS,xfeas,M1);
data.M2   =  dec2mat(LMISYS,xfeas,M2);
data.N1   =  dec2mat(LMISYS,xfeas,N1);
data.N2   =  dec2mat(LMISYS,xfeas,N2);
data.N3   =  dec2mat(LMISYS,xfeas,N3);
data.N4   =  dec2mat(LMISYS,xfeas,N4);
data.Q1   =  dec2mat(LMISYS,xfeas,Q1);
data.Q2   =  dec2mat(LMISYS,xfeas,Q2);
data.Q3   =  dec2mat(LMISYS,xfeas,Q3);
data.L1   =  dec2mat(LMISYS,xfeas,L1);
data.L2   =  dec2mat(LMISYS,xfeas,L2);
data.L3   =  dec2mat(LMISYS,xfeas,L3);
data.G1   =  dec2mat(LMISYS,xfeas,G1);
data.G2   =  dec2mat(LMISYS,xfeas,G2);
data.G3   =  dec2mat(LMISYS,xfeas,G3);
data.G4   =  dec2mat(LMISYS,xfeas,G4);
data.G5   =  dec2mat(LMISYS,xfeas,G5);

if D ~= zeros(m)
    disp('D not equal to zero. Park criterion may not be applied!');
    alpha     =  nan;
    dec       =  nan;
    data.P11  =  nan;
    data.P12  =  nan;
    data.P13  =  nan;
    data.P14  =  nan;
    data.P22  =  nan;
    data.P23  =  nan;
    data.P24  =  nan;
    data.P33  =  nan;
    data.P34  =  nan;
    data.P44  =  nan;
    data.M1   =  nan;
    data.M2   =  nan;
    data.N1   =  nan;
    data.N2   =  nan;
    data.N3   =  nan;
    data.N4   =  nan;
    data.Q1   =  nan;
    data.Q2   =  nan;
    data.Q3   =  nan;
    data.L1   =  nan;
    data.L2   =  nan;
    data.L3   =  nan;
    data.G1   =  nan;
    data.G2   =  nan;
    data.G3   =  nan;
    data.G4   =  nan;
    data.G5   =  nan;
end

end
