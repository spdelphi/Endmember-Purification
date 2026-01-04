function [e0,B,R] = accsearch(Y, e0,varargin)
% Endmember purification affine convex cone exploration
% ----------------------
% Reference(Please cite): 
% Luo, W.F., Gao, L.R., Hong, D.F., Chanussot, J., 2022. Endmember Purification With Affine Simplicial Cone Model. IEEE Transactions on Geoscience and Remote Sensing 60, 23.
% ----------------------
% Email: luowenfei@m.scnu.edu.cn
% ----------------------

% Inputs:

% Y (Observation): L*N L-bands,N-the number of pixels

% w_init(L*m):initialization for other endmembers (except for e0).Use to
% initialize the rays.

% m: the number of endmembers (If w_init is  assigned, m can be detected by
% w_init automaticly. So it can be ignored. Otherwise, it should be assigned)

% Outputs:
% e(L*1): endmember e0;
% R (L*p): rays;
% B (P*N): coefficents;

%e.g. [e0,R,B] = accsearch(Y, e0,'w_init',w);
%or [e0,R,B] = accsearch(Y, e0,'m',m)

%local vars
% m: the number of endmembers



%Convergent parameters

E_init = entryvalue(varargin,'w_init',[]);


if isempty(E_init)
    m = entryvalue(varargin,'m',0);
else
    m = size(E_init,2);
end
if m == 0
    error( 'please set the number of endmembers(m) or the initialization of endmembers');
end



Z = Y - e0;
if isempty(E_init)
    % VCA Initialization
    [E_init, ~] = VCA(Y,'init',e0,'Endmembers', m,'verbose','on');
end

R = E_init-e0;
R = normSpectral(R);
m = size( R,2);

max_iter = entryvalue(varargin,'max_iter',1000);
lamda_r = entryvalue(varargin,'lamda_r',1e-4);
lamda_ro = entryvalue(varargin,'lamda_ro',0);
lamda_b = entryvalue(varargin,'lamda_b',1e-4);

lamda = entryvalue(varargin,'lamda',1e-3);%for sparse (default: 1e-4*sqrt(N/m))

iter = 0;

%initialize
B = sunsal(R,Z,'lambda', lamda , 'POSITIVITY' , 'yes', 'ADDONE', 'no','AL_ITERS',500);

eR_org = [e0,R];
eR_old = eR_org;
one = ones(1,size(B,2));

while iter < max_iter
    
    B_ = [one;B];
    M = (Y*B_'+lamda_r*eR_old + lamda_ro*eR_org)/(B_*B_'+lamda_r*eye(m+1));
    M = M./(sqrt(S.* (S1*(M.*M))) + S2 );%It can be replaced by R = normSpectral(R);
    e0=  M(:,1);
    R=M(:,2:end);
    R = normSpectral(R);
    eR_old = [e0,R];
    Z = Y - e0;

    B = nmf_acc_update_B(Z,R,B,'constraints','L1','lamda',lamda,'lamda_b',lamda_b);
    
    iter = iter + 1;
end
