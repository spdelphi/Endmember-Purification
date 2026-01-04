function newB = nmf_acc_update_B(Z,R,B,varargin)
% Endmember purification affine convex cone exploration: calculate B
% ----------------------
% Reference(Please cite): 
% Luo, W.F., Gao, L.R., Hong, D.F., Chanussot, J., 2022. Endmember Purification With Affine Simplicial Cone Model. IEEE Transactions on Geoscience and Remote Sensing 60, 23.
% ----------------------
% Email: luowenfei@m.scnu.edu.cn
% ----------------------

    constraints = entryvalue(varargin,'constraints','L1/2');
    delta = entryvalue(varargin,'deltab',1e-50);

    B(B<delta) = delta;%project to ball
    [m,N] = size(B);
    
    if strcmpi(constraints,'fast')
        newB = (B.* (R'* Z)) ./ ((R' * R * B ));
    elseif strcmpi(constraints,'L1')
        lamda = entryvalue(varargin,'lamda',1e-5);%);%1.5e-2
        lamda_b = entryvalue(varargin,'lamda_b',1e-5);%1e-5);%0
        AL_ITERS = entryvalue(varargin,'AL_ITERS',100);% default is 1000
        lamda_b = sqrt(lamda_b);
        Z_ = [ Z ; lamda_b*B ];
        R_ = [ R ; lamda_b*eye(m)];
        newB = sunsal(R_,Z_,'lambda', lamda, 'POSITIVITY' , 'yes', 'ADDONE', 'no','AL_ITERS',AL_ITERS);
    end
    
    newB(newB<delta) = delta;%project to ball

end

