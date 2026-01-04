function [ newE ] = purifyall( Endmember , Image , mask , varargin )
% This code is written by Wenfei Luo. Please cite the following reference if
% you use this code.
% ----------------------
% Reference(Please cite): 
% Luo, W.F., Gao, L.R., Hong, D.F., Chanussot, J., 2022. Endmember Purification With Affine Simplicial Cone Model. IEEE Transactions on Geoscience and Remote Sensing 60, 23.

%mask can be made by the user

    newE = [];
    EmNum = size( Endmember , 2 );

    if isempty(mask)
        mask = ones(EmNum,size(Image,2));
    end
    for i = 1:EmNum
        disp(['processing endmember:',num2str(i)]);
        [e,a] = purify( Endmember , i , Image , mask(i,:) , varargin{:} );
        newE = [ newE , e(:,1) ];
    end

end