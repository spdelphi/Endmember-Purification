function [ spectral , Abundance,R,record ] = purify( Endmember , idx , Image , mask , varargin )%,sel, Abundance , varargin )
% Endmember purification affine convex cone exploration
% ----------------------
% Reference(Please cite): 
% Luo, W.F., Gao, L.R., Hong, D.F., Chanussot, J., 2022. Endmember Purification With Affine Simplicial Cone Model. IEEE Transactions on Geoscience and Remote Sensing 60, 23.
% ----------------------
% Email: luowenfei@m.scnu.edu.cn
% ----------------------

%1. find the abundance >30% as well as other endmembers

    method = entryvalue(varargin,'method');
    Abundance = entryvalue(varargin,'abundance');

    use_pj = entryvalue( varargin,'use_pj',false);%use dats projection
    fullrank = entryvalue( varargin,'fullrank',false);%if use_pj, then full rank means not truncte
    spherize = entryvalue( varargin,'spherize',true);%if use_pj, then set spherize or not
    proj_type = entryvalue( varargin,'proj_type','affine_no_noise');

    [bands,N] = size(Image);
    mask_pre = ones(1,N);

    if isempty(mask)
        mask = ones( 1,N);
    end
    newmask = mask & mask_pre;
    %2. unmixing
    oldImage = domask( Image , newmask , 1 );
    if size(oldImage , 2 ) <= 1
        disp('not enough samples, use old mask...');
        Image = domask( Image,mask,1);
    else
        Image = oldImage;
    end



    if isempty(Image) || size(Image , 2 ) <= 1
        disp('not enough samples, skip...');
        spectral = Endmember(:,idx);
        return;
    end


    newE = Endmember;
    if idx ~= 1
        e0 = Endmember(:,idx);
        newE(:,idx) = Endmember(:,1);
        newE(:,1) = e0;
    end



    %data projection
    if use_pj
        EmNum = getEmNum(Image);
        drn = EmNum;
        if fullrank
            drn = bands+1;
        end
        [Image,S] = dr_s(Image,drn,'spherize',spherize,'proj_type',proj_type);

        newE = redr_s(newE,S);
    else

    end

    if guessrays
        disp('guess rays...');
        newE = accselrays_EEA(Image,newE,varargin{:});
        disp('done!');
    end
    EmNum = size(newE,2);

    %initialize abundance or rearrange the abundance if necessary (not for acc model)
    if ~isempty(Abundance)
        if Abundance == 0
            rndab = unifrnd(0,1,EmNum,size(Image,2));
            rndab = rndab ./ sum(rndab);
            Abundance = rndab;
        else
            Abundance = domask(Abundance,mask,1);
            newA = Abundance(idx,:);
            for j = 1:EmNum
                if sel(j)
                    newA = [ newA ; Abundance(j,:) ];
                end
            end
            Abundance = newA;
        end
    end


    record = [];

    if strcmpi(method,'accsearch')
        [Endmember,Abundance,R] = accsearch(Image,newE(:,1),'w_init',newE(:,2:end),varargin{:});
    elseif strcmpi(method,'accsearch_HALS')
        [Endmember,Abundance,R,record] = accsearch_HALS(Image,newE(:,1),'w_init',newE(:,2:end),varargin{:});
    elseif strcmpi(method,'accsearch_sep')
        [Endmember,Abundance,R] = accsearch_sep(Image,newE(:,1),'w_init',newE(:,2:end),'dark',dark,varargin{:});%'using rb');
    end

    % [Endmember,Abundance,rp,record,prefixname]=nlunmixing(Image,Par,num,al,modelType);
    if use_pj
        spectral=undr_s(Endmember,S);
        R = undr_s(R,S);
        if ~isempty(record)
            a = undr_s(record(2:end,:),S);
            record(2:end,:) = [];
            record = [ record;a];
        end
    else
        spectral = Endmember;
    end


end