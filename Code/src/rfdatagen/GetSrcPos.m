% Sua Bae
% 2016-06-16
function mSrcPos = GetSrcPos(pTxTransducer, nTransEleNum_x, nTransEleNum_y, nSubDivNum_x, nSubDivNum_y, sType, nTransRadius)
    
    mSrcPos = zeros(4, nSubDivNum_x*nSubDivNum_y*nTransEleNum_x*nTransEleNum_y); % (x,y,z,eleidx)x(SrcNum)
    
    mTxApertureInfo = xdc_get(pTxTransducer,'rect');
    
    for eidx_y = 1:nTransEleNum_y            
        for eidx_x = 1:nTransEleNum_x
            for sidx_y = 1:nSubDivNum_y
                for sidx_x = 1:nSubDivNum_x        
                    idx =   (eidx_y-1)*nSubDivNum_x*nSubDivNum_y*nTransEleNum_x ...
                          + (eidx_x-1)*nSubDivNum_x*nSubDivNum_y ...
                          + (sidx_y-1)*nSubDivNum_x ...
                          + sidx_x; % index of each source in the matrix 'mSrcPos'
                    mSrcPos(1,idx) = mTxApertureInfo(8, idx); % x poistion of source
                    mSrcPos(2,idx) = mTxApertureInfo(9, idx); % y poistion of source
                    switch sType
                        case 'Linear'
                            mSrcPos(3,idx) = mTxApertureInfo(10,idx);% z poistion of source
                        case 'Convex'
                            mSrcPos(3,idx) = mTxApertureInfo(10,idx) + nTransRadius;% z poistion of source  (Origin of Field II is at immediately below the convex xdcr)
                        otherwise
                            error('undefined');
                    end
                    
                    mSrcPos(4,idx) = (eidx_y-1)*nTransEleNum_x + eidx_x; % index of elements which the source belongs to
                    mSrcPos(5,idx) = eidx_x;  % x index of elements which the source belongs to
                    mSrcPos(6,idx) = eidx_y;  % y index of elements which the source belongs to
                    mSrcPos(7,idx) = sidx_x;  % sub x index within a element
                    mSrcPos(8,idx) = sidx_y;  % sub y index within a element
                end
            end
        end
    end
end