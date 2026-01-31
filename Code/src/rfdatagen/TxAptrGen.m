% 2019-04-01
%   update 'CF/Linear'
%   stTrans  must have   mElePos (x,y,z,theta)
%   stTx     must have   mFocalPos (x,y,z,theta)  when stTx.sType='CF'
%
% 2019-01-27
%   Apod -> Aperture
%
% 2017-03-04
%
% Sua Bae

function [aTxAptr] = TxAptrGen(stTrans, stTx, txidx)

    aTxAptr = zeros(1, stTrans.nNumEle); % dim = (pRFInfo.stTrans.nEleNum)x(1) 
    
    switch stTx.sType
        
        case 'CF'
            switch stTrans.sType
                case 'Linear'
                    nSclinePitch = stTrans.nNumEle_x*stTrans.nPitch_x/stTx.nNum;
                    lft_x = nSclinePitch*(txidx-(stTx.nNum+1)/2) - (stTx.mFocalPos(txidx,3)/stTx.nFnum)/2;
                    rgt_x = nSclinePitch*(txidx-(stTx.nNum+1)/2) + (stTx.mFocalPos(txidx,3)/stTx.nFnum)/2;                 
                    aTxAptr((lft_x < stTrans.mElePos(:,1)) &(stTrans.mElePos(:,1) < rgt_x)) = 1;
                case 'Convex'
                    a = 4*stTx.nFnum*stTx.nFnum + 1;
                    b = -2*stTx.nFnum*(stTx.mFocalPos(txidx,4));
                    c = stTx.mFocalPos(txidx,4)^2 - stTrans.nRadius^2;
                    nAptLen = 2*(-b -sqrt(b*b-a*c))/a; % lateral width of arc aperture
                    lft = round(( stTx.mFocalPos(txidx,5) - asind(nAptLen/2/stTrans.nRadius) )/stTrans.nPitchAngle_x + (stTrans.nNumEle_x+1)/2);
                    rt  = round(( stTx.mFocalPos(txidx,5) + asind(nAptLen/2/stTrans.nRadius) )/stTrans.nPitchAngle_x + (stTrans.nNumEle_x+1)/2);          
                    if lft < 1, lft = 1; end;
                    if rt > stTrans.nNumEle, rt = stTrans.nNumEle; end;
                    aTxAptr(lft:rt) = 1; % activate!
                otherwise
                    error('undefined');
            end
            
        case 'PW'
            switch stTrans.sType
                case 'Linear'
                    aTxAptr(:) = 1; % all elements activated!
                case 'Convex'
                    nTxAngle = stTx.aPwAngle_deg(txidx);
                    aTxAptr((abs(stTrans.mElePos(:,4) - nTxAngle) < stTx.nAptrAngle)) = 1;
%                     nLeftMostEle_theta = -stTrans.nMaxTheta; % theta position of the leftmost element in array[degree] 
%                     nTxAptrNumElHalf = round((stTrans.nRadius*stTx.nAptrAngle/180*pi)/stTrans.nPitch_x); % no. of elements in 1/2 aperture
%                     nTxAptrCenterEle = round(stTrans.nRadius*(nTxAngle-nLeftMostEle_theta)/180*pi/stTrans.nPitch_x+1);
%                     lft = max(nTxAptrCenterEle - nTxAptrNumElHalf, 1);
%                     rt  = min(nTxAptrCenterEle + nTxAptrNumElHalf, stTrans.nNumEle);
%                     aTxApod(lft:rt) = 1; % activate!
                    
                otherwise
                    error('undefined');
            end
            
        case 'DW'
            switch stTrans.sType
                case 'Linear'
                    error('undefined');
                case 'Convex'
                    aAng_VS2Ele = atand((stTrans.mElePos(:,1)-stTx.mVSPos(txidx,1))./(stTrans.mElePos(:,3)-stTx.mVSPos(txidx,3)));                
                    aTxAptr((abs(stTrans.mElePos(:,4) - aAng_VS2Ele) < stTx.nAptrAngle)) = 1;
                    
                otherwise
                    error('undefined');
            end
        otherwise
            error('undefined');
    end
    
end