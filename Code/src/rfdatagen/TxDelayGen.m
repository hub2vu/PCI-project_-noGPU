% 2017-08-26
%  updatate: PW
%
% 2017-06-12
%   update: CF/linear
%   nScIdx: 1 ~ nNumScline
% 2017-03-04
%
% Sua Bae

function  [aTxDelay, nDelayOffset] = TxDelayGen(stTrans, stTx, mSrcPos, nSrcNum, nTxAngle, nNumScline, nScIdx, aTxAperture, nSoundSpeed, bStartAcqTime)

    if any(aTxAperture) 

        switch stTx.sType
            
            case 'PW'  %%%%%%% 1. Plane wave transmission %%%%%%%%%%%%%%%
                
                %%%%%%% 1) FOR LINEAR ARRAY & 2) CONVEX ARRAY %%%%%%%
                aTxDelayDist_src = CalcDelayDist(stTrans.sType, mSrcPos, 'PW', nTxAngle, 0, 0); % [meter]
                aTxDelayDist_ele = aTxDelayDist_src(1:nSrcNum:end);
                aTxDelay_tmp =  aTxDelayDist_ele/nSoundSpeed; % full aperture delay [sec] <- [meter]
                switch bStartAcqTime                    
                    case 'first firing time'
                        nDelayOffset = min(aTxDelay_tmp(logical(aTxAperture))); % minimum delay of elements which are excited (included in the Tx aperture).
                        aTxDelay = aTxDelay_tmp - nDelayOffset; % [sec] % Delays are adjusted so that minimum delay can be zero
                                                                           % and t = 0 is the excitation time of the first excited element.
                        aTxDelay(~logical(aTxAperture)) = 0; % [sec] % Set unused delay values to zero.
                    case 'plane wave placed at origin'
                        aTxDelay =  aTxDelay_tmp; % [sec] <- [meter]
                    otherwise
                        error('undefined');
                end
                
                
            case 'CF' %%%%%%% 2. Focused beam transmission %%%%%%%%%%%%%%%
                
                switch stTrans.sType
                    case 'Linear' %%%%%%% 1) FOR LINEAR ARRAY

                        nSclinePitch = stTrans.nNumEle_x*stTrans.nPitch_x/nNumScline;
                        % aFocalPos_x  = nSclinePitch*(nScIdx-(nNumScline+1)/2);
                        % 
                        % aFocalPos = [aFocalPos_x, 0, stTx.nFocalDepth]; % (x,y,z) [m]
                        % 
                        % aTxDelayDist_src = CalcDelayDist(stTrans.sType, mSrcPos, 'CF', 0, aFocalPos, 0); % [meter]
                        % aTxDelayDist_ele = aTxDelayDist_src(1:nSrcNum:end);
                        % aTxDelay_tmp =  aTxDelayDist_ele/nSoundSpeed; % full aperture delay [sec] <- [meter]
                        % 
                        % nMinDelay = min(aTxDelay_tmp(logical(aTxAperture))); % minimum delay of elements which are excited (included in the Tx aperture).
                        % aTxDelay = aTxDelay_tmp - nMinDelay; % [sec] % Delays are adjusted so that minimum delay can be zero
                        % aTxDelay(~logical(aTxAperture)) = 0; % [sec] % Set unused delay values to zero.
                        % nDelayOffset = max(aTxDelay); % [sec] % time between 'first firing time' and 'when all transmissions completed'

                        a = sqrt( ( ((1:(nNumScline))-nScIdx)*nSclinePitch ).^2 + stTx.nFocalDepth^2)/nSoundSpeed;
                        b = stTx.nFocalDepth/nSoundSpeed;
                        c = b-a;
                        m = min(c(logical(aTxAperture)));
                        d = c-m;
                        d(~logical(aTxAperture)) = 0;

                        aTxDelay = d;                        
                        nDelayOffset = max(aTxDelay);

                    case 'Convex' %%%%%%% 2) FOR CONVEX ARRAY   
                        aFocalPos = (stTrans.nRadius + stTx.nFocalDepth)*[sind(nTxAngle), 0, cosd(nTxAngle)]; % (x,y,z) [m]

                        aTxDelayDist_src = CalcDelayDist(stTrans.sType, mSrcPos, 'CF', 0, aFocalPos, 0); % [meter]
                        aTxDelayDist_ele = aTxDelayDist_src(1:nSrcNum:end);
                        aTxDelay_tmp =  aTxDelayDist_ele/nSoundSpeed; % full aperture delay [sec] <- [meter]

                        nMinDelay = min(aTxDelay_tmp(logical(aTxAperture))); % minimum delay of elements which are excited (included in the Tx aperture).
                        aTxDelay = aTxDelay_tmp - nMinDelay; % [sec] % Delays are adjusted so that minimum delay is zero
                        aTxDelay(~logical(aTxAperture)) = 0; % [sec] % Set unused delay values to zero.
                        nDelayOffset = max(aTxDelay); % [sec] % time between 'first firing time' and 'when all transmissions completed'                        

                        nSclinePitchAngle = stTrans.nNumEle_x*stTrans.nPitchAngle_x/nNumScline;
                        aa = ((1:(nNumScline))-nScIdx)*nSclinePitchAngle;
                        a = sqrt ( (stTrans.nRadius*sind(aa)).^2 + (stTrans.nRadius + stTx.nFocalDepth - stTrans.nRadius*cosd(aa)).^2 )/nSoundSpeed;
                        b = stTx.nFocalDepth/nSoundSpeed;
                        c = b-a;
                        m = min(c(logical(aTxAperture)));
                        d = c-m;
                        d(~logical(aTxAperture)) = 0;

                        aTxDelay = d;                        
                        nDelayOffset = max(aTxDelay);

                    otherwise
                        error('undefined');
                end

        end
        
    else
        % No excited element
        aTxDelay = zeros(1,128);
    end
            
end