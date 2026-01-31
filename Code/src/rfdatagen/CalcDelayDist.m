%
% 2016-09-30
%  modify variable names
% 
% 2015-09-10
%
% mSrcPos       : source poisition and element index, dim = (xpos, ypos, zpos, eleidx)x(nRxAptrEleNum*nSrcNumPerTransEle)
% sTXType       : 'PW' - plane wave
%                 'CF' - conventional focusing
% nPWAngle      : Steering angle of plane wave, [theta, psi] when 'Matrix' [degree]
% aFocalPointPos    : Position of focal point [meter], dim: 3x1 (1st row: x position, 2nd row: y position, 3rd row: z position)
%
% bIndependentSrc = 0; % 1: all sources have their own delays, 0: sources belonging to the same element have the same delays
%
% aSrcDelayDist : Delay distance at each source point [meter] (dim: 1xSrcNum)
%
function [aSrcDelayDist] = CalcDelayDist(sTransType, mSrcPos, sTXType, nPWAngle, aFocalPointPos, bIndependentSrc)

    
    
    switch sTXType
        
        case 'PW'         
            
            if(~bIndependentSrc)
                
                % Calculate Element Position
                nEleNum = max(mSrcPos(4,:));
                for eidx = 1:nEleNum
                    aEle_x(eidx) = mean(mSrcPos(1,(mSrcPos(4,:)==eidx))); % x-pos of the eidx-th element (= mean x-pos of elements belong to the eidx-th element)
                    aEle_y(eidx) = mean(mSrcPos(2,(mSrcPos(4,:)==eidx))); % y-pos of the eidx-th element (= mean y-pos of elements belong to the eidx-th element)
                    aEle_z(eidx) = mean(mSrcPos(3,(mSrcPos(4,:)==eidx))); % z-pos of the eidx-th element (= mean z-pos of elements belong to the eidx-th element)
                end
                
                % Calculate Delay Distance at Each Element
                switch sTransType
                    case 'Linear'
                        aEleToPW = aEle_x.*sind(nPWAngle); % Distance from the PW at origin to target [meter] % aEleToPW: dim = 1x(nEleNum)        
                    case 'Convex'
                        %%%%%%%%%%% Ignoring aEle_y 
                        aEle_r = sqrt( aEle_x.^2 + aEle_z.^2);  % r position of element [meter] (r,theta)
                        aEle_theta = atand(aEle_x./aEle_z);      % theta position of element [meter] (r,theta)
                        aEleToPW = aEle_r.*cosd(aEle_theta-nPWAngle); % Distance from the PW at origin to target [meter] % aEleToPW: dim = 1x(nEleNum) 
                    case 'Matrix'
                        nPWAngle_az = nPWAngle(1);
                        nPWAngle_el = nPWAngle(2);
                        aEleToPW = aEle_x*sind(nPWAngle_az) + aEle_y*cosd(nPWAngle_az)*sind(nPWAngle_el) + aEle_z*cosd(nPWAngle_az)*cosd(nPWAngle_el);
                end
                   
                % Set Delay Distance of Each Source
                nSrcNum = size(mSrcPos,2);
                for sidx = 1:nSrcNum
                    aSrcDelayDist(sidx) = aEleToPW(mSrcPos(4,sidx)); % aSrcDelayDist: dim = 1x(nSrcNum)
                end                    
                    
            else
                switch sTransType
                    case 'Linear'
                        %%%%%%%%%%%% Ignoring y position of sources 
                        aSrc_r = sqrt( mSrcPos(1,:).^2 + mSrcPos(3,:).^2 );  % r position of source [meter] (r,theta)
                        aSrc_theta = atand(mSrcPos(1,:)./mSrcPos(3,:)); % theta position of source [meter] (r,theta)
                        aSrc_theta(isnan(aSrc_theta)) = 90;

                        aSrcToPW = aSrc_r.*cosd(aSrc_theta-nPWAngle); % Distance from the PW at origin to target [meter] 
                        aSrcDelayDist = aSrcToPW;
                    case 'Convex'
                        %%%%%%%%%%%% Ignoring y position of sources 
                        aSrc_r = sqrt( mSrcPos(1,:).^2 + mSrcPos(3,:).^2 );  % r position of source [meter] (r,theta)
                        aSrc_theta = atand(mSrcPos(1,:)./mSrcPos(3,:)); % theta position of source [meter] (r,theta)
                        aSrc_theta(isnan(aSrc_theta)) = 90;

                        aSrcToPW = aSrc_r.*cosd(aSrc_theta-nPWAngle); % Distance from the PW at origin to target [meter] 
                        aSrcDelayDist = aSrcToPW;
                    case 'Matrix'
                        nPWAngle_az = nPWAngle(1);
                        nPWAngle_el = nPWAngle(2);
                        aSrcToPW = mSrcPos(1,:)*sind(nPWAngle_az) + mSrcPos(2,:)*cosd(nPWAngle_az)*sind(nPWAngle_el) + mSrcPos(3,:)*cosd(nPWAngle_az)*cosd(nPWAngle_el);
                        aSrcDelayDist = aSrcToPW;
                end
            end
            
        case 'CF'  
                        
            if(~bIndependentSrc)
                nEleNum = max(mSrcPos(4,:)); % 4th row : element index
                for eidx = 1:nEleNum
                    aEle_x(eidx) = mean(mSrcPos(1,(mSrcPos(4,:)==eidx))); % x-pos of the eidx-th element (= mean x-pos of elements belong to the eidx-th element)
                    aEle_y(eidx) = mean(mSrcPos(2,(mSrcPos(4,:)==eidx))); % z-pos of the eidx-th element (= mean y-pos of elements belong to the eidx-th element)
                    aEle_z(eidx) = mean(mSrcPos(3,(mSrcPos(4,:)==eidx))); % z-pos of the eidx-th element (= mean z-pos of elements belong to the eidx-th element)
                end

                aEleToFocalPoint = sqrt( (aEle_x-aFocalPointPos(1)).^2 + (aEle_y-aFocalPointPos(2)).^2 + (aEle_z-aFocalPointPos(3)).^2 ); % Distance from the PW at origin to target [meter]                 
                aEleDelayDist = max(aEleToFocalPoint) - aEleToFocalPoint; % aEleDelayDist: dim = 1x(nEleNum)
                
                nSrcNum = size(mSrcPos,2);
                for sidx = 1:nSrcNum
                    aSrcDelayDist(sidx) = aEleDelayDist(mSrcPos(4,sidx)); % aSrcDelayDist: dim = 1x(nSrcNum)
                end                    
                    
            else
                aSrcToFocalPoint = sqrt( (mSrcPos(1,:)-aFocalPointPos(1)).^2 + (mSrcPos(2,:)-aFocalPointPos(2)).^2 + (mSrcPos(3,:)-aFocalPointPos(3)).^2 ); 
                aSrcDelayDist = max(aSrcToFocalPoint) - aSrcToFocalPoint; 
            end
    end
    
end