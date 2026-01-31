% Sua Bae
% 2019-04-01
%
classdef clDir
   
    properties (SetAccess = public)
        % common attributes
        exp_idx         % experiment index
        phantom_idx     % phantom index
        acqType_idx     % acquisition type index
        ang_deg         % data number index 
        
        sExp            % experiment type
        sPhantom        % phantom name
        sAcqType        % acquisition type
        
        sExpDir         % experiment directory
        sPhtmDir        % phantom directory
        sAcqDir         % acquisition directory
        sAngDir         % angle directory (only for 'PW')
        sDataDir        % data directory
    end
    
    %-- Constructor
    methods (Access = public)
        
        function h = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, sTopDir)   
            % ang_deg will be ignored if acqType_idx = 1 (CF)
            if acqType_idx == 1
                h.exp_idx           = exp_idx;
                h.acqType_idx       = acqType_idx;
            elseif  acqType_idx == 2||acqType_idx == 3
                h.exp_idx           = exp_idx;
                h.acqType_idx       = acqType_idx;
                h.ang_deg           = ang_deg;
            else 
                error('undefined');
            end
            
            %%--Experiment Type: Simulation, Phantom, In vivo
            if     exp_idx == 1; h.sExp = 'Simulation';
            elseif exp_idx == 2; h.sExp = 'Phantom';
            elseif exp_idx == 3; h.sExp = 'In vivo';
            else   error('undefined');
            end
            h.sExpDir = [sTopDir '\' num2str(exp_idx) '. ' h.sExp '\'];
            
            %%--Phantom idx
            if exp_idx == 1 %  'Simulation'
                if     phantom_idx == 1; h.sPhantom = 'point 10 mm';
                elseif phantom_idx == 2; h.sPhantom = 'point 20 mm';
                elseif phantom_idx == 3; h.sPhantom = 'point 30 mm';
                elseif phantom_idx == 4; h.sPhantom = 'point 40 mm';
                elseif phantom_idx == 5; h.sPhantom = 'point 50 mm';
                elseif phantom_idx == 6; h.sPhantom = 'point 60 mm';
                elseif phantom_idx == 7; h.sPhantom = '6 point targets';
                else   error('undefined');
                end
            elseif exp_idx == 2 % 'Phantom'
                if     phantom_idx == 1; h.sPhantom = 'ATS';
                else   error('undefined');
                end
            elseif exp_idx == 3 % 'In vivo'
                if     phantom_idx == 1; h.sPhantom = 'Carotid 1';
                elseif     phantom_idx == 2; h.sPhantom = 'Carotid 2';
                elseif     phantom_idx == 3; h.sPhantom = 'Carotid 3';
                elseif     phantom_idx == 4; h.sPhantom = 'Carotid 4';
                else   error('undefined');
                end
            else   error('undefined');
            end
            h.sPhtmDir = [h.sExpDir num2str(phantom_idx) '. ' h.sPhantom '\'];
            
            %%--Acquisition Type
            if     acqType_idx == 1; h.sAcqType = 'CF';
            elseif acqType_idx == 2; h.sAcqType = 'PW';
            elseif acqType_idx == 3; h.sAcqType = 'PW tx apodized';
            else   error('undefined');
            end
            h.sAcqDir = [h.sPhtmDir num2str(acqType_idx) '. ' h.sAcqType '\'];
            
            %%--Angle 
            % ang_deg will be ignored if acqType_idx = 1 (CF)
            if acqType_idx == 1
                h.sDataDir = h.sAcqDir;
            elseif acqType_idx == 2||acqType_idx == 3
                h.sAngDir = [h.sAcqDir num2str(ang_deg) ' deg' '\' ];                  
                h.sDataDir = h.sAngDir;
            else   error('undefined'); 
            end
            
        end
    end
    
    
    methods (Access = public)
        function CreateFolders(h)

            % %%%%% Exp Select %%%%% %
            if ~exist(h.sExpDir,'dir')
                display(['creating folder:: ' h.sExpDir]);
                mkdir(h.sExpDir); 
            end

            % %%%%% Acquisition Type Select (acq dir is the final dir for 'CF') %%%%% %
            if ~exist(h.sAcqDir,'dir')
                display(['creating folder:: ' h.sAcqDir]);
                mkdir(h.sAcqDir); 
            end

            if h.acqType_idx == 2                
                % %%%%% Angle Select (angle dir is the final dir for 'PW') %%%%% %            
                if ~exist(h.sDataDir,'dir')
                    display(['creating folder:: ' h.sDataDir]);
                    mkdir(h.sDataDir); 
                end
            end
            
            % %%%%% Final Data Directory %%%%% %
            display(['creating RcvData folder at ' h.sDataDir]);

            % %%%%% Create folders %%%% %    
            sRFDir = [h.sDataDir 'RcvData\'];
            if(~isdir(sRFDir)); 
                mkdir(sRFDir); 
            end  
            sRFDir = [h.sDataDir 'BfData\'];
            if(~isdir(sRFDir)); 
                mkdir(sRFDir); 
            end  

        end
    end

    
    
end