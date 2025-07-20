function Pars = parsSet(dataTyp, stdGau)
%PARASET Summary of this function goes here
%   Detailed explanation goes here


switch dataTyp
    case 'video'
        % refeNbh.m parameters
        Pars.patSize = [6, 2];  % first is for R and C; second is for P
        Pars.nbhRad = [6, 4];  % bottom and right are nbhSize(1), forward is nbhSize(2)
        Pars.rptStep = [3, 1];
        Pars.cptStep = [1, 1];
        
        % simiMat.m parameters
        Pars.simNum = (Pars.patSize(1)^2)*Pars.patSize(2);  % matching number
        
        % MPCA.m parameters
        Pars.balChl = [3, 1, 1];
        Pars.balDiv = [128, Pars.simNum, Pars.simNum];
        
    case 'vlum'
        Pars.patSize = [4, 4];
        Pars.nbhRad = [5, 5];
        Pars.rptStep = [3, 3];
        Pars.cptStep = [1, 1];
        
        Pars.simNum = (Pars.patSize(1)^2)*Pars.patSize(2);
        
        Pars.balChl = 1.5;
        Pars.balWie = 48; 
        Pars.balDiv = 112; 
        
        if (stdGau>=35)
            Pars.patSize = [5, 5];
            Pars.simNum = (Pars.patSize(1)^2)*Pars.patSize(2);
            Pars.balChl = 3;
            Pars.balWie = 80;
            Pars.balDiv = 184;
        end
        
end

% image composition parameter, for pilot matching
Pars.balCps = 0.7;  % Im = balCps*Im + (1-balCps)*Jm
end

