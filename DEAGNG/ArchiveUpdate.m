function [Data,nND] = ArchiveUpdate(Data,N,Z1,Z2,Zmin)
% Update Input Signal Archive in DEAGNG

%--------------------------------------------------------------------------
% Copyright 2018-2019 Yiping Liu
% This is the code of DEA-GNG proposed in "Yiping Liu, Hisao Ishibuchi, 
% Naoki Masuyama, and Yusuke Nojima, Adapting reference vectors and 
% scalarizing functions by growing neural gas to handle irregular Pareto 
% fronts, IEEE Transactions on Evolutionary Computation, 2019, Early 
% Access, DOI: 10.1109/TEVC.2019.2926151".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete duplicate
    Data = unique(Data,'rows');

    %% Non-dominated sorting
    FrontNo = NDSort(Data,N);
    Data = Data((FrontNo == 1)',:);
    nND = size(Data,1);
    
    %% 
    if nND > N        
       Choose = LastSelection(Data,N,Z1,Z2,Zmin);
        Data = Data(Choose',:); 
        nND = N;      
    end      
    
end

function Choose = LastSelection(PopObj,K,Z1,Z2,Zmin)
% Select part of the solutions in the last front


    %% Initialize
    Z = [Z1;Z2];    
    [N,M]  = size(PopObj);
    NZ     = size(Z,1);
    NZ1     = size(Z1,1);
    
    Zmax = max(PopObj,[],1);
    PopObj = (PopObj - repmat(Zmin,N,1))./repmat(Zmax-Zmin,N,1);
%     PopObj = PopObj - repmat(Zmin,N,1); 
       
        
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Distance   = pdist2(PopObj,Z,'cosine');  
    % Associate each solution with its nearest reference point
    [d2,pi] = min(Distance',[],1);  
    
    %% Calculate the number of associated solutions of each reference point
%     rho = hist(pi,1:NZ); 
    
    cn = zeros(1,NZ);
    Choose  = false(1,N);
    Zchoose = true(1,NZ);
    while sum(Choose) < K
        Temp = find(Zchoose);
%         Jmin = find(cn(Temp)==min(cn(Temp)));
%         j = Temp(Jmin(randi(length(Jmin))));%randomly choose a ref point with min solutions among Z
        Jmin = cn(Temp) == min(cn(Temp));
        Temp = Temp(Jmin);
        Temp1 = Temp(Temp<=NZ1);
        if sum(Temp1) > 0
            j = Temp1(randi(length(Temp1)));
        else
            j = Temp(randi(length(Temp)));
        end
        
        
        
        I = find(Choose==0 & pi==j);%the solutions relate to j
        % Then delete one solution associated with this reference point
        if ~isempty(I)
            if cn(j) == 0
                [~,s] = min(d2(I));
%                 s = randi(length(I));
            elseif cn(j) < M
                [~,s] = max(d2(I));
%                 s = randi(length(I));
            else
                s = randi(length(I));
%                 [~,s] = max(d2(I)); 
            end
            Choose(I(s)) = true;
            cn(j) = cn(j) + 1;
        else
            Zchoose(j) = false;
        end
        
    end
    
%     %% Environmental selection
%     Choose  = false(1,N);
%     Zchoose = true(1,NZ);
%     % Select K solutions one by one
%     while sum(Choose) < K
%         % Select the least crowded reference point
%         Temp = find(Zchoose);      
%         Jmin = find(rho(Temp)==min(rho(Temp)));
%         j = Temp(Jmin(randi(length(Jmin))));%randomly choose a ref point with min solutions among Z
%         I = find(Choose==0 & pi==j);%the solutions in pop2 relate to j
%         % Then select one solution associated with this reference point
%         if ~isempty(I)
%             s = randi(length(I));
%             Choose(I(s)) = true;
%             rho(j) = rho(j) + 1;
%         else
%             Zchoose(j) = false;
%         end
%     end   
    
end

