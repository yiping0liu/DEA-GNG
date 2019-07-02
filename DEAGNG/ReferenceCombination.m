function [Ruq,net] = ReferenceCombination(Ru,net)
% Combine uniform reference vectors and nodes in GNG (Algorihtm 4)

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

    numNode = size(net.NodeS,1);    
    net.NodeP = net.NodeS;
    
    %% Map nodes to hyperplane
    for i = 1:numNode
       x = 1./sum(net.NodeS(i,:));
       net.NodeP(i,:) = net.NodeS(i,:)*x; 
    end    
    
    %% Average distance ammong nodes
    Distance1 = pdist2(net.NodeP,net.NodeP).*net.edge;   
    AvgDis = sum(Distance1(:))./sum(net.edge(:)); 
    
    %% Min distance among Ru
    Distance2 = pdist2(Ru,Ru);
    Distance2(Distance2 == 0) = 1;   
    MinDis = min(min(Distance2));  
  
    %% Choose the smaller one
    AvgDis = min(AvgDis,MinDis);

    %% Remove some reference vectors in Ru which are too close to nodes
    Distance3 = pdist2(Ru,net.NodeP);      
    Choose = all(Distance3 > AvgDis,2)==1;   
    Ruq = Ru(Choose,:);   
end

