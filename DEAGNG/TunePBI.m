function theta = TunePBI(net,eps)
% Adapt scalarizing functions by tune theta in PBI

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
 
    N = size(net.NodeS,1);
    theta=zeros(1,N);
    for i=1:N
        Neighbor = find(net.edge(i,:)==1);
        N1 = length(Neighbor);
        if N1==0
            theta(i)=Inf;
            continue;
        end
        EdgeVector = repmat(net.NodeS(i,:),N1,1)-net.NodeS(Neighbor,:);      
        flag=0;
        for j=1:N1
            if all(EdgeVector(j,:)==0)
                flag=1;
                break;
            end
        end
        if flag==1
            theta(i)=Inf;
            continue; 
        end
        Cosine = 1 - pdist2(net.NodeS(i,:),EdgeVector,'cosine');
        Cosine = max(Cosine);
        Angle = acos(Cosine);
        Angle =  Angle - eps;
        if Angle<0
            Angle=0;
        end
        theta(i)=1./tan(Angle);
        if theta(i)<0
            theta(i)=0;
        end
    end   
end