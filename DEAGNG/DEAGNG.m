function DEAGNG(Global)
% <algorithm> <A-G>
% Adapting reference vectors and scalarizing functions by growing neural 
% gas to handle irregular Pareto fronts 
% aph --- 0.1 --- alpha
% eps --- 0.314 --- epsilon

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
    
    %% Parameter setting
    [aph,eps] = Global.ParameterSet(.1,.314); 
                           
    %% DEA Initialization
    [Ru,Global.N] = UniformPoint(Global.N,Global.M); % Uniform reference vectors
    Population    = Global.Initialization(); % Random population
    Zmin          = min(Population.objs,[],1); %Ideal Point  
    AS            = []; %Input Signal Archive
    Ruq           = Ru; %Refernce vectors in Ru for selection
    [FrontNo,~] = NDSort(Population.objs,Global.N); % Fitness for the first mating selection
    crd = zeros(1,Global.N); % Fitness for the first mating selection
    
    %% GNG Initialization
    ArchiveSize = Global.M*Global.N;    % Size of Input Signal Archive
    NoG = aph*Global.maxgen;            % Number of generations of Not Training GNG
    GNGnet.maxIter = 1;                 % Number of iterations to train GNG per Generation  
    GNGnet.maxAge = Global.N;           % Maximum cluster age     
    GNGnet.maxNode = Global.N;          % Max number of nodes
    GNGnet.lambda = 0.2*Global.N;       % Cycle for topology reconstruction   
    GNGnet.hp = [];                     % Hit point of node 
    GNGnet.maxHP = 2*ArchiveSize;       % Max HP of node 
    GNGnet.Node = [];                   % Node
    GNGnet.NodeS = [];                  % Expanded node 
    GNGnet.NodeP = [];                  % Node mapped to hyperplane 
    GNGnet.Err = [];                    % Error 
    GNGnet.edge = zeros(2,2);           % Edge between nodes 
    GNGnet.age = zeros(2,2);            % Age of edge 
    GNGnet.epsilon_a = 0.2;             % Learning coefficient
    GNGnet.epsilon_nb = 0.01;           % Learning coefficient of neighbor
    GNGnet.alpha = 0.5;                 % Nodes r1max and r2max error reduction constant
    GNGnet.delta = 0.9;                 % Error reduction coefficient 

    %% Optimization
    while Global.NotTermination(Population)        
               
        MatingPool = TournamentSelection(2,Global.N,FrontNo,crd);
        Offspring  = Global.Variation(Population(MatingPool));       
        Zmin       = min([Zmin;Offspring.objs],[],1);
                
        %% GNG-based adaptation
        if Global.gen <= Global.maxgen - NoG   
              
            % Input Signal Archive Update
            AS = ArchiveUpdate([AS;Offspring.objs],ArchiveSize,Ruq,GNGnet.NodeS,Zmin);
            nAS = length(AS);      
            
            % GNG Update (and Algorithm 3)
            GNGnet.maxNode = min(Global.N,floor(nAS/2)); % paramter reset 
            GNGnet.maxHP = 2*nAS; % paramter reset
            GNGnet = GNGUpdate(AS,GNGnet);
            
            % Reference Vector Adaptation (Algorithm 4) 
            if size(GNGnet.NodeS,1)>2
                [Ruq,GNGnet] = ReferenceCombination(Ru,GNGnet);
            end
            
            % Scalarizing Function Adaptation
            theta = TunePBI(GNGnet,eps); % Tune theta in PBI function                    
        end
                            
       %% Environmental Selection                 
       [Population,FrontNo,crd] = ESelection([Population,Offspring],Global.N,Ruq,GNGnet.NodeS,theta,Zmin);  
       
    end           
end