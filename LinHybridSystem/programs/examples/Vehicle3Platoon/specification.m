tau1 = 4;
tau2 = 20;
Ac =   [0 1 0 0 0 0 0 0 0;
        0 0 -1 0 0 0 0 0 0;
        1.6050 4.8680 -3.5754 -0.8198 0.4270 -0.0450 -0.1942 0.3626 -0.0946;
        0 0 0 0 1 0 0 0 0;
        0 0 1 0 0 -1 0 0 0;
        0.8718 3.8140 -0.0754 1.1936 3.6258 -3.2396 -0.5950 0.1294 -0.0796;
        0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 1 0 0 -1;
        0.7132 3.5730 -0.0964 0.8472 3.2568 -0.0876 1.2726 3.0720 -3.1356];
 
Bc = [0 1 zeros(1,7)]';

Ac1_dis = expm(Ac*tau1);
fun = @(t) expm(Ac*(tau1-t))*Bc;
Bc1_dis = integral( fun,0,tau1,'ArrayValued',true );


Ac2_dis = expm(Ac*tau2);
fun = @(t) expm(Ac*(tau2-t))*Bc;
Bc2_dis = integral( fun,0,tau2,'ArrayValued',true );

An =   [0 1.0000 0 0 0 0 0 0 0;
        0 0 -1.0000 0 0 0 0 0 0;
        1.6050 4.8680 -3.5754 0 0 0 0 0 0;
        0 0 0 0 1.0000 0 0 0 0;
        0 0 1.0000 0 0 -1.0000 0 0 0;
        0 0 0 1.1936 3.6258 -3.2396 0 0 0;
        0 0 0 0 0 0 0 1.0000 0;
        0 0 0 0 0 1.0000 0 0 -1.0000;
        0.7132 3.5730 -0.0964 0.8472 3.2568 -0.0876 1.2726 3.0720 -3.1356];
    
Bn = Bc;

An1_dis = expm(An*tau1);
fun = @(t) expm(An*(tau1-t))*Bn;
Bn1_dis = integral( fun,0,tau1,'ArrayValued',true );


An2_dis = expm(An*tau2);
fun = @(t) expm(An*(tau2-t))*Bn;
Bn2_dis = integral( fun,0,tau2,'ArrayValued',true );

template = [1 0 0 0 0 0 0 0 0
            0 0 0 1 0 0 0 0 0
            0 0 0 0 0 0 1 0 0];
orth = null(template);
e = 1;

[Ec,~] = eig(Ac1_dis);
[En,~] = eig(An1_dis);

[Ec1,~] = eig(Ac1_dis*An2_dis);
[En1,~] = eig(An1_dis*Ac2_dis);

Apr = Ac; Apr([1,4,7],:) = zeros(3,9);
[Epr,~] = eig(Apr);

B_dis{1} = Bc1_dis; B_dis{2} = Bn1_dis;

system.NumLocs = 2;
system.dim = 9;

%% Location 1
system.Ptemplate{1} = zeros(1,9);
system.V{1} = [Ec,En,eye(9),Bc1_dis,Bn1_dis,Bc2_dis,Bn2_dis];
system.stay{1}.l = -inf*e;
system.stay{1}.u = inf*e;
system.map{1} = Ac1_dis;
system.input{1}.V = ones(9,1);
system.input{1}.c = zeros(9,1);
system.input{1}.s = 0;
system.input{1}.W = Bc1_dis;
system.input{1}.l =  -9;
system.input{1}.u = 1;

%% Location 2
system.Ptemplate{2} = system.Ptemplate{1};
system.V{2} = [Ec,En,eye(9),Bc1_dis,Bn1_dis,Bc2_dis,Bn2_dis];
system.stay{2}.l = -inf*e;
system.stay{2}.u = inf*e;
system.map{2} = An1_dis;
system.input{2}.V = ones(9,1);
system.input{2}.c = zeros(9,1);
system.input{2}.s = 0;
system.input{2}.W = Bn1_dis;
system.input{2}.l = -9;
system.input{2}.u = 1;

%% edge 1 from 1 to 2
system.edges{1}.loc1 = 1;
system.edges{1}.loc2 = 2;
system.edges{1}.u = inf*e;
system.edges{1}.l = -inf*e;
system.edges{1}.map = Ac2_dis;
system.edges{1}.input.V = ones(9,1);
system.edges{1}.input.c = zeros(9,1);
system.edges{1}.input.s = 0;
system.edges{1}.input.W = Bc2_dis;
system.edges{1}.input.l = -9;
system.edges{1}.input.u = 1;
system.edges{1}.lambda = (1-10^-5);
system.edges{1}.type = 1;

%% edge 2 from 2 to 1
system.edges{2}.loc1 = 2;
system.edges{2}.loc2 = 1;
system.edges{2}.u = inf*e;
system.edges{2}.l = -inf*e;
system.edges{2}.map = An2_dis;
system.edges{2}.input.V = ones(9,1);
system.edges{2}.input.c = zeros(9,1);
system.edges{2}.input.s = 0;
system.edges{2}.input.W = Bn2_dis;
system.edges{2}.input.l = -9;
system.edges{2}.input.u = 1;
system.edges{2}.lambda = (1-10^-5);
system.edges{2}.type = 1;

%% edge 3 from 1 to 2
system.edges{3}.loc1 = 1;
system.edges{3}.loc2 = 2;
system.edges{3}.u = inf*e;
system.edges{3}.l = -inf*e;
system.edges{3}.map = Ac2_dis;
system.edges{3}.input.V = ones(9,1);
system.edges{3}.input.c = zeros(9,1);
system.edges{3}.input.s = 0;
system.edges{3}.input.W = Bc2_dis;
system.edges{3}.input.l = 0;
system.edges{3}.input.u = 0;
system.edges{3}.lambda = (10^-5)/2.5;
system.edges{3}.type = 0;

%% edge 4 from 2
system.edges{4}.loc1 = 2;
system.edges{4}.loc2 = 1;
system.edges{4}.u = inf*e;
system.edges{4}.l = -inf*e;
system.edges{4}.map = An2_dis;
system.edges{4}.input.V = ones(9,1);
system.edges{4}.input.c = zeros(9,1);
system.edges{4}.input.s = 0;
system.edges{4}.input.W = Bn2_dis;
system.edges{4}.input.l = 0;
system.edges{4}.input.u = 0;
system.edges{4}.lambda = (1-10^-5);
system.edges{4}.type = 0;


%% Initial set
for i = 1:system.NumLocs
    system.initial{i}.V = ones(9,1);
    system.initial{i}.s = 0;
    system.initial{i}.c = zeros(9,1);
    system.initial{i}.W = B_dis{i};
    system.initial{i}.l = -9;
    system.initial{i}.u = 1;
end


%% safety property
system.safe{2}.T = -1*[zeros(1,6) 1 zeros(1,2)];
system.safe{2}.d = 1;
system.safe{1} = system.safe{2};

