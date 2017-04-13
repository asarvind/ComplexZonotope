clear 'system'
Initialization;
system.NumLocs = 1;
system.dim = 10;

%% Location (with no dynamics, identity map, zero input)
system.Ptemplate{1} = g_hat;
system.stay{1}.u = inf*[1;1];
system.stay{1}.l = -inf*[1;1];
% reference vectors
[E,~] = eig(A_hat); 
[E1,~] = eig(Aedge_hat1);
orth = null(g_hat);
% Prim template
[E2,~] = eig(Aedge_hat2);
[E3,~] = eig(Aedge_hat3);
orthvect = [E1(:,[1:7,10]),E2(:,10),E3(:,10),E(:,10),orth,orth*orth'*[E2,E]];
template = [E1(:,[8,9]),E(:,[1,4:6,8:10]),orthvect];
system.V{1} = template;
% dynamics
system.map{1} = 0*eye(10,10);
system.input{1}.V = ones(10,1);
system.input{1}.c = zeros(size(A_hat,1),1);
system.input{1}.s = 0;
system.input{1}.W = ones(10,1);
system.input{1}.l = 0;
system.input{1}.u = 0;

%% edge 1
system.edges{1}.loc1 = 1;
system.edges{1}.loc2 = 1;
system.edges{1}.l = [-v;-v];
system.edges{1}.u = [v;v];
system.edges{1}.map = A_hat;
system.edges{1}.input.V = ones(10,1);
system.edges{1}.input.c = zeros(size(A_hat,1),1);
system.edges{1}.input.s = 0;
system.edges{1}.input.W = B_hat;
system.edges{1}.input.l = -eps*ones(4,1);
system.edges{1}.input.u = eps*ones(4,1);

%% edge 2
system.edges{2}.loc1 = 1;
system.edges{2}.loc2 = 1;
system.edges{2}.l = [v;v];
system.edges{2}.u = [inf;inf];
system.edges{2}.map = Aedge_hat1;
system.edges{2}.input = system.edges{1}.input;
system.edges{2}.input.c = input1_hat+input2_hat;

%% edge 3
system.edges{3}.loc1 = 1;
system.edges{3}.loc2 = 1;
system.edges{3}.l = [-inf;-inf];
system.edges{3}.u = [-v;-v];
system.edges{3}.map = Aedge_hat1;
system.edges{3}.input = system.edges{1}.input;
system.edges{3}.input.c = -1*(input1_hat+input2_hat);

%% edge 4
system.edges{4}.loc1 = 1;
system.edges{4}.loc2 = 1;
system.edges{4}.l = [v;-inf];
system.edges{4}.u = [inf;-v];
system.edges{4}.map = Aedge_hat1;
system.edges{4}.input = system.edges{1}.input;
system.edges{4}.input.c = input1_hat-input2_hat;

%% edge 5
system.edges{5}.loc1 = 1;
system.edges{5}.loc2 = 1;
system.edges{5}.l = [-inf;v];
system.edges{5}.u = [-v;inf];
system.edges{5}.map = Aedge_hat1;
system.edges{5}.input = system.edges{1}.input;
system.edges{5}.input.c = input2_hat-input1_hat;

%% edge 6
system.edges{6}.loc1 = 1;
system.edges{6}.loc2 = 1;
system.edges{6}.l = [-inf;-v];
system.edges{6}.u = [-v;v];
system.edges{6}.map = Aedge_hat2;
system.edges{6}.input = system.edges{1}.input;
system.edges{6}.input.c = -input1_hat;

%% edge 7
system.edges{7}.loc1 = 1;
system.edges{7}.loc2 = 1;
system.edges{7}.l = [-v;-inf];
system.edges{7}.u = [v;-v];
system.edges{7}.map = Aedge_hat3;
system.edges{7}.input = system.edges{1}.input;
system.edges{7}.input.c = -input2_hat;

%% edge 8
system.edges{8}.loc1 = 1;
system.edges{8}.loc2 = 1;
system.edges{8}.l = [-v;v];
system.edges{8}.u = [v;inf];
system.edges{8}.map = Aedge_hat3;
system.edges{8}.input = system.edges{1}.input;
system.edges{8}.input.c = input2_hat;

%% edge 9
system.edges{9}.loc1 = 1;
system.edges{9}.loc2 = 1;
system.edges{9}.l = [v;-v];
system.edges{9}.u = [inf;v];
system.edges{9}.map = Aedge_hat2;
system.edges{9}.input = system.edges{1}.input;
system.edges{9}.input.c = input1_hat;




%% Initialization
system.initial{1}.V = eye(10,1);
system.initial{1}.s = 0;
system.initial{1}.c = zeros(10,1);
system.initial{1}.W = pinv(T_hat);
system.initial{1}.l = 0;
system.initial{1}.u = 0;

%% Safety bounds
system.safe{1}.T = T_hat;
system.safe{1}.d = 1;

