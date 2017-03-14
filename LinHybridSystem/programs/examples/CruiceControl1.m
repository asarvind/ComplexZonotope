%% Parameters
tau = 10^-3; % discretization parameter
amin = -1; % -1* maximum comfort deceleration
amax = 1; % maximum comfort acceleration
D = 0; % safety distance between cars

%% Matirces
Ac = [1 tau; 0 1]; 
Bc = [tau^2/2; tau];
Bd = -1*Bc;
v1 = [0 -0.4];
v2 = [-0.25 -1];
b = -0.25*(5+D)+(amax+amin)/2;

%% Hybrid system specification
system.NumLocs = 2;
system.dim = 6;

%% Location 1 dynamics 
% parallelotope template 
r1 = [zeros(1,2) 1 zeros(1,3)];
r2 = [zeros(1,5) 1];
r3 = [zeros(1,3) (v1-v2) 0];
r4 = [zeros(1,3) v1 0];
system.Ptemplate{1} = [r1; r2; r3; r4];
% map
system.map{1} = [zeros(2,3) eye(2,2) zeros(2,1); zeros(1,3) v1 0; zeros(2,3) Ac+Bc*v1 zeros(2,1); zeros(1,3) (v1-v2) 0];
% input
system.input{1}.V = [zeros(3,1); Bd; 0];
system.input{1}.s = ((amin-amax)/2);
system.input{1}.c = system.input{1}.V*(amax+amin)/2-[zeros(5,1); b];
system.input{1}.W = zeros(6,1);
system.input{1}.l = 0;
system.input{1}.u = 0;
% stay constriant
system.stay{1}.l = [amin; -inf; -inf; -inf];
system.stay{1}.u = [inf; b; inf; inf];

%% Location 2 dynamics
% parallelotope template
r1 = [zeros(1,2) 1 zeros(1,3)];
r2 = [zeros(1,5) 1];
r3 = [zeros(1,3) (v1-v2) 0];
r4 = [zeros(1,3) v2 0];
system.Ptemplate{2} = [r1; r2; r3; r4];
% map
system.map{2} = [zeros(2,3) eye(2,2) zeros(2,1); zeros(1,3) v2 0; zeros(2,3) Ac+Bc*v2 zeros(2,1); zeros(1,3) (v1-v2) 0];
% input
system.input{2} = system.input{1};
system.input{2}.c = system.input{2}.c + [zeros(2,1); b; Bc*b; 0];
% stay constraint
system.stay{2}.l = [amin; b; -inf; -inf];
system.stay{2}.u = [inf; inf; inf; inf];

%% edge 1 from location 1 to 1
system.edges{1}.loc1 = 1;
system.edges{1}.loc2 = 1;
system.edges{1}.l = [-inf -inf -inf -inf]';
system.edges{1}.u = [inf inf  b amin]';
system.edges{1}.map = [zeros(2,3) eye(2,2) zeros(2,1); zeros(1,6); zeros(2,3) Ac+Bc*v1 zeros(2,1); zeros(1,3) (v1-v2) 0];
system.edges{1}.input = system.input{1};
system.edges{1}.input.c = system.edges{1}.input.c+[zeros(2,1); amin; zeros(3,1)];
system.edges{1}.input.W = zeros(6,1);
system.edges{1}.input.l = 0;
system.edges{1}.input.u = 0;

%% edge 2 from location 2 to 2
system.edges{2}.loc1 = 2;
system.edges{2}.loc2 = 2;
system.edges{2}.l = [-inf -inf b -inf]';
system.edges{2}.u = [inf inf inf amin]';
system.edges{2}.map = [zeros(2,3) eye(2,2) zeros(2,1); zeros(1,6); zeros(2,3) Ac+Bc*v2 zeros(2,1); zeros(1,3) (v1-v2) 0];
system.edges{2}.input = system.input{1};
system.edges{2}.input.c = system.edges{2}.input.c+[zeros(2,1); amin; Bc*b; 0];
system.edges{2}.input.W = zeros(6,1);
system.edges{2}.input.l = 0;
system.edges{2}.input.u = 0;

%% edge 3 from location 1 to 2
system.edges{3}.loc1 = 1;
system.edges{3}.loc2 = 2;
system.edges{3}.l = [-inf -inf b amin]';
system.edges{3}.u = [inf inf inf inf]';
system.edges{3}.map = system.map{2};
system.edges{3}.input = system.input{2};
system.edges{3}.input.W = zeros(6,1);
system.edges{3}.input.l = 0;
system.edges{3}.input.u = 0;


%% edge 4 from location 1 to 2
system.edges{4}.loc1 = 1;
system.edges{4}.loc2 = 2;
system.edges{4}.l = [-inf -inf b -inf]';
system.edges{4}.u = [inf inf inf amin]';
system.edges{4}.map = system.edges{2}.map;
system.edges{4}.input = system.edges{2}.input;
system.edges{4}.input.W = zeros(6,1);
system.edges{4}.input.l = 0;
system.edges{4}.input.u = 0;

%% edge 5 from location 2 to 1
system.edges{5}.loc1 = 2;
system.edges{5}.loc2 = 1;
system.edges{5}.l = [-inf -inf -inf -inf]';
system.edges{5}.u = [inf inf b amin]';
system.edges{5}.map = system.map{1};
system.edges{5}.input = system.input{1};
system.edges{5}.input.W = zeros(6,1);
system.edges{5}.input.l = 0;
system.edges{5}.input.u = 0;

%% edge 6 from location 2 to 1
system.edges{6}.loc1 = 1;
system.edges{6}.loc2 = 2;
system.edges{6}.l = [-inf -inf -inf amin]';
system.edges{6}.u = [inf inf b inf]';
system.edges{6}.map = system.edges{1}.map;
system.edges{6}.input = system.edges{1}.input;
system.edges{6}.input.W = zeros(6,1);
system.edges{6}.input.l = 0;
system.edges{6}.input.u = 0;


%% Initial set
system.initial{1}.V = [zeros(3,1);Bd;0];
system.initial{1}.s = (amax-amin)/2;
system.initial{1}.c = [-20; 0; 0; (Ac*[-20; 0]+[(amax+amin)/2; (amax+amin)/2]); 0];
system.initial{1}.W = zeros(6,1);
system.initial{1}.l = 0;
system.initial{1}.u = 0;
system.initial{2} = [];

%% safety of system
system.safe{1}.T = [1 zeros(1,5); zeros(1,2) 1 zeros(1,3)];
system.safe{1}.d = [0; 10^10];
system.safe{2} = system.safe{1};

%% Direction vectors for synthesizing template
[E1,~] = eig(system.map{1});
[E2,~] = eig(system.map{2});
E = [E1 E2];

