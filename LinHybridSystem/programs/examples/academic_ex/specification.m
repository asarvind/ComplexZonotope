parameters;

system.NumLocs = 2;
system.dim = 3;
[E{1},~] = eig(A{1}); [E{2},~] = eig(A{2});
E{3} = eig(A{1}*A{2}); E{4} = eig(A{2}*A{1});
ptemplate = [1 0 0];
orth = null(ptemplate);
prim_template = [E{1} E{2} E{3} E{4} orth orth*orth'*[E{1} E{2} E{3} E{4}]];
b{1} = [-inf 0]; b{2} = [0 inf];
a = input('\n bound on input magnitude\n');

%% Locations
for i = 1:system.NumLocs
    system.Ptemplate{i} = ptemplate;
    system.V{i} = prim_template;
    system.stay{i}.l = b{i}(1);
    system.stay{i}.u = b{i}(2);
    system.map{i} = A{i};
    system.input{i}.V = ones(3,1);
    system.input{i}.s = 0;
    system.input{i}.c = zeros(3,1);
    system.input{i}.W = B{i};
    system.input{i}.u = a;
    system.input{i}.l = -a;
end

for i = 1:system.NumLocs
    for j = 1: system.NumLocs
        if i~=j
            system.edges{i}.loc1 = i;
            system.edges{i}.loc2 = j;
            system.edges{i}.u = inf;
            system.edges{i}.l = -inf;
            system.edges{i}.map = system.map{i};
            system.edges{i}.input = system.input{i};
        end
    end
end

%% initial condition
for i = 1: system.NumLocs
    system.initial{i}.V = ones(3,1);
    system.initial{i}.c = zeros(3,1);
    system.initial{i}.s = 0;
    system.initial{i}.W = ones(3,1);
    system.initial{i}.l = 0;
    system.initial{i}.u = 0;
end

%% safety
for i = 1:system.NumLocs
    system.safe{i}.T = [1 0 0];
    system.safe{i}.d = 1;
end