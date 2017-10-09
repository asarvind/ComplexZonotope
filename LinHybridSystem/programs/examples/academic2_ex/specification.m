parameters;
system.dim = 2;
system.NumLocs = 4;

[E1,~] = eig(M{1});
[E2,~] = eig(M{3});
[E3,~] = eig(M{1}*M{3});
[E4,~] = eig(M{3}*M{1});

prim_template = [E1,E2,E3,E4];
Ptemplate = eye(2);

%% Locations
for i = 1:system.NumLocs
    system.V{i} = prim_template;
    system.Ptemplate{i} = Ptemplate;
    system.stay{i}.u = u{i};
    system.stay{i}.l = l{i};
    system.map{i} = M{i};
    system.input{i}.V = ones(2,1);
    system.input{i}.c = zeros(2,1);
    system.input{i}.s = 0;
    system.input{i}.W = eye(2);
    system.input{i}.u = [0.2;0.2];
    system.input{i}.l = -1*[0.2;0.2];
end

%% Edges
for i = 1:system.NumLocs
    for j = 1:system.NumLocs
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

%% Initial set
for i = 1:system.NumLocs
    system.initial{i}.V = ones(2,1);
    system.initial{i}.c = zeros(2,1);
    system.initial{i}.s = 0;
    system.initial{i}.W = eye(2);
    system.initial{i}.l = [0; 0];
    system.initial{i}.u = [0; 0];
end

%% Safety
for i = 1:system.NumLocs
    system.safe{i}.T = [0 1;0 -1];
    system.safe{i}.d = [0.36;0.36];
end
