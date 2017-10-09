parameters;
system.dim = 2;
system.NumLocs = 4;
system.K = eye(2);

[E1,~] = eig(M{1});
[E2,~] = eig(M{3});
[E3,~] = eig(M{1}*M{3});
[E4,~] = eig(M{3}*M{1});

system.V = [E1,E2,E3,E4,eye(2),M{1},M{3},M{1}*M{3},M{3}*M{1}];
system.input.V = eye(2);
system.Vinit = ones(2,1);

c = input('\n enter input bound\n');

%% Locations
for i = 1:system.NumLocs
    system.stay{i}.b = u{i};
    system.stay{i}.a = l{i};
    system.map{i} = M{i};
    system.input.b{i} = [c;c];
    system.input.a{i} = -1*[c;c];
end

%% Edges
for i = 1:system.NumLocs
    for j = 1:system.NumLocs
        if i~=j
            system.edges{i}.loc1 = i;
            system.edges{i}.loc2 = j;
            system.edges{i}.b = inf;
            system.edges{i}.a = -inf;
            system.edges{i}.map = system.map{i};
            system.edges{i}.input.a = system.input.a{i};
            system.edges{i}.input.b = system.input.b{i};
        end
    end
end

%% Initial set
for i = 1:system.NumLocs
    system.initial{i}.a = 0;
    system.initial{i}.b = 0;
end

%% Safety
for i = 1:system.NumLocs
    system.safe{i}.T = [1 0];
    system.safe{i}.d = 1;
end
