delete('test.m')
clear
pause
run('examples/CruiceControl1')
%load('ex_system');
Obj = Invariant(system,E,'test.m','feasible');
