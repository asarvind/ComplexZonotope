 delete('test.m')
pause
clear 'system'
clearvars all
clear -global
clear -global


run('examples/academic_ex/specification')

rng(0,'twister')
timervar = tic;
Obj = Invariant(system,'test.m','optimize','mosek','default');
fprintf( 'Computation time for verifying property = %f \n\n\n',toc(timervar) );
