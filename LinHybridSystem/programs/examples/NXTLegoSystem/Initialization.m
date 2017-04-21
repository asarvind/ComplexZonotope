clearvars all
tau = 4*10^-3;

Ap = [-162.1272879 0 162.1272879 -409.7184124 0 0;
      1 0 0 0 0 0;
      78.14959139 0 -78.14959139 269.6273393 0 0;
      0 0 1 0 0 0;
      0 0 0 0 -92.4134822 0;
      0 0 0 0 1 0];
  
Bp = [157.5798419 157.5798419;
      0 0;
      -75.95760351 -75.95760351;
      0 0;
      -51.3265211 51.3265211;
      0 0];
  
Cp =   [0 0 1 0 0 0;
        0 57.29577951 0 -57.29577951 0 -100.2676141;
        0 57.29577951 0 -57.29577951 0 100.2676141];

Dp = [0.08087, 0;
     0, 0.08087];
 
Ac =   [0 1.0 0 -1.0 0;
        0 0 0.996 0 0;
        0 0 -0.9999999525 0 0;
        0 0 0 0 0;
        0 0 0 49.99999763 -49.99999763];
 
 Bc =  [0 -0.00872664619 -0.00872664619 0 0;
        0 0 0 0.0003 0;
        0 0 0 0.07499999644 0;
        0.999999992 0 0 0 0;
        0 0.4363322888 0.4363322888 0 0];
    
    
        
Cc = [-5.529986398 -10.31821442 -13.529195 1161.547315 -679.1764238
      -5.529986398 -10.31821442 -13.529195 1161.547315 -679.1764238];
  
Dc =   [43.95659666 6.016975757 6.016975757 -0.004075058736 0.25;
        43.95659666 6.016975757 6.016975757 -0.004075058736 -0.25];

 sys = ss(Ac,Bc,Cc,Dc);
 sys_dis = c2d(sys,tau);
 Ac_dis = sys_dis.a;
 Bc_dis = sys_dis.b;
 
Bc1 = Bc(:,1:3); 
Bc1_dis = Bc_dis(:,1:3);
Bc2 = Bc(:,4:5); 
Bc2_dis = Bc_dis(:,4:5); 

Dc1 = Dc(:,1:3);
Dc2 = Dc(:,4:5);



M = [Ap zeros(6,5) Bp;
     Bc1*Cp Ac zeros(5,2)
     zeros(2,13)];
 
R = [eye(6,6) zeros(6,7);
     zeros(5,6) eye(5,5) zeros(5,2);
     Dp*Dc1*Cp Dp*Cc zeros(2,2)];
                                            
B = [zeros(6,4);
     Bc2_dis zeros(5,2);
     zeros(2,2) Dp*Dc2];
 A = R*expm(M*tau);
 
   A = round(A,6);
 
Aedge = A; Aedge([12,13],:) = zeros(2,13);

[V1,D] = jordan(A);
[V1,Apr] = cdf2rdf(V1,D); 
V2 = [pinv(V1([12:13,4],:)),null(V1([12:13,4],:),'r')];
Apr = V2\Apr*V2;
V = V1*V2;

Bpr = V\B;

A_hat = Apr(1:10,1:10);
B_hat = Bpr(1:10,:);
Aedge_hat1 = A_hat; 
Aedge_hat1(:,1:2) = zeros(10,2);

Aedge_hat2 = A_hat;
Aedge_hat2(:,1) = zeros(10,1);

Aedge_hat3 = A_hat;
Aedge_hat3(:,2) = zeros(10,1);

saturation = 100;
eps = 100;
v = 0.08087*saturation;



 input1_hat = A_hat(:,1)*v;
 input2_hat = A_hat(:,2)*v;

g_hat = [eye(2,2),zeros(2,8)];
T_hat = round(V(4,1:10));

 
