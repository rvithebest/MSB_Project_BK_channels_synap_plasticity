function ddt = eqn_solver(t,z,params)
BASAL= params.BASAL;
BASAL1= params.BASAL1;
BASAL2= params.BASAL2;
k1= params.k1;
k_1= params.k_1;
k_2= params.k_2;
k_4= params.k_4;
k6= params.k6;
k_6= params.k_6;
k7= params.k7;
k_8= params.k_8;
k9= params.k9;
k10= params.k10;
k_10= params.k_10;
k11= params.k11;
kSYN1= params.kSYN1;
k_12= params.k_12;
k8= params.k8;
lambda1= params.lambda1;
P= params.P;
YT= params.YT;
T= params.T;
%%%%%%%%%%%%%%%%
I2= params.I2;
X   = z(1);
XA  = z(2);
XAP = z(3);
YP  = z(4);
C1 = z(5);
C2 = z(6);
C3 = z(7);
C4 = z(8);
C5 = z(9);
C6 = z(10);
t_axis=params.t_axis;
Ca_conc=params.Ca;
if t<=15
   % Find the idx in t_axis that is closest to the current time t
   [~, idx] = min(abs(t_axis - t*1000)); % units of t_axis is ms, so multiply t by 1000
   CaCaM=Ca_conc(idx);
else 
    CaCaM = 2.94*I2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%
kSYN2 = 0.08*500*2.8;
k2 = 0.008*2.8;
k4 = 0.5*2.8;
k3 = 0.0075*2.8;
k5 = 0.022*2.8;

Y  = YT - YP;
lambda2 = lambda1;
lambda3 = lambda1;


ddt1  = kSYN2*C6 - k1*X*CaCaM + k_1*XA - lambda1*(X-BASAL);  % CaMKII Dynamics in inactive state

ddt2 = -k4*XA*XAP + k_4*C2 + k7*C3 + k1*X*CaCaM - k_1*XA - 2*k2*XA*XA + 2*k_2*C1 + k3*C1 - lambda2*(XA-BASAL1); % CaMKII Dynamics in active state

ddt3 =  k3*C1 - k4*XA*XAP + k_4*C2  + 2*k5*C2  - k6*XAP*P + k_6*C3 - k8*Y*XAP   + k_8*C4  + k9*C4  - lambda3*(XAP-BASAL2); % CaMKII Dynamics in active and phosphorylated state

ddt4 = k9*C4 - k10*YP*P + k_10*C5 - kSYN1*YP*T  + k_12*C6 + kSYN2*C6; % CPEBP- DYNAMICS

ddt5  = k2*XA*XA - k_2*C1 - k3*C1;

ddt6 = k4*XA*XAP - k_4*C2 - k5*C2;

ddt7 = k6*P*XAP - k_6*C3 - k7*C3;

ddt8 = k8*(Y)*XAP - k_8*C4 - k9*C4;

ddt9 = k10*YP*P - k_10*C5 - k11*C5;

ddt10 = kSYN1*YP*T - k_12*C6 - kSYN2*C6;

ddt = [ddt1;ddt2;ddt3;ddt4;ddt5;ddt6;ddt7;ddt8;ddt9;ddt10];

