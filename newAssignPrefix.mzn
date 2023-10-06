% Declare Transitions Variable
var 1..10: N1N2; 
var 1..10: P1Q2; 
var 1..10: N2O3; 
var 1..10: Q2S3; 
var 1..10: O3P4; 
var 1..10: S3S4; 
var 0..10: EXC; 
% Initial-Final Conditions
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint N1N2 +  0 == 6; 
constraint  0 == 0; 
constraint P1Q2 +  0 == 4; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
% Conservation Law
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint N2O3 +  0 == N1N2 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint Q2S3 +  0 == P1Q2 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint O3P4 +  0 == N2O3 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint S3S4 +  0 == Q2S3 +  0; 
% Goal Function
constraint EXC == max( 0, N1N2 - 10 ) + max( 0, P1Q2 - 5 ) + max( 0, N2O3 - 10 ) + max( 0, Q2S3 - 10 ) + max( 0, O3P4 - 5 ) + max( 0, S3S4 - 10 ) +  0; 
solve minimize EXC;
output ["[step 1 2 ] KK: 0 KL: 0 KM: 0 KO: 0 LK: 0 LL: 0 LQ: 0 MK: 0 MM: 0 MR: 0 NN: \(N1N2) NO: 0 NR: 0 OK: 0 ON: 0 OO: 0 OP: 0 PO: 0 PP: 0 PQ: \(P1Q2) QL: 0 QP: 0 QQ: 0 QS: 0 RM: 0 RN: 0 RR: 0 RS: 0 SQ: 0 SR: 0 SS: 0 \n"]; 
output ["[step 2 3 ] KK: 0 KL: 0 KM: 0 KO: 0 LK: 0 LL: 0 LQ: 0 MK: 0 MM: 0 MR: 0 NN: 0 NO: \(N2O3) NR: 0 OK: 0 ON: 0 OO: 0 OP: 0 PO: 0 PP: 0 PQ: 0 QL: 0 QP: 0 QQ: 0 QS: \(Q2S3) RM: 0 RN: 0 RR: 0 RS: 0 SQ: 0 SR: 0 SS: 0 \n"]; 
output ["[step 3 4 ] KK: 0 KL: 0 KM: 0 KO: 0 LK: 0 LL: 0 LQ: 0 MK: 0 MM: 0 MR: 0 NN: 0 NO: 0 NR: 0 OK: 0 ON: 0 OO: 0 OP: \(O3P4) PO: 0 PP: 0 PQ: 0 QL: 0 QP: 0 QQ: 0 QS: 0 RM: 0 RN: 0 RR: 0 RS: 0 SQ: 0 SR: 0 SS: \(S3S4) \n"]; 
