% Declare Transitions Variable
var 1..10: M1R2; 
var 1..10: Q1Q2; 
var 1..10: Q2P3; 
var 1..10: R2S3; 
var 0..10: EXC; 
% Initial-Final Conditions
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint M1R2 +  0 == 4; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint Q2P3 +  0 == 6; 
constraint Q1Q2 +  0 == 6; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint R2S3 +  0 == 4; 
% Conservation Law
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint Q2P3 +  0 == Q1Q2 +  0; 
constraint R2S3 +  0 == M1R2 +  0; 
constraint  0 ==  0; 
% Goal Function
constraint EXC == max( 0, M1R2 - 10 ) + max( 0, Q1Q2 - 10 ) + max( 0, Q2P3 - 10 ) + max( 0, R2S3 - 10 ) +  0; 
solve minimize EXC;
output ["[step 1 2 ] KK: 0KL: 0KM: 0KO: 0LK: 0LL: 0LQ: 0MK: 0MM: 0MR: \(M1R2) NN: 0NO: 0NP: 0NR: 0OK: 0ON: 0OO: 0OP: 0PN: 0PO: 0PP: 0PQ: 0QL: 0QP: 0QQ: \(Q1Q2) QS: 0RM: 0RN: 0RR: 0RS: 0SQ: 0SR: 0SS: 0\n"]; 
output ["[step 2 3 ] KK: 0KL: 0KM: 0KO: 0LK: 0LL: 0LQ: 0MK: 0MM: 0MR: 0NN: 0NO: 0NP: 0NR: 0OK: 0ON: 0OO: 0OP: 0PN: 0PO: 0PP: 0PQ: 0QL: 0QP: \(Q2P3) QQ: 0QS: 0RM: 0RN: 0RR: 0RS: \(R2S3) SQ: 0SR: 0SS: 0\n"]; 
