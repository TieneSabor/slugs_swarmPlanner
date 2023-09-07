% Declare Transitions Variable
var 1..10: L1Q2; 
var 1..10: O1N2; 
var 1..10: N2P3; 
var 1..10: Q2S3; 
var 0..10: EXC; 
% Initial-Final Conditions
constraint  0 == 0; 
constraint L1Q2 +  0 == 4; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint O1N2 +  0 == 6; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
% Conservation Law
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint N2P3 +  0 == O1N2 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint Q2S3 +  0 == L1Q2 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
% Goal Function
constraint EXC == max( 0, L1Q2 - 10 ) + max( 0, O1N2 - 5 ) + max( 0, N2P3 - 5 ) + max( 0, Q2S3 - 10 ) +  0; 
solve minimize EXC;
output ["[step 1 2 ] KK: 0KL: 0KM: 0KO: 0LK: 0LL: 0LQ: \(L1Q2) MK: 0MM: 0MR: 0NN: 0NO: 0NP: 0NR: 0OK: 0ON: \(O1N2) OO: 0OP: 0PN: 0PO: 0PP: 0PQ: 0QL: 0QP: 0QQ: 0QS: 0RM: 0RN: 0RR: 0RS: 0SQ: 0SR: 0SS: 0\n"]; 
output ["[step 2 3 ] KK: 0KL: 0KM: 0KO: 0LK: 0LL: 0LQ: 0MK: 0MM: 0MR: 0NN: 0NO: 0NP: \(N2P3) NR: 0OK: 0ON: 0OO: 0OP: 0PN: 0PO: 0PP: 0PQ: 0QL: 0QP: 0QQ: 0QS: \(Q2S3) RM: 0RN: 0RR: 0RS: 0SQ: 0SR: 0SS: 0\n"]; 
