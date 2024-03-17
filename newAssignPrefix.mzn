% Declare Transitions Variable
var 1..10: A1C2; 
var 1..10: B1C2; 
var 1..10: C2C3; 
var 1..10: C2E3; 
var 0..10: EXC; 
% Initial-Final Conditions
constraint A1C2 +  0 == 7; 
constraint B1C2 +  0 == 3; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
% Conservation Law
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint C2C3 + C2E3 +  0 == A1C2 + B1C2 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
% Goal Function
constraint EXC == max( 0, A1C2 - 5 ) + max( 0, B1C2 - 5 ) + max( 0, C2C3 - 5 ) + max( 0, C2E3 - 5 ) +  0; 
solve minimize EXC;
output ["[step 1 2 ] AA: 0 AB: 0 AC: \(A1C2) AD: 0 BA: 0 BB: 0 BC: \(B1C2) BE: 0 CA: 0 CB: 0 CC: 0 CD: 0 CE: 0 DA: 0 DC: 0 DD: 0 DE: 0 EB: 0 EC: 0 ED: 0 EE: 0 \n"]; 
output ["[step 2 3 ] AA: 0 AB: 0 AC: 0 AD: 0 BA: 0 BB: 0 BC: 0 BE: 0 CA: 0 CB: 0 CC: \(C2C3) CD: 0 CE: \(C2E3) DA: 0 DC: 0 DD: 0 DE: 0 EB: 0 EC: 0 ED: 0 EE: 0 \n"]; 
