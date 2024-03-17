% Declare Transitions Variable
var 1..10: C1C2; 
var 1..10: E1C2; 
var 1..10: C2B3; 
var 1..10: C2D3; 
var 0..10: EXC; 
% Initial-Final Conditions
constraint  0 == 0; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint C2B3 +  0 == 4; 
constraint C1C2 +  0 == 5; 
constraint  0 == 0; 
constraint  0 == 0; 
constraint C2D3 +  0 == 6; 
constraint E1C2 +  0 == 5; 
constraint  0 == 0; 
% Conservation Law
constraint  0 ==  0; 
constraint  0 ==  0; 
constraint C2B3 + C2D3 +  0 == C1C2 + E1C2 +  0; 
constraint  0 ==  0; 
constraint  0 ==  0; 
% Goal Function
constraint EXC == max( 0, C1C2 - 5 ) + max( 0, E1C2 - 5 ) + max( 0, C2B3 - 5 ) + max( 0, C2D3 - 5 ) +  0; 
solve minimize EXC;
output ["[step 1 2 ] AA: 0 AB: 0 AC: 0 AD: 0 BA: 0 BB: 0 BC: 0 BE: 0 CA: 0 CB: 0 CC: \(C1C2) CD: 0 CE: 0 DA: 0 DC: 0 DD: 0 DE: 0 EB: 0 EC: \(E1C2) ED: 0 EE: 0 \n"]; 
output ["[step 2 3 ] AA: 0 AB: 0 AC: 0 AD: 0 BA: 0 BB: 0 BC: 0 BE: 0 CA: 0 CB: \(C2B3) CC: 0 CD: \(C2D3) CE: 0 DA: 0 DC: 0 DD: 0 DE: 0 EB: 0 EC: 0 ED: 0 EE: 0 \n"]; 
