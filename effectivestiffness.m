l1 = 1.1E-3; 
w1 = 1.1E-3; 
h = 1.1E-3; 
L = 6.05E-3; 
l = 1.725E-3; 
E = 500E3; 
G = 179E3; 
C = 1.2; 
L2 = 4.6E-3; 

F = 5E-3; 
delta = F*(C*l/(l1*w1*G) + 4*l^3/(E*l1*w1^3) + 2*C*h/(l1*L*G) + 32*((l-h)^3-l^3)/(E*l1*L)); 
theta = F*(6*l^2/(E*l1*w1^3) + 48*((l+h)^2-l^2)/(E*l1*L^3)); 

delta2  = L2*sin(theta);
k = F/delta2; 