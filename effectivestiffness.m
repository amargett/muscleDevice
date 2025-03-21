clear; close all

N=5; 

E_vals = linspace(0.5, 1.5, N); % MPa

K_vals = zeros(size(E_vals)); 

for i =1:N
    K_vals(i) = get_skel_k(E_vals(i)); 
end

plot(E_vals, K_vals)
xlabel('E [MPa]'); ylabel('K [N/m]'); 


l2_vals = linspace(1,5, N); %mm

figure(2)
for E = E_vals
    k_vals = zeros(size(l2_vals)); 
    for i =1:N
        k_vals(i) = get_new_k(l2_vals(i), E);
    end
    plot(l2_vals, k_vals); 
    hold on
end
legend
hold off

function k = get_new_k(L2, E)
L2 = L2*10^-3; 
E = E*10^6; 
l1 = 1.1E-3; 
w1 = 1.1E-3; 
h = 1.1E-3; 
L = 6.05E-3; 
l = 0.625E-3; 
v = 0.4; 
G = E*(2*(1+v)); 
C = 1.2; 
F = 3E-4; % 300 uN active tension

I = l1*w1^3/12;

delta_shear = F*(l+h)/(l1*w1*G);
delta_F = F*(l+h)^2/(6*E*I)*(3*L2-h-l) ;

delta_tip = delta_shear + delta_F; 

k = F/delta_tip; 
end



function k = get_skel_k(E)
E = E*10^6; 
l1 = 1.1E-3; 
w1 = 1.1E-3; 
h = 1.1E-3; 
L = 6.05E-3; 
l = 0.625E-3; 
v = 0.4; 
G = E*(2*(1+v)); 
C = 1.2; 
L2 = 4.6E-3; 

F = 3E-4; % 300 uN active tension

% Pure Bending
I = l1*w1^3/12;
M = F*(l+ h/2); 
theta = M*L/(2*E*I); 
delta_theta = L2*theta ; % small angle

delta_shear = F*l/(l1*w1*G) + 2*F*C*h/(l1*L*G);
delta_F = F*l^2/(6*E*I)*(3*(L2-h)-l) ;

delta_tip = delta_shear + delta_F + delta_theta; 

k = F/delta_tip; 
end