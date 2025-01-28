clear all; close all; clc;

% x vector is of the form: [L1; L2; w1; w2];
x0 = [0.005; 0.005; 0.005; 0.005];
lb = [0.001; 0.001; 0.001; 0.001];
ub = [0.01; 0.01; 0.01; 0.01];

x_max = fmincon(@(x) get_diff(x(1), x(2), x(3), x(4)), x0, [], [], [], [], lb, ub);

% Display optimal values
disp('Optimal solution :');
disp(x_max);

k1 = get_stiffness(x_max(1), x_max(2), x_max(3), x_max(4), 0)
k2 = get_stiffness(x_max(1), x_max(2), x_max(3), x_max(4), pi/2)

%%
function diff = get_diff(L1, L2, w1, w2)
    k1 = get_stiffness(L1, L2, w1, w2, 0); 
    k2 = get_stiffness(L1, L2, w1, w2, pi/2); 
    diff = k2 -k1; 
end
function k = get_stiffness(L1, L2, w1, w2, theta)
    % resulting angles
    theta_A = -(pi/2-theta); theta_B = pi/2-theta; 
    theta_E1 = 0; theta_E2 = 0; theta_F1 = 0; theta_F2 = 0; 

    % define beams and their properties
    % 'nodes' array order matters -- first entry is left node of beam
    % pre-rotation. second entry is the right node pre-rotation.
    beam_A = create_beam(L1, w1, theta_A, [1,3]); 
    beam_B = create_beam(L1, w1, theta_B, [3,2]); 
    beam_E1 = create_beam(L2, w2, theta_E1, [4,3]); 
    beam_E2 = create_beam(L2, w2, theta_E2, [3,5]); 
    beam_F1 = create_beam(L2, w2, theta_F1, [4,6]); 
    beam_F2 = create_beam(L2, w2, theta_F2, [6,5]); 
    
    all_beams = [beam_A, beam_B, beam_E1, beam_E2, beam_F1, beam_F2];
    
    %symbolic global stiffness matrix:
    K = generate_global_stiffness_matrix(all_beams, 6); %11 nodes, each with 3 d.o.f. since elements are beam elements
    alldofs = 1:18; 
    dofspec = [1 2 4 5]; 
    doffree = alldofs; doffree(dofspec) = []; 
    u = zeros(18, 1); 
    f = zeros(18, 1); 
    
    f(17) = 5; % add 5 N load
    
    u(doffree) = K(doffree, doffree)\f(doffree) ;
    f(dofspec) = K(dofspec,:)*u ;
    
    k = f(17)/u(17);
end

function beam = create_beam(L, w, theta, nodes)
    t = 0.125*0.0254; %[m] thickness of LDPE
    E = 2E8;
    I = t*w^3/12; 
    A = t*w; 
    beam = struct('E', E, 'L', L, 'w', w, 'theta', theta,'nodes', nodes, 'A', A, 'I', I);
end

function global_k_mat = generate_global_stiffness_matrix(beams, n)
   % beam is a struct representing a beam element with properties:
   % A, E, L, W, theta, nodes. n is the total number of nodes
   % returns a local stiffness matrix, k_mat, for the beam element input
   global_k_mat = sym(zeros(3*n, 3*n));
   for beam=beams
       A = beam.A; E = beam.E; L = beam.L; theta = beam.theta; I = beam.I;
       first_node = beam.nodes(1); second_node = beam.nodes(2);
       
       k_mat = [A*E/L,    0,              0,          -A*E/L,     0,              0;
              0,        12*E*I/L^3,     6*E*I/L^2,  0,          -12*E*I/L^3,    6*E*I/L^2;
              0,        6*E*I/L^2,      4*E*I/L,    0,          -6*E*I/L^2,     2*E*I/L;
              -A*E/L,   0,              0,          A*E/L,      0,              0;
              0,        -12*E*I/L^3,   -6*E*I/L^2,  0,           12*E*I/L^3,   -6*E*I/L^2;
              0,        6*E*I/L^2,      2*E*I/L,    0,          -6*E*I/L^2,     4*E*I/L;
              ];
       % create rotation matrix % entries 
       c = cos(theta); s = sin(theta);
       rot_mat = [c,  -s,  0,  0,  0, 0;
                  s,   c,  0,  0,  0, 0;
                  0,   0,  1,  0,  0, 0;
                  0,   0,  0,  c,  s, 0;
                  0,   0,  0, -s,  c, 0;
                  0,   0,  0,  0,  0, 1];
       % apply rotation
       k_mat = rot_mat'*k_mat*rot_mat;

       local_k = zeros(size(global_k_mat)); 
       ind1 = (first_node-1)*3+1; ind2 = ind1+2; 
       ind3 = (second_node-1)*3+1; ind4 = ind3+2;
       local_k(ind1:ind2, ind1:ind2) = k_mat(1:3, 1:3); 
       local_k(ind1:ind2, ind3:ind4) = k_mat(1:3, 4:6); 
       local_k(ind3:ind4, ind3:ind4) = k_mat(4:6, 4:6);
       local_k(ind3:ind4, ind1:ind2) = k_mat(4:6, 1:3); 
       
       global_k_mat = global_k_mat + local_k; 
   end
end