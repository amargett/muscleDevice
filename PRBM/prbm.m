function rtip = prbm(theta)
    % Initial Variables for CAD model
    l = 4.07; 
    lrigid = 1.68; 
    ltip = 0.68; 
    ltop = 1.68; 
    theta0 = 78.9*pi/180;
    
    % Find Model Variables from CAD vars
    y = 0.8517; 
    lc = (1-y)*l/2;
    
    theta2 = theta + pi - theta0; 
    theta3 = theta + pi/2 -theta0; 
    
    ihat = [1 0]; jhat = [0 1]; 
    rA = (lc+ y*l*sin(theta))*jhat + y*l*cos(theta)*ihat; 
    rtip = rA + ltop/2*(cos(theta2)*ihat + sin(theta2)*jhat) ...
    + (ltip + lc + lrigid/2) * (cos(theta3)*ihat+ sin(theta3)*jhat);
end