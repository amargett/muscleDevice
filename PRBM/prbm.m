function rtip = prbm(theta)
    % Initial Variables for CAD model
    lfull = 50; 
    ltip = 3.25; 
    lrigid = 5; 
    ltop = 17; 
    theta0 = (90- 17/2)*pi/180;
    
    % Find Model Variables from CAD vars
    l = lfull - lrigid/2; 
    y = 0.8517; 
    lc = (1-y)*l/2;
    
    theta2 = theta + pi - theta0; 
    theta3 = theta + pi/2 -theta0; 
    
    ihat = [1 0]; jhat = [0 1]; 
    rA = (lc+ y*l*sin(theta))*jhat + y*l*cos(theta)*ihat; 
    rtip = rA + ltop/2*(cos(theta2)*ihat + sin(theta2)*jhat) ...
    + (ltip + lc + lrigid/2) * (cos(theta3)*ihat+ sin(theta3)*jhat);

    % r0 = (2*lc +ltip)*jhat -ltop/2*ihat + y*l*(cos(theta0)*ihat + sin(theta0)*jhat);
    % 
    % rtiprel = rtip-r0; 
end