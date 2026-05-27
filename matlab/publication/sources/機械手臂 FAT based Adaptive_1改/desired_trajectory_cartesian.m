% Desired trajectory in the Cartesian space

function [Xd,Xd_dot,Xd_ddot,qd,qd_dot,qd_ddot] = desired_trajectory_cartesian(t,l1,l2)
    Xd = [1+0.1*sin(0.1*2*pi*t) 0.5-0.1*cos(0.1*2*pi*t)]';
    Xd_dot = [0.02*pi*cos(0.2*pi*t) 0.02*pi*sin(0.2*pi*t)]';
    Xd_ddot = [-0.004*pi*pi*sin(0.2*pi*t) 0.004*pi*pi*cos(0.2*pi*t)]';

qd = zeros(2,1);
qd(2) = acos((Xd(1)^2 + Xd(2)^2 - l1^2 - l2^2)/(2*l1*l2));
qd(1) = atan(Xd(2)/Xd(1)) - atan(l2*sin(qd(2))/(l1 + l2*cos(qd(2))));

J = [-l1*sin(qd(1))-l2*sin(qd(1)+qd(2)) -l2*sin(qd(1)+qd(2)); l1*cos(qd(1))+l2*cos(qd(1)+qd(2)) l2*cos(qd(1)+qd(2))];
qd_dot = J\Xd_dot;

J_dot = [-l1*cos(qd(1))*qd_dot(1)-l2*cos(qd(1)+qd(2))*(qd_dot(1)+qd_dot(2)) -l2*cos(qd(1)+qd(2))*(qd_dot(1)+qd_dot(2));
    -l1*sin(qd(1))*qd_dot(1)-l2*sin(qd(1)+qd(2))*(qd_dot(1)+qd_dot(2)) -l2*sin(qd(1)+qd(2))*(qd_dot(1)+qd_dot(2))];
qd_ddot =J\(Xd_ddot-J_dot*qd_dot) ;

end