function [J,J_dot] = Jacobian_matrix(l1,l2,q,q_dot)
J = [-l1*sin(q(1))-l2*sin(q(1)+q(2)) -l2*sin(q(1)+q(2)); l1*cos(q(1))+l2*cos(q(1)+q(2)) l2*cos(q(1)+q(2))];
J_dot = [-l1*cos(q(1))*q_dot(1)-l2*cos(q(1)+q(2))*(q_dot(1)+q_dot(2)) -l2*cos(q(1)+q(2))*(q_dot(1)+q_dot(2));
    -l1*sin(q(1))*q_dot(1)-l2*sin(q(1)+q(2))*(q_dot(1)+q_dot(2)) -l2*sin(q(1)+q(2))*(q_dot(1)+q_dot(2))];
end