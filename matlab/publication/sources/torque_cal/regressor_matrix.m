function [Y,p] = regressor_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q,q_dot,q_ddot)
p = [m1*lc1^2; m2*l1^2; m2*lc2^2; m2*l1*lc2; I1; I2; m1*lc1*g; m2*l1*g; m2*lc2*g];
y14 = 2*cos(q(2))*q_ddot(1)+cos(q(2))*q_ddot(2)-2*sin(q(2))*q_dot(1)*q_dot(2)-sin(q(2))*q_dot(2)^2;
y24 = cos(q(2))*q_ddot(1)+sin(q(2))*q_dot(1)^2;
Y = [q_ddot(1) q_ddot(1) q_ddot(1)+q_ddot(2) y14 q_ddot(1) q_ddot(1)+q_ddot(2) cos(q(1)) cos(q(1)) cos(q(1)+q(2));
    0 0 q_ddot(1)+q_ddot(2) y24 0 q_ddot(1)+q_ddot(2) 0 0 cos(q(1)+q(2))];
end