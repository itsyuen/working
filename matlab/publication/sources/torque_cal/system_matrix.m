function [D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q,q_dot,J,J_dot)
d11 = m1*(lc1^2) + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(q(2))) + I1 + I2;
d12 = m2*(lc2^2 + l1*lc2*cos(q(2))) + I2;
d21 = d12;
d22 = m2*(lc2^2) + I2;
c11 = -m2*l1*lc2*sin(q(2))*q_dot(2);
c12 = -m2*l1*lc2*sin(q(2))*(q_dot(1) + q_dot(2));
c21 = m2*l1*lc2*sin(q(2))*q_dot(1);
c22 = 0;
g1 = (m1*lc1 + m2*l1)*g*cos(q(1)) + m2*lc2*g*cos(q(1) + q(2));
g2 = m2*lc2*g*cos(q(1) + q(2));
D = [d11 d12; d21 d22];
C = [c11 c12; c21 c22];
G = [g1; g2];
Dx = J'\D/J;
Cx = J'\(C-D/J*J_dot)/J;
Gx = J'\G;
end