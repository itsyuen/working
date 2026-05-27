function [Yx] = regressor_matrix(l1,l2,X,X_dot,v,v_dot)

q = zeros(2,1);
q(2) = acos((X(1)^2 + X(2)^2 - l1^2 - l2^2)/(2*l1*l2));
q(1) = atan(X(2)/X(1)) - atan(l2*sin(q(2))/(l1 + l2*cos(q(2))));

J = [-l1*sin(q(1))-l2*sin(q(1)+q(2)) -l2*sin(q(1)+q(2)); l1*cos(q(1))+l2*cos(q(1)+q(2)) l2*cos(q(1)+q(2))];
q_dot = zeros(2,1);
q_dot = J\X_dot;

Yx = [v_dot(1) v_dot(1)+v_dot(2) (2*v_dot(1)+v_dot(2))*cos(q(2))-(v(1)*q_dot(2)+v(2)*(q_dot(1)+q_dot(2)))*sin(q(2)) cos(q(1)) cos(q(1)+q(2));
    0 v_dot(1)+v_dot(2) v_dot(1)*cos(q(2))+v(1)*sin(q(2))*q_dot(1) 0 cos(q(1)+q(2))];
end