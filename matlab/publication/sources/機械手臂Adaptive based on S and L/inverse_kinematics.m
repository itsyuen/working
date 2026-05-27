function q = inverse_kinematics(l1,l2,X)
q = zeros(2,1);
q(2) = acos((X(1)^2 + X(2)^2 - l1^2 - l2^2)/(2*l1*l2));
q(1) = atan(X(2)/X(1)) - atan(l2*sin(q(2))/(l1 + l2*cos(q(2))));
end