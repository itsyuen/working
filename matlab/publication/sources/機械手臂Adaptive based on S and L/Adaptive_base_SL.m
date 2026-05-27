function Yp = Adaptive_base_SL(t,Y,lambda,gamma,Kd)
%% initialization
Yp = zeros(9,1);

[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();
[Xd,Xd_dot,Xd_ddot,qd,qd_dot,qd_ddot] = desired_trajectory_cartesian(t,l1,l2);

q = [Y(1);Y(2)]; % robot trajectory in the joint space
q_dot = [Y(3);Y(4)];

p_hat = [Y(5);Y(6);Y(7);Y(8);Y(9)];

[X,X_dot] = forward_kenimatics(l1,l2,q,q_dot); % robot trajectory in the Cartesian space
[J,J_dot] = Jacobian_matrix(l1,l2,q,q_dot);
[D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q,q_dot,J,J_dot);


%% control method
e = X-Xd;
e_dot = X_dot - Xd_dot;
s = e_dot + lambda*e;
v = X_dot - lambda*e;
v_dot = Xd_ddot - lambda*e_dot;

[Yx] = regressor_matrix(l1,l2,X,X_dot,v,v_dot);

tau = J'*(Yx*p_hat-Kd*s);% controller

q_ddot = D\(-C*q_dot - G + tau); % closed-loop system

p_hat_dot = -gamma*Yx'*s;

%% return data
Yp = [q_dot; q_ddot;p_hat_dot];
end