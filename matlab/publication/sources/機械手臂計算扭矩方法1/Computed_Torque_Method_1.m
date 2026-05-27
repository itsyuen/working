function Yp = Computed_Torque_Method_1(t, Y)
%% initialization
Yp = zeros(4,1);

Kv = diag([100, 100]);
Kp = diag([1000, 1000]);

[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();
[Xd,Xd_dot,Xd_ddot] = desired_trajectory_cartesian(t,l1,l2);

q = [Y(1);Y(2)]; % robot trajectory in the joint space
q_dot = [Y(3);Y(4)];

[X,X_dot] = forward_kenimatics(l1,l2,q,q_dot); % robot trajectory in the Cartesian space
[J,J_dot] = Jacobian_matrix(l1,l2,q,q_dot);
[D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q,q_dot,J,J_dot);


%% control method

tau = J'*(Cx*X_dot+Gx+Dx*(Xd_ddot-Kv*(X_dot-Xd_dot)-Kp*(X-Xd)));% controller

q_ddot = D\(-C*q_dot - G + tau); % closed-loop system


%% return data
Yp = [q_dot; q_ddot];
end