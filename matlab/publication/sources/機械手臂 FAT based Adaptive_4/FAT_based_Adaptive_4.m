function Yp = FAT_based_Adaptive_4(t,Y,Kd,Kp,Gamma,B,P,l,T)
%% initialization
Yp = zeros(4+2*l,1);

[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();
[Xd,Xd_dot,Xd_ddot,qd,qd_dot,qd_ddot] = desired_trajectory_cartesian(t,l1,l2);

q = [Y(1);Y(2)]; % robot trajectory in the joint space
q_dot = [Y(3);Y(4)];

[X,X_dot] = forward_kenimatics(l1,l2,q,q_dot); % robot trajectory in the Cartesian space
[J,J_dot] = Jacobian_matrix(l1,l2,q,q_dot);
[D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q,q_dot,J,J_dot);

W_hat = reshape(Y(5:5+2*l-1),[l,2]); % estimation of weighting matrix W
Z = generate_basis(t,l,T); % basis vector z
Psi_hat = W_hat'*Z; % estimation of lumped uncertainty Psi

%% control method
e = X-Xd;
e_dot = X_dot - Xd_dot;
Xe = [e; e_dot]; % define a state vector
tau =  (Psi_hat +  Xd_ddot - Kd*e_dot - Kp*e);% controller
q_ddot = D\(-C*q_dot - G + tau); % closed-loop system
W_hat_dot = -Gamma*Z*Xe'*P*B; % update law

%% return data
Yp = [q_dot; q_ddot];
Yp(5:5+2*l-1) = reshape(W_hat_dot,[2*l,1]);
end