function Yp = FAT_based_Adaptive_1(t,Y,lambda,Gamma_D,Gamma_C,Gamma_g,Kd,l,T)
%% initialization
Yp = zeros(4+10*l,1);

[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();
[Xd,Xd_dot,Xd_ddot,qd,qd_dot,qd_ddot] = desired_trajectory_cartesian(t,l1,l2);

q = [Y(1);Y(2)]; % robot trajectory in the joint space
q_dot = [Y(3);Y(4)];

[X,X_dot] = forward_kenimatics(l1,l2,q,q_dot); % robot trajectory in the Cartesian space
[J,J_dot] = Jacobian_matrix(l1,l2,q,q_dot);
[D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q,q_dot,J,J_dot);

Z = generate_basis(t,l,T); % basis vector z

WD1_hat = reshape(Y(5:5+l-1),[l,1]); % estimation of weighting matrix W
WD2_hat = reshape(Y(5+l:5+2*l-1),[l,1]);
WD3_hat = reshape(Y(5+2*l:5+3*l-1),[l,1]);
WD4_hat = reshape(Y(5+3*l:5+4*l-1),[l,1]);

WC1_hat = reshape(Y(5+4*l:5+5*l-1),[l,1]);
WC2_hat = reshape(Y(5+5*l:5+6*l-1),[l,1]);
WC3_hat = reshape(Y(5+6*l:5+7*l-1),[l,1]);
WC4_hat = reshape(Y(5+7*l:5+8*l-1),[l,1]);

Wg1_hat = reshape(Y(5+8*l:5+9*l-1),[l,1]);
Wg2_hat = reshape(Y(5+9*l:5+10*l-1),[l,1]);


WD_hat = [WD1_hat zeros(l,1); zeros(l,1) WD3_hat ; WD2_hat zeros(l,1) ;zeros(l,1) WD4_hat];
WC_hat = [WC1_hat zeros(l,1); zeros(l,1) WC3_hat ; WC2_hat zeros(l,1) ;zeros(l,1) WC4_hat];
Wg_hat = [Wg1_hat zeros(l,1); zeros(l,1) Wg2_hat];

ZD = [Z zeros(l,1); Z zeros(l,1); zeros(l,1) Z;zeros(l,1) Z];
ZC = [Z zeros(l,1); Z zeros(l,1); zeros(l,1) Z;zeros(l,1) Z];
Zg = [Z;Z];


%% control method
e = X-Xd;
e_dot = X_dot - Xd_dot;
s = e_dot + lambda*e;
v = X_dot - lambda*e;
v_dot = Xd_ddot - lambda*e_dot;

tau = J'*(WD_hat'*ZD*v_dot +WC_hat'*ZC*v +Wg_hat'*Zg -Kd*s);% controller

q_ddot = D\(-C*q_dot - G + tau); % closed-loop system

WD_hat_dot = -Gamma_D*ZD*v_dot*s'; % update law
WC_hat_dot = -Gamma_C*ZC*v*s';
Wg_hat_dot = -Gamma_g*Zg*s';

WD_hat_dot = reshape(WD_hat_dot,[8*l,1]);
WC_hat_dot = reshape(WC_hat_dot,[8*l,1]);
Wg_hat_dot = reshape(Wg_hat_dot,[4*l,1]);

%% return data
Yp = [q_dot; q_ddot;WD_hat_dot(1:l);WD_hat_dot(1+2*l:3*l);WD_hat_dot(1+4*l:5*l);WD_hat_dot(1+7*l:8*l);
     WC_hat_dot(1:l);WC_hat_dot(1+2*l:3*l);WC_hat_dot(1+4*l:5*l);WC_hat_dot(1+7*l:8*l);
     Wg_hat_dot(1:l);Wg_hat_dot(1+3*l:4*l)];

end