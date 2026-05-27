% FAT-based force sensorless design without S&L and MRC

clc;clf;clear;

%% initialization
[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();

l = 13; % number of terms of z
T = 17; % period

I = eye(2);
Kd = 300*I;
Kp = 500*I;
Gamma = 1000*eye(l);

Q = diag([4500 4500 500 500]);
A = [zeros(2,2) eye(2);
    -Kp -Kd];
B = [zeros(2,2) eye(2)]';
P = lyap(A',Q); % Lyapunov function

tspan = [0 10]; % simulation time

% initial condition 
X0 = [1 0.35]'; % robot (Cartesian space)
q0 = inverse_kinematics(l1,l2,X0); % transfer into the joint space

W0_reshaped = zeros(1,2*l);
Y0 = [q0(1) q0(2) 0 0 W0_reshaped];

%% ode89

[t,Y] = ode89(@(t,Y) FAT_based_Adaptive_2(t,Y,Kd,Kp,Gamma,B,P,l,T),tspan,Y0);
q = [Y(:,1) Y(:,2)]'; % robot trajectory in the joint space
q_dot = [Y(:,3) Y(:,4)]';

%% recall data
n = length(t);
Xd = zeros(2,n); % desired trajectory in the Cartesian space
X = zeros(2,n);  % robot trajectory in the Cartesian space
qd = zeros(2,n); % desired trajectory in the joint space
tau = zeros(2,n); % control input

for i = 1:n
    [Xd(:,i),Xd_dot,Xd_ddot,qd(:,i),qd_dot,qd_ddot] = desired_trajectory_cartesian(t(i),l1,l2);
    [X(:,i),X_dot] = forward_kenimatics(l1,l2,q(:,i),q_dot(:,i));
    [J,J_dot] = Jacobian_matrix(l1,l2,q(:,i),q_dot(:,i));
    [D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q(:,i),q_dot(:,i),J,J_dot);    
    
    e(:,i) = X(:,i)-Xd(:,i);
    e_dot = X_dot - Xd_dot;
    
    Z = generate_basis(t(i),l,T);
    W_hat = reshape(Y(i,5:5+2*l-1),[l,2]);
    Psi_hat(:,i) = W_hat'*Z;
    tau(:,i) = J'*(Psi_hat(:,i) + Xd_ddot - Kd*e_dot - Kp*e(:,i));
    X_ddot = Dx\(-Cx*X_dot - Gx + J'\tau(:,i));
    Psi(:,i) = Dx*X_ddot + Cx*X_dot + Gx - X_ddot;

end


%% plot figures
figure(1);
plot(X(1,1), X(2,1), 'bx','LineWidth',1.5);
hold;
plot(X(1,:), X(2,:),'k','LineWidth',1.5);
plot(Xd(1,:), Xd(2,:),'r--','LineWidth',1.5);
% grid;
axis equal;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Initial Position','Robot','Desired', 'Interpreter','latex','FontSize', 14);
xlabel('x(m)', 'Interpreter','latex','FontSize', 14);
ylabel('z(m)', 'Interpreter','latex','FontSize', 14);
%title('Tracking Result by FAT-based Design', 'Interpreter','latex','FontSize', 24);
ylim([0.34 0.62]);

figure(2);
subplot(2,1,1);
plot(t, q(1,:),'k','LineWidth',1.5);
hold;
plot(t, qd(1,:),'r--','LineWidth',1.5);

% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Actual','Desired', 'Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$q_1$(rad)', 'Interpreter','latex','FontSize', 14);
%title('First joint angle', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, q(2,:),'k','LineWidth',1.5);
hold;
plot(t, qd(2,:),'r--','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Actual','Desired', 'Interpreter','latex','FontSize', 14,'Location','southeast');
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$q_2$(rad)', 'Interpreter','latex','FontSize', 14);
%title('Second joint angle', 'Interpreter','latex','FontSize', 24);
%% 

figure(3);
subplot(2,1,1);
plot(t, tau(1,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_1$(N$\cdot$m)', 'Interpreter','latex','FontSize', 14);
%title('Control Input $\tau_1$', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, tau(2,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_2$(N$\cdot$ m)', 'Interpreter','latex','FontSize', 14);
%title('Control Input $\tau_2$', 'Interpreter','latex','FontSize', 24);


%% 
figure(4);

subplot(2,1,1);
plot(t, Psi_hat(1,:),'k','LineWidth',1.5);
hold;
plot(t, Psi(1,:),'r--','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\Psi_1$', 'Interpreter','latex','FontSize', 14);
legend('FAT','Actual','FontSize', 14,'Location','southeast')

subplot(2,1,2);
plot(t, Psi_hat(2,:),'k','LineWidth',1.5);
hold;
plot(t, Psi(2,:),'r--','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\Psi_2$', 'Interpreter','latex','FontSize', 14);
legend('FAT','Actual','FontSize', 14)