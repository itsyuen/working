% FAT-based force sensorless design without S&L and MRC

clc;clf;clear;

%% initialization
[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();

tspan = [0 10]; % simulation time

% initial condition 
X0 = [1 0.35]'; % robot (Cartesian space)
q0 = inverse_kinematics(l1,l2,X0); % transfer into the joint space
Y0 = [q0(1) q0(2) 0 0];

%% ode89

[t,Y] = ode89(@(t,Y) Computed_Torque_Method_1(t,Y),tspan,Y0);
q = [Y(:,1) Y(:,2)]'; % robot trajectory in the joint space
q_dot = [Y(:,3) Y(:,4)]';

%% recall data
n = length(t);
Xd = zeros(2,n); % desired trajectory in the Cartesian space
X = zeros(2,n);  % robot trajectory in the Cartesian space
qd = zeros(2,n); % desired trajectory in the joint space
tau = zeros(2,n); % control input
taux = zeros(2,n); % control input
q_ddot = zeros(2,n);

for i = 1:n
    [Xd(:,i),Xd_dot,Xd_ddot] = desired_trajectory_cartesian(t(i),l1,l2);
    [X(:,i),X_dot] = forward_kenimatics(l1,l2,q(:,i),q_dot(:,i));
    qd(:,i) = inverse_kinematics(l1,l2,Xd(:,i));
    [J,J_dot] = Jacobian_matrix(l1,l2,q(:,i),q_dot(:,i));
    [D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q(:,i),q_dot(:,i),J,J_dot);

    e = X(:,i) - Xd(:,i);
    e_dot = X_dot - Xd_dot;
    
    Kv = diag([100, 100]);
    Kp = diag([1000, 1000]);

    tau(:,i) = J'*(Cx*X_dot+Gx+Dx*(Xd_ddot-Kv*(X_dot-Xd_dot)-Kp*(X(:,i)-Xd(:,i))));
    q_ddot(:,i) = D\(-C*q_dot(:,i) - G + tau(:,i));
    X_ddot = J_dot*q_dot(:,i)+J*q_ddot(:,i);

    [Y,p] = regressor_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q(:,i),q_dot(:,i),q_ddot(:,i));
    tauy(:,i) = Y*p;
    
    taux(:,i) = J'*(Dx*X_ddot+Cx*X_dot+Gx);

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

ylim([0.34 0.62]);
% title('Tracking Result by FAT-based Force Sensorless Design', 'Interpreter','latex','FontSize', 24);

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
% title('First joint angle', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, q(2,:),'k','LineWidth',1.5);
hold;
plot(t, qd(2,:),'r--','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Actual','Desired', 'Interpreter','latex','FontSize', 14,'Location','southeast');
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$q_2$(rad)', 'Interpreter','latex','FontSize', 14);
% title('Second joint angle', 'Interpreter','latex','FontSize', 24);

figure(3);
subplot(2,1,1);
plot(t, taux(1,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_1$(N$\cdot$m)', 'Interpreter','latex','FontSize', 14);
% title('Control Input $\tau_1$', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, taux(2,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_2$(N$\cdot$ m)', 'Interpreter','latex','FontSize', 14);
% title('Control Input $\tau_2$', 'Interpreter','latex','FontSize', 24);


