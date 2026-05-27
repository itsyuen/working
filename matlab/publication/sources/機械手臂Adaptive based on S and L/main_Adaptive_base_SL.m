% FAT-based force sensorless design without S&L and MRC

clc;clf;clear;

%% initialization
[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();
p = [m1*lc1^2+m2*l1^2+I1 m2*lc2^2+I2 m2*l1*lc2 (m1*lc1+m2*l1)*g m2*lc2*g]';
lambda = diag([10 5]);
gamma = 10000*diag([1 1 1 1 1]);
Kd = diag([100,200]);

tspan = [0 10]; % simulation time

% initial condition 
X0 = [1 0.35]'; % robot (Cartesian space)
q0 = inverse_kinematics(l1,l2,X0); % transfer into the joint space
Y0 = [q0(1) q0(2) 0 0 0 0 0 0 0];

%% ode23s

[t,Y] = ode23s(@(t,Y) Adaptive_base_SL(t,Y,lambda,gamma,Kd),tspan,Y0);
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
    s(:,i) = e_dot + lambda*e(:,i);
    v(:,i) = X_dot - lambda*e(:,i);
    v_dot(:,i) = Xd_ddot - lambda*e_dot;

    [Yx] = regressor_matrix(l1,l2,X(:,i),X_dot,v(:,i),v_dot(:,i));
    p_hat = Y(i,5:9)';
    tau(:,i) = J'*(Yx*p_hat -Kd*s(:,i));

    
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
% title('Tracking Result by FAT-based Force Sensorless Design', 'Interpreter','latex','FontSize', 24);
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
%% 

figure(3);
subplot(2,1,1);
plot(t, tau(1,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_1$(N$\cdot$m)', 'Interpreter','latex','FontSize', 14);
% title('Control Input $\tau_1$', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, tau(2,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_2$(N$\cdot$ m)', 'Interpreter','latex','FontSize', 14);
% title('Control Input $\tau_2$', 'Interpreter','latex','FontSize', 24);

%% 

figure(4);
subplot(2,1,1);
plot(t, p(1)*ones(n,1),'r--','LineWidth',1.5);
hold;
plot(t, Y(:,5),'k','LineWidth',1.5);
grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('$p_{1}$','$\hat{p}_{1}$', 'Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$p_{1}$', 'Interpreter','latex','FontSize', 14);
%title('The estimation of $p_{1}$', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, p(2)*ones(n,1),'r--','LineWidth',1.5);
hold;
plot(t, Y(:,6),'k','LineWidth',1.5);
grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('$p_{2}$','$\hat{p}_{2}$', 'Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$p_{2}$', 'Interpreter','latex','FontSize', 14);
%title('The estimation of $p_{2}$', 'Interpreter','latex','FontSize', 24);

%% 

figure(5);
subplot(2,1,1);
plot(t, p(3)*ones(n,1),'r--','LineWidth',1.5);
hold;
plot(t, Y(:,7),'k','LineWidth',1.5);
grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('$p_{3}$','$\hat{p}_{3}$', 'Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$p_{3}$', 'Interpreter','latex','FontSize', 14);
%title('The estimation of $p_{3}$', 'Interpreter','latex','FontSize', 24)
subplot(2,1,2);
plot(t, p(4)*ones(n,1),'r--','LineWidth',1.5);
hold;
plot(t, Y(:,8),'k','LineWidth',1.5);
grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('$p_{4}$','$\hat{p}_{4}$', 'Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$p_{4}$', 'Interpreter','latex','FontSize', 14);
%title('The estimation of $p_{4}$', 'Interpreter','latex','FontSize', 24);
%% 

figure(6);
subplot(2,1,1);
plot(t, p(5)*ones(n,1),'r--','LineWidth',1.5);
hold;
plot(t, Y(:,9),'k','LineWidth',1.5);
grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('$p_{5}$','$\hat{p}_{5}$', 'Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$p_{5}$', 'Interpreter','latex','FontSize', 14);
%title('The estimation of $p_{5}$', 'Interpreter','latex','FontSize', 24);
