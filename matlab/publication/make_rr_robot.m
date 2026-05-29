%% make a rigidBodyTree of rr-robot on 2d plane with gravity on y

robot = rigidBodyTree("DataFormat","column",'MaxNumBodies',3);

% using same idea in parameters() under `\source\cal_torqe` directory

addpath(genpath('C:\Users\itsyu\GitHub\working\matlab\publication\sources\torque_cal'));

[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();    % from old script

% link1 assignment with joint1 attachment using MDH convention
body = rigidBody('link1');
joint = rigidBodyJoint('joint1', 'revolute');
setFixedTransform(joint,trvec2tform([0 0 0]));  % mdh transform frame
joint.JointAxis = [0 0 1];
body.Joint = joint;
body.CenterOfMass = [lc1 0 0];
body.Mass = m1;
body.Inertia = [0 0 I1 0 0 0];
cyl = collisionCylinder(0.01, l1);
%But collisionCylinder extends along z-axis by default.
%So you typically need BOTH:
%translation
%rotation
pose = trvec2tform([l1/2 0 0])*axang2tform([0 1 0 pi/2]);
cyl.Pose = pose;
addCollision(body,cyl)
addBody(robot, body, 'base');

% link22 assignment with joint2 attachment
body = rigidBody('link2');
joint = rigidBodyJoint('joint2','revolute');
setFixedTransform(joint, trvec2tform([l1,0,0]));
joint.JointAxis = [0 0 1];
body.Joint = joint;
body.CenterOfMass = [lc2 0 0];
body.Mass = m2;
body.Inertia = [0 0 I2 0 0 0];
cyl = collisionCylinder(0.01, l2);
%But collisionCylinder extends along z-axis by default.
%So you typically need BOTH:
%translation
%rotation
pose = trvec2tform([l2/2 0 0])*axang2tform([0 1 0 pi/2]);
cyl.Pose = pose;
addCollision(body,cyl)
addBody(robot, body, 'link1')

% tool assignment with fixed frame as EOF (no mass, no COM, no I)
body = rigidBody('tool');
joint = rigidBodyJoint('fix1','fixed');
setFixedTransform(joint, trvec2tform([l2, 0, 0]));
body.Joint = joint;
sph = collisionSphere(0.025);
sph.Pose = trvec2tform([0.025/2 0 0]);
addCollision(body,sph)
addBody(robot, body, 'link2');

% configure home position from the initial work space's condition
% X_0 = [1; 0.35] and such that q_0 = [-0.4498;1.5730]

q_0 = inverse_kinematics(l1, l2, [1; 0.35]);  % from old script
%homeConfig = homeConfiguration(robot);

% assign home configuration and gravity

robot.Gravity = [0 -g 0];    % g is in the parameters

showdetails(robot)

%% add collision box for figured output
homeConfig = q_0;
show(robot,homeConfig,"Collisions","on");
axis equal
view(135,25)
grid on
rotate3d on

save my_rr_robot.mat robot homeConfig