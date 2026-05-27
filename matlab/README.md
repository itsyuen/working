# working matlab's repo

This repo is dedicated to working projects that will provide matlab/simulink code that is ready to publish or code-gen for product development.

## repo directory

1st directory dedicated to adaptive controller deployment on RST (Robotics System Toolbox) for MATLAB/simulink environment simulation

## publication goal

This journal publication is trageted towards a parallel highlevel dynamcis adaptive controller running in simulink, to drive a full modeled standard ROS type rigidBodyTree object with full rigidBody and rigidBodyJoint, under the RST setup and simulated fully with simulink.

The matrix and all the state estimation will be included inside the adaptive folder from `learning`.

### robot setup

1. using modified d-h table for joint and body frame that is natrual to rigid body tree object. i.e. $f_{0} \sim f_{3}$

2. notice that $^0f_1$ defines the frame of body 1 (using the same joint 1 frame), with respect to world frame $\{0\}$. And it can be always set to $^1f_1 = eye(4)$.

3. 