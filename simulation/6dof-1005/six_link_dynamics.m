function [M, C, G] = six_link_dynamics(q, dq)
    persistent robot
    if isempty(robot)
        robot = build_robot_model();
    end
    M = massMatrix(robot, q);
    C = velocityProduct(robot, q, dq);
    G = gravityTorque(robot, q);
    
end


function robot = build_robot_model()
    robot = rigidBodyTree('DataFormat','column','MaxNumBodies',6);
    a = [-0.046, 0.22, 0, 0.11, -0.04, 0];
    alpha = [0, pi/2, deg2rad(7.5), deg2rad(-66), pi/2, -pi/2];
    d = [1.1923, 0.34, 0.008, -0.5, 0.63, 0.93];

    for i = 1:6
        body = rigidBody(['body' num2str(i)]);
        joint = rigidBodyJoint(['joint' num2str(i)], 'revolute');
        setFixedTransform(joint, [a(i), alpha(i), d(i), 0], 'dh');
        body.Joint = joint;
        
        
        if i == 1
            addBody(robot, body, 'base');
        else
            addBody(robot, body, ['body' num2str(i-1)]);
        end
    end
end