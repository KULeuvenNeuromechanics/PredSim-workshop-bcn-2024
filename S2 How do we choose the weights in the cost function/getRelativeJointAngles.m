function relativeJointAngles = getRelativeJointAngles(q1,q2,q3,q4,q5)
%GETRELATIVEJOINTANGLES
%    RELATIVEJOINTANGLES = GETRELATIVEJOINTANGLES(Q1,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    21-Aug-2024 16:54:40

t2 = -q3;
relativeJointAngles = [q1;q1-q2;q2+t2;q4+t2;-q4+q5];
