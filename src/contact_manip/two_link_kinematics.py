#!/usr/bin/env python

import sympy as sp
import numpy as np
import math
import rospy
from sensor_msgs.msg import JointState
from visualization_msgs.msg import Marker
import tf


def hat(a):
    """
    Construct skew-symmetric 3x3 matrix from 3 dim'l vector
    
    :param a: 3 dim'l vector
    :return: 
    """
    return [[0, -a[2], a[1]],
            [a[2], 0, -a[0]],
            [-a[1], a[0], 0]]


def vee(a):
    """
    Extract 6 dim'l twist vector (linear, angular) from 4x4 matrix 
    :param a: 
    :return: 
    """
    twist = np.zeros(6)

    w = a[0:3, 0:3]
    twist[3] = w[1, 2]
    twist[4] = w[0, 2]
    twist[5] = w[1, 0]

    twist[0:3] = a[0:3, 3]

    return twist


def axis_angle_to_rot(axis, angle):
    return sp.eye(3) + sp.sin(angle) * sp.Matrix(hat(axis)) + (1 - sp.cos(angle)) * sp.Matrix(hat(axis))**2


def transformation(twist, theta):
    rot = axis_angle_to_rot(twist[1], theta)
    g = sp.zeros(4, 4)
    g[0:3, 0:3] = rot
    g[0:3, 3] = (sp.eye(3) - rot) * sp.Matrix(np.array(hat(twist[1])).dot(twist[0]).reshape(3, 1)) + \
                sp.Matrix(np.outer(twist[1], twist[1]).dot(twist[0])).reshape(3, 1) * theta
    g[3, 3] = 1
    return g


def adjoint(g):
    rot = g[0:3, 0:3]
    p = g[0:3, 3]

    ad_g = sp.zeros(6, 6)

    ad_g[0:3, 0:3] = rot
    ad_g[0:3, 3:6] = sp.Matrix(hat(p)) * rot
    ad_g[3:6, 3:6] = rot

    return ad_g


class TwoLinkKinematics:
    def __init__(self, l0, l1):
        """
        
        :param l0: length of link 0
        :param l1: length of link 1
        """
        self.l0 = l0
        self.l1 = l1

        self.forward_kinematics = self._forward_kinematics()
        self.spatial_manip_jacobian = self._spatial_manip_jacobian()

    def _forward_kinematics(self):
        theta1, theta2 = sp.symbols('theta1 theta2')

        # First joint
        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, 0])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        # Second joint
        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, self.l0])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        g2 = transformation(twist2, theta2)

        g0 = sp.eye(4)
        g0[2, 3] = self.l0 + self.l1

        g_st = g1 * g2 * g0

        def eval(th1, th2):
            # print("th1 = {}, th2 = {}".format(th1, th2))
            return g_st.evalf(subs={'theta1': th1, 'theta2': th2})

        return eval

    def _spatial_manip_jacobian(self):
        theta1, theta2 = sp.symbols('theta1 theta2')

        jac = sp.zeros(6, 2)

        # First twist (first column of manipulator Jacobian)
        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, 0])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        # Second twist (second column of manipulator Jacobian)
        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, self.l0])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)

        jac[:, 0] = sp.Matrix(np.concatenate(twist1).reshape(6, 1))
        jac[:, 1] = adjoint(g1) * sp.Matrix(np.concatenate(twist2).reshape(6, 1))

        def eval(th1, th2):
            return jac.evalf(subs={'theta1': th1, 'theta2': th2})

        return eval


class JointStatesHolder:
    def __init__(self):
        self.joints_msg = JointState()
        self.joints_msg.name = ['joint0', 'joint1']
        self.joints_msg.position = [0.0, 0.0]

    def callback(self, msg):
        self.joints_msg = msg


def quaternion_to_axis_angle(q):
    angle = 2.0 * math.acos(q[3])
    denom = math.sqrt(1.0 - math.pow(q[3], 2))
    axis = (q[0] / denom,
            q[1] / denom,
            q[2] / denom)
    return axis, angle


def main():
    kin = TwoLinkKinematics(l0=0.4, l1=0.4)
    holder = JointStatesHolder()

    rospy.init_node('two_link_kinematics')
    rate = rospy.Rate(10)

    joint_state_sub = rospy.Subscriber('/joint_states', JointState, holder.callback)
    marker_pub = rospy.Publisher('/visualization_marker', Marker, queue_size=10)
    joint_state_pub = rospy.Publisher('/joint_states', JointState, queue_size=5)

    p_goal = (0, -0.496276260471099, 0.554429385958024)
    q_goal = (0.525309397897, 0.0, 0.0, 0.850911297658)

    while not rospy.is_shutdown():
        g_st = kin.forward_kinematics(holder.joints_msg.position[0],
                                      holder.joints_msg.position[1])
        print(g_st)

        jac_st = kin.spatial_manip_jacobian(holder.joints_msg.position[0],
                                            holder.joints_msg.position[1])
        print(jac_st)

        tool_marker = Marker()
        tool_marker.header.frame_id = 'root'
        tool_marker.ns = 'tool_frame'
        tool_marker.id = 0
        tool_marker.type = Marker.SPHERE
        tool_marker.pose.position.x = g_st[0, 3]
        tool_marker.pose.position.y = g_st[1, 3]
        tool_marker.pose.position.z = g_st[2, 3]
        q = tf.transformations.quaternion_from_matrix(np.array(g_st[0:4, 0:4]))
        tool_marker.pose.orientation.x = q[0]
        tool_marker.pose.orientation.y = q[1]
        tool_marker.pose.orientation.z = q[2]
        tool_marker.pose.orientation.w = q[3]
        tool_marker.scale.x = 0.15
        tool_marker.scale.y = 0.3
        tool_marker.scale.z = 0.15
        tool_marker.color.r = 0.5
        tool_marker.color.g = 0.0
        tool_marker.color.b = 0.0
        tool_marker.color.a = 0.35
        marker_pub.publish(tool_marker)

        goal_marker = Marker()
        goal_marker.header.frame_id = 'root'
        goal_marker.ns = 'goal_frame'
        goal_marker.id = 0
        goal_marker.type = Marker.SPHERE
        goal_marker.pose.position.x = p_goal[0]
        goal_marker.pose.position.y = p_goal[1]
        goal_marker.pose.position.z = p_goal[2]
        goal_marker.pose.orientation.x = q_goal[0]
        goal_marker.pose.orientation.y = q_goal[1]
        goal_marker.pose.orientation.z = q_goal[2]
        goal_marker.pose.orientation.w = q_goal[3]
        goal_marker.scale.x = 0.15
        goal_marker.scale.y = 0.3
        goal_marker.scale.z = 0.15
        goal_marker.color.r = 0.0
        goal_marker.color.g = 0.5
        goal_marker.color.b = 0.0
        goal_marker.color.a = 0.35
        marker_pub.publish(goal_marker)

        print("Current pos: ({}, {}, {}), orien: ({}, {}, {}, {})".format(tool_marker.pose.position.x,
                                                                          tool_marker.pose.position.y,
                                                                          tool_marker.pose.position.z,
                                                                          tool_marker.pose.orientation.x,
                                                                          tool_marker.pose.orientation.y,
                                                                          tool_marker.pose.orientation.z,
                                                                          tool_marker.pose.orientation.w))

        p_diff = (p_goal[0] - g_st[0, 3],
                  p_goal[1] - g_st[1, 3],
                  p_goal[2] - g_st[2, 3])
        q_diff = tf.transformations.quaternion_multiply(tf.transformations.quaternion_inverse(q),
                                                        q_goal)
        axis_diff, angle_diff = quaternion_to_axis_angle(q_diff)
        print("Diff axis: ({}, {}, {}), diff angle: {}".format(axis_diff[0],
                                                               axis_diff[1],
                                                               axis_diff[2],
                                                               angle_diff))
        print("Pos diff: ({}, {}, {})".format(p_diff[0],
                                              p_diff[1],
                                              p_diff[2]))

        twist = np.array([p_diff[0], p_diff[1], p_diff[2],
                          axis_diff[0] * angle_diff, axis_diff[1] * angle_diff, axis_diff[2] * angle_diff]) / 4.0
        jac = np.array(jac_st.tolist())
        joint_vel = np.linalg.pinv(jac.astype(np.float32)).dot(twist)
        print("joint vel = {}".format(joint_vel))

        joints_msg = JointState()
        joints_msg.header.stamp = rospy.Time.now()
        joints_msg.name = ['joint0', 'joint1']
        joints_msg.position = [holder.joints_msg.position[0] + 0.1 * joint_vel[0],
                               holder.joints_msg.position[1] + 0.1 * joint_vel[1]]
        joint_state_pub.publish(joints_msg)

        rate.sleep()


if __name__ == '__main__':
    main()
