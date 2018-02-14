#!/usr/bin/env python

import sympy as sp
import numpy as np
import rospy
from sensor_msgs.msg import JointState
from visualization_msgs.msg import Marker
import tf


def hat(a):
    """
    Construct skew-symmetric matrix from 3 dim'l vector
    
    :param a: 3 dim'l vector
    :return: 
    """
    return [[0, -a[2], a[1]],
            [a[2], 0, -a[0]],
            [-a[1], a[0], 0]]


def axis_angle_to_rot(axis, angle):
    return sp.eye(3) + sp.sin(angle) * sp.Matrix(hat(axis)) + (1 - sp.cos(angle)) * sp.Matrix(hat(axis))**2


class TwoLinkKinematics:
    def __init__(self, l0, l1):
        """
        
        :param l0: length of link 0
        :param l1: length of link 1
        """
        self.l0 = l0
        self.l1 = l1

        self.forward_kinematics = self._forward_kinematics()

    def _forward_kinematics(self):
        theta1, theta2 = sp.symbols('theta1 theta2')

        # First joint
        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, 0])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        rot1 = axis_angle_to_rot(twist1[1], theta1)
        tmp1 = (sp.eye(3) - rot1) * sp.Matrix(np.array(hat(twist1[1])).dot(twist1[0]).reshape(3, 1)) + \
               sp.Matrix(np.outer(twist1[1], twist1[1]).dot(twist1[0])).reshape(3, 1) * theta1
        g1 = sp.zeros(4, 4)
        g1[0:3, 0:3] = rot1
        g1[0:3, 3] = tmp1
        g1[3, 3] = 1
        # print(g1)

        # Second joint
        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, self.l0])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        rot2 = axis_angle_to_rot(twist2[1], theta2)
        tmp2 = (sp.eye(3) - rot2) * sp.Matrix(np.array(hat(twist2[1])).dot(twist2[0]).reshape(3, 1)) + \
               sp.Matrix(np.outer(twist2[1], twist2[1]).dot(twist2[0])).reshape(3, 1) * theta2
        g2 = sp.zeros(4, 4)
        g2[0:3, 0:3] = rot2
        g2[0:3, 3] = tmp2
        g2[3, 3] = 1
        # print(g2)

        g0 = sp.eye(4)
        g0[2, 3] = self.l0 + self.l1

        g_st = g1 * g2 * g0

        def eval(th1, th2):
            print("th1 = {}, th2 = {}".format(th1, th2))
            return g_st.evalf(subs={'theta1': th1, 'theta2': th2})

        return eval


class JointStatesHolder:
    def __init__(self):
        self.joints_msg = JointState()
        self.joints_msg.name = ['joint0', 'joint1']
        self.joints_msg.position = [0.0, 0.0]

    def callback(self, msg):
        self.joints_msg = msg


def main():
    kin = TwoLinkKinematics(l0=0.4, l1=0.4)
    holder = JointStatesHolder()

    rospy.init_node('two_link_kinematics')
    rate = rospy.Rate(10)

    joint_state_sub = rospy.Subscriber('/joint_states', JointState, holder.callback)
    marker_pub = rospy.Publisher('/visualization_marker', Marker, queue_size=10)

    while not rospy.is_shutdown():
        g_st = kin.forward_kinematics(holder.joints_msg.position[0],
                                      holder.joints_msg.position[1])
        print(g_st)

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
        tool_marker.scale.x = 0.25
        tool_marker.scale.y = 0.25
        tool_marker.scale.z = 0.25
        tool_marker.color.r = 0.5
        tool_marker.color.g = 0.0
        tool_marker.color.b = 0.0
        tool_marker.color.a = 0.35
        marker_pub.publish(tool_marker)

        rate.sleep()


if __name__ == '__main__':
    main()
