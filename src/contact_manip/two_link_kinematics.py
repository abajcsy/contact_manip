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


def adjoint_inv(g):
    rot = g[0:3, 0:3]
    p = g[0:3, 3]

    ad_inv_g = sp.zeros(6, 6)

    ad_inv_g[0:3, 0:3] = rot.transpose()
    ad_inv_g[0:3, 3:6] = -rot.transpose() * sp.Matrix(hat(p))
    ad_inv_g[3:6, 3:6] = rot.transpose()

    return ad_inv_g


class TwoLinkKinematics:
    def __init__(self, l1, l2):
        """
        
        :param l1: length of link 1
        :param l2: length of link 2
        """
        self.l1 = l1
        self.l2 = l2

        self.forward_kinematics = self._forward_kinematics()
        self.spatial_manip_jacobian = self._spatial_manip_jacobian()

        # self.joint_accs = self._joint_accs()
        # self._dynamics_sym()
        # self._forward_kinematics_sym()
        # self._spatial_manip_jacobian_sym()
        self._body_manip_jacobian_sym()

    def _dynamics_sym(self):
        theta1, theta2, \
        dtheta1, dtheta2, \
        g, m1, m2, \
        i1x, i1y, i1z, \
        i2x, i2y, i2z, \
        l1, l2 = sp.symbols('\\theta_1 \\theta_2 \\dot{\\theta}_1 \\dot{\\theta}_2 g m_1 m_2 I_{1x} I_{1y} I_{1z} I_{2x} I_{2y} I_{2z} l_1 l_2')

        z_eps = 0.0  # 1e-2
        z_offset = l1 + l2 + z_eps

        # Link 1 inertia matrix
        M1 = sp.zeros(6, 6)
        M1[0:3, 0:3] = sp.eye(3) * m1
        M1[3, 3] = i1x
        M1[4, 4] = i1y
        M1[5, 5] = i1z

        # Link 2 inertia matrix
        M2 = sp.zeros(6, 6)
        M2[0:3, 0:3] = sp.eye(3) * m2
        M2[3, 3] = i2x
        M2[4, 4] = i2y
        M2[5, 5] = i2z

        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, z_offset])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, z_offset + l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        g2 = transformation(twist2, theta2)

        # Forward kinematics for position of link 1
        g_s1_0 = sp.eye(4)
        g_s1_0[2, 3] = z_offset + 0.5 * l1

        # Forward kinematics for position of link 2
        g_s2_0 = sp.eye(4)
        g_s2_0[2, 3] = z_offset + l1 + 0.5 * l2

        jac_1 = sp.zeros(6, 2)
        jac_1[:, 0] = sp.simplify(adjoint_inv(g1 * g_s1_0) * sp.Matrix(np.concatenate(twist1).reshape(6, 1)))

        jac_2 = sp.zeros(6, 2)
        jac_2[:, 0] = sp.simplify(adjoint_inv(g1 * g2 * g_s2_0) * sp.Matrix(np.concatenate(twist1).reshape(6, 1)))
        jac_2[:, 1] = sp.simplify(adjoint_inv(g2 * g_s2_0) * sp.Matrix(np.concatenate(twist2).reshape(6, 1)))

        M = sp.simplify(jac_1.transpose() * M1 * jac_1 + jac_2.transpose() * M2 * jac_2)

        # Coriolis matrix
        n = 2
        thetas = [theta1, theta2]
        dthetas = [dtheta1, dtheta2]
        C = sp.zeros(n, n)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    C[i, j] += (sp.diff(M[i, j], thetas[k]) +
                                sp.diff(M[i, k], thetas[j]) -
                                sp.diff(M[k, k], thetas[i])) * dthetas[k]

                C[i, j] *= 0.5

        C = sp.simplify(C)

        # Potential forces
        g_s1 = g1 * g_s1_0
        g_s2 = g1 * g2 * g_s2_0

        V = g * (m1 * g_s1[2, 3] + m2 * g_s2[2, 3])
        N = sp.zeros(2, 1)
        N[0, 0] = sp.diff(V, theta1)
        N[1, 0] = sp.diff(V, theta2)

        N = sp.simplify(N)

        print("M = {}".format(sp.latex(M)))
        print("C = {}".format(sp.latex(C)))
        print("N = {}".format(sp.latex(N)))

    def _joint_accs(self):
        g = 9.8

        z_offset = self.l1 + self.l2 + 0.01

        theta1, theta2, dtheta1, dtheta2 = sp.symbols('theta1 theta2 dtheta1 dtheta2')

        # Mass of link 1
        m1 = 1.0
        # Moment of inertia for link 1 (modeled as a rod)
        i1 = (m1 * self.l1**2) / 12.0

        # Mass of link 2
        m2 = 1.0
        # Moment of inertia for link 2 (modeled as a rod)
        i2 = (m2 * self.l2**2) / 12.0

        # Link 1 inertia matrix
        M1 = sp.zeros(6, 6)
        M1[0:3, 0:3] = sp.eye(3) * m1
        M1[3, 3] = i1
        M1[4, 4] = i1
        M1[5, 5] = 0.0

        # Link 2 inertia matrix
        M2 = sp.zeros(6, 6)
        M2[0:3, 0:3] = sp.eye(3) * m2
        M2[3, 3] = i2
        M2[4, 4] = i2
        M2[5, 5] = 0.0

        print("M1 = {}, M2 = {}".format(M1, M2))

        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, z_offset])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, z_offset + self.l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        g2 = transformation(twist2, theta2)

        # Forward kinematics for position of link 1
        g_s1_0 = sp.eye(4)
        g_s1_0[2, 3] = z_offset + 0.5 * self.l1

        # Forward kinematics for position of link 2
        g_s2_0 = sp.eye(4)
        g_s2_0[2, 3] = z_offset + self.l1 + 0.5 * self.l2

        jac_1 = sp.zeros(6, 2)
        # TODO change to (g1 * g_s1_0) (doesn't make sense to have dependency on theta2)
        # jac_1[:, 0] = adjoint_inv(g1 * g2 * g_s1_0) * sp.Matrix(np.concatenate(twist1).reshape(6, 1))
        print(np.concatenate(twist1))
        jac_1[:, 0] = sp.simplify(adjoint_inv(g1 * g_s1_0) * sp.Matrix(np.concatenate(twist1).reshape(6, 1)))

        jac_2 = sp.zeros(6, 2)
        jac_2[:, 0] = sp.simplify(adjoint_inv(g1 * g2 * g_s2_0) * sp.Matrix(np.concatenate(twist1).reshape(6, 1)))
        jac_2[:, 1] = sp.simplify(adjoint_inv(g2 * g_s2_0) * sp.Matrix(np.concatenate(twist2).reshape(6, 1)))

        M = sp.simplify(jac_1.transpose() * M1 * jac_1 + jac_2.transpose() * M2 * jac_2)
        # print("M = {}".format(M.shape))
        print("jac_1 = {}".format(jac_1))
        print("jac_2 = {}".format(jac_2))
        print("M = {}".format(M))
        # print("jac_1 = {}".format(jac_1.evalf(subs={'theta1': 2.0, 'theta2': 2.0})))
        # print("jac_1 = {}\njac_2 = {}".format(jac_1, jac_2))
        # print("M = {}".format(M))

        # Coriolis matrix
        n = 2
        thetas = [theta1, theta2]
        dthetas = [dtheta1, dtheta2]
        C = sp.zeros(n, n)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    C[i, j] += (sp.diff(M[i, j], thetas[k]) +
                                sp.diff(M[i, k], thetas[j]) -
                                sp.diff(M[k, k], thetas[i])) * dthetas[k]

                C[i, j] *= 0.5

        C = sp.simplify(C)
        print("C = {}".format(C))

        # Potential forces
        g_s1 = g1 * g_s1_0
        g_s2 = g1 * g2 * g_s2_0

        V = g * (m1 * g_s1[2, 3] + m2 * g_s2[2, 3])
        N = sp.zeros(2, 1)
        N[0, 0] = sp.diff(V, theta1)
        N[1, 0] = sp.diff(V, theta2)

        print("N = {}".format(N))

        th = sp.zeros(2, 1)
        th[0, 0] = theta1
        th[1, 0] = theta2

        # TODO add in applied torques
        ddtheta = M.inv() * (-N - C * th)

        print(ddtheta)
        print("done1")

        def eval(th1, th2, dth1, dth2):
            return ddtheta.evalf(subs={'theta1': th1,
                                       'theta2': th2,
                                       'dtheta1': dth1,
                                       'dtheta2': dth2})

        print("done2")

        return eval

    def _forward_kinematics(self):
        theta1, theta2 = sp.symbols('theta1 theta2')

        # First joint
        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, 0])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        # Second joint
        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, self.l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        g2 = transformation(twist2, theta2)

        g0 = sp.eye(4)
        g0[2, 3] = self.l1 + self.l2

        g_st = g1 * g2 * g0

        def eval(th1, th2):
            # print("th1 = {}, th2 = {}".format(th1, th2))
            return g_st.evalf(subs={'theta1': th1, 'theta2': th2})

        return eval

    def _forward_kinematics_sym(self):
        theta1, theta2, l1, l2 = sp.symbols('\\theta_1 \\theta_2 l_1 l_2')

        # First joint
        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, 0])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        g0 = sp.eye(4)
        g0[2, 3] = l1

        print("g_sl1 = {}".format(sp.latex(sp.simplify(g1 * g0))))

        # Second joint
        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        g2 = transformation(twist2, theta2)

        g0 = sp.eye(4)
        g0[2, 3] = l1 + l2

        print("g_st = {}".format(sp.latex(sp.simplify(g1 * g2 * g0))))

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
        q2 = np.array([0, 0, self.l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)

        jac[:, 0] = sp.Matrix(np.concatenate(twist1).reshape(6, 1))
        jac[:, 1] = adjoint(g1) * sp.Matrix(np.concatenate(twist2).reshape(6, 1))

        print("jac = {}".format(jac))

        def eval(th1, th2):
            return jac.evalf(subs={'theta1': th1, 'theta2': th2})

        return eval

    def _spatial_manip_jacobian_sym(self):
        theta1, theta2, l1, l2 = sp.symbols('\\theta_1 \\theta_2 l_1 l_2')

        z_eps = 0.0  # 1e-2
        z_offset = l1 + l2 + z_eps

        jac = sp.zeros(6, 2)

        # First twist (first column of manipulator Jacobian)
        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, z_offset])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        # Second twist (second column of manipulator Jacobian)
        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, z_offset + l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)

        jac[:, 0] = sp.Matrix(np.concatenate(twist1).reshape(6, 1))
        jac[:, 1] = adjoint(g1) * sp.Matrix(np.concatenate(twist2).reshape(6, 1))

        g_st_0 = sp.eye(4)
        g_st_0[2, 3] = z_offset + l1 + l2

        print("J^s = {}".format(sp.latex(sp.simplify(jac))))
        print("J^b = {}".format(sp.latex(sp.simplify(adjoint_inv(g_st_0) * jac))))

    def _body_manip_jacobian_sym(self):
        theta1, theta2, l1, l2 = sp.symbols('\\theta_1 \\theta_2 l_1 l_2')

        z_eps = 0.0  # 1e-2
        # z_offset = l1 + l2 + z_eps
        z_offset = 0

        w1 = np.array([1, 0, 0])
        q1 = np.array([0, 0, z_offset])
        twist1 = (-np.array(hat(w1)).dot(q1), w1)
        g1 = transformation(twist1, theta1)

        w2 = np.array([1, 0, 0])
        q2 = np.array([0, 0, z_offset + l1])
        twist2 = (-np.array(hat(w2)).dot(q2), w2)
        g2 = transformation(twist2, theta2)

        # Forward kinematics for initial pose of link 2
        g_s2_0 = sp.eye(4)
        g_s2_0[2, 3] = z_offset + l1 + l2
        print(g1)
        print(g2)
        print(g_s2_0)

        jac = sp.zeros(6, 2)
        jac[:, 0] = sp.simplify(adjoint(g1 * g2 * g_s2_0).inv() * sp.Matrix(np.concatenate(twist1).reshape(6, 1)))
        jac[:, 1] = sp.simplify(adjoint(g2 * g_s2_0).inv() * sp.Matrix(np.concatenate(twist2).reshape(6, 1)))
        # print(jac)

        # TODO is the above correct??
        print("J^b = {}".format(sp.latex(sp.simplify(jac))))


class JointStatesHolder:
    def __init__(self, j0=0.0, j1=0.0):
        self.joints_msg = JointState()
        self.joints_msg.name = ['joint0', 'joint1']
        self.joints_msg.position = [j0, j1]

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
    kin = TwoLinkKinematics(l1=0.4, l2=0.4)
    exit(0)
    holder = JointStatesHolder(j1=0.01)

    rospy.init_node('two_link_kinematics')
    rate = rospy.Rate(10)
    # dt = 1e-3
    dt = 1e-2

    joint_state_sub = rospy.Subscriber('/joint_states', JointState, holder.callback)
    marker_pub = rospy.Publisher('/visualization_marker', Marker, queue_size=10)
    joint_state_pub = rospy.Publisher('/joint_states', JointState, queue_size=5)

    p_goal = (0, -0.496276260471099, 0.554429385958024)
    q_goal = (0.525309397897, 0.0, 0.0, 0.850911297658)

    joint_vel = np.zeros(2)

    print("Starting...")

    while not rospy.is_shutdown():
        g_st = kin.forward_kinematics(holder.joints_msg.position[0],
                                      holder.joints_msg.position[1])
        print(g_st)

        jac_st = kin.spatial_manip_jacobian(holder.joints_msg.position[0],
                                            holder.joints_msg.position[1])
        print(jac_st)

        joint_accs = kin.joint_accs(holder.joints_msg.position[0],
                                    holder.joints_msg.position[1],
                                    joint_vel[0],
                                    joint_vel[1])
        # joint_accs = np.zeros(2)
        print("joint accs = {}".format(joint_accs))
        joint_vel[0] += dt * joint_accs[0]
        joint_vel[1] += dt * joint_accs[1]

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

        # print("Current pos: ({}, {}, {}), orien: ({}, {}, {}, {})".format(tool_marker.pose.position.x,
        #                                                                   tool_marker.pose.position.y,
        #                                                                   tool_marker.pose.position.z,
        #                                                                   tool_marker.pose.orientation.x,
        #                                                                   tool_marker.pose.orientation.y,
        #                                                                   tool_marker.pose.orientation.z,
        #                                                                   tool_marker.pose.orientation.w))
        #
        # p_diff = (p_goal[0] - g_st[0, 3],
        #           p_goal[1] - g_st[1, 3],
        #           p_goal[2] - g_st[2, 3])
        # q_diff = tf.transformations.quaternion_multiply(tf.transformations.quaternion_inverse(q),
        #                                                 q_goal)
        # axis_diff, angle_diff = quaternion_to_axis_angle(q_diff)
        # print("Diff axis: ({}, {}, {}), diff angle: {}".format(axis_diff[0],
        #                                                        axis_diff[1],
        #                                                        axis_diff[2],
        #                                                        angle_diff))
        # print("Pos diff: ({}, {}, {})".format(p_diff[0],
        #                                       p_diff[1],
        #                                       p_diff[2]))
        #
        # twist = np.array([p_diff[0], p_diff[1], p_diff[2],
        #                   axis_diff[0] * angle_diff, axis_diff[1] * angle_diff, axis_diff[2] * angle_diff]) / 4.0
        # jac = np.array(jac_st.tolist())
        # joint_vel = np.linalg.pinv(jac.astype(np.float32)).dot(twist)
        # print("joint vel = {}".format(joint_vel))

        joints_msg = JointState()
        joints_msg.header.stamp = rospy.Time.now()
        joints_msg.name = ['joint0', 'joint1']
        joints_msg.position = [holder.joints_msg.position[0] + 0.1 * joint_vel[0],
                               holder.joints_msg.position[1] + 0.1 * joint_vel[1]]
        joint_state_pub.publish(joints_msg)

        rate.sleep()


if __name__ == '__main__':
    main()
