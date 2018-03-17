#! /usr/bin/env python
import numpy as np
import rospy
from sensor_msgs.msg import JointState
from visualization_msgs.msg import Marker
import tf
import sys, select, os
import sensor_msgs.msg

class TwoLinkArm:
	"""
	Base information for simulation of a two-link arm.
	"""

	def __init__(self, init_q=[.75613, 1.8553], init_dq=[0.0,0.0], l1=0.4, l2=0.4, dt=1e-5):
		self.DOF = 2
		self.init_q = np.zeros(self.DOF) if init_q is None else init_q
		self.init_dq = np.zeros(self.DOF) if init_dq is None else init_dq
		self.l1 = l1
		self.l2 = l2
		self.L = np.array([self.l1, self.l2])
		# mass of links
		self.m1 = 1.98
		self.m2 = 1.32
		# z axis inertia moment of links
		izz1 = 15.0; izz2 = 8.0
		# create mass matrices at COM for each link
		self.M1 = np.zeros((6,6))
		self.M2 = np.zeros((6,6))
		self.M1[0:3,0:3] = np.eye(3)*self.m1
		self.M1[3:,3:] = np.eye(3)*izz1
		self.M2[0:3,0:3] = np.eye(3)*self.m2
		self.M2[3:,3:] = np.eye(3)*izz2

		self.rest_angles = np.array([np.pi/4.0, np.pi/4.0])

		if dt is not None:
			self.dt = dt

		# compute non changing constants
		self.K1 = (self.m1 + self.m2) * self.l1**2.0 + self.m2 * self.l2**2.0
		self.K2 = 2.0 * self.m2 * self.l1 * self.l2
		self.K3 = self.m2 * self.l2**2.0
		self.K4 = self.m2 * self.l1 * self.l2

		# force at end-effecor
		self.fEE = np.array([0.0, 0.0])

		# setup initial configurations
		self.q = self.init_q
		self.dq = self.init_dq
		self.t = 0.0

		# make temporary default message
		self.joints_msg = None

	def gen_jacCOM1(self, q=None):
		"""
		Generates the Jacobian from the COM of the first
		link to the origin frame
		"""
		q = self.q if q is None else q

		JCOM1 = np.zeros((6,2))
		#JCOM1[0,0] = self.l1 / 2. * -np.sin(q[0])
		#JCOM1[1,0] = self.l1 / 2. * np.cos(q[0])
		JCOM1[0,0] = self.l1 / 2. * -np.sin(q[0])
		JCOM1[1,0] = self.l1 / 2. * np.cos(q[0])
		JCOM1[5,0] = 1.0

		return JCOM1

	def gen_jacCOM2(self, q=None):
		"""
		Generates the Jacobian from the COM of the second
		link to the origin frame
		"""
		q = self.q if q is None else q

		JCOM2 = np.zeros((6,2))
		# define column entries right to left
		#JCOM2[0,1] = self.l2 / 2. * -np.sin(q[0]+q[1])
		#JCOM2[1,1] = self.l2 / 2. * np.cos(q[0]+q[1])
		JCOM2[0,1] = self.l2 / 2. * -np.sin(q[0]+q[1])
		JCOM2[1,1] = self.l2 / 2. * np.cos(q[0]+q[1])
		JCOM2[5,1] = 1.0

		JCOM2[0,0] = self.l1 * -np.sin(q[0]) + JCOM2[0,1]
		JCOM2[1,0] = self.l1 * np.cos(q[0]) + JCOM2[1,1]
		JCOM2[5,0] = 1.0

		return JCOM2

	def gen_jacEE(self, q=None):
		"""
		Generates the Jacobian from end-effector to the origin frame
		"""
		q = self.q if q is None else q

		JEE = np.zeros((2,2))
		# define column entries right to left
		JEE[0,1] = self.l2 * -np.sin(q[0]+q[1])
		JEE[1,1] = self.l2 * np.cos(q[0]+q[1])

		JEE[0,0] = self.l1 * -np.sin(q[0]) + JEE[0,1]
		JEE[1,0] = self.l1 * np.cos(q[0]) + JEE[1,1]

		return JEE

	def gen_Mq(self, q=None):
		"""
		Generates the mass matrix for the arm in joint space
		"""
		# get the instantaneous Jacobians
		JCOM1 = self.gen_jacCOM1(q=q)
		JCOM2 = self.gen_jacCOM2(q=q)
		# generate the mass matrix in joint space
		Mq = np.dot(JCOM1.T, np.dot(self.M1, JCOM1)) + \
			np.dot(JCOM2.T, np.dot(self.M2, JCOM2))

		return Mq

	def gen_Gq(self, q=None):
		"""
		Generates the gravity matrix for the arm in joint space
		"""
		q = self.q if q is None else q

		C1 = np.cos(self.q[0])
		C12 = np.cos(self.q[0] + self.q[1])
		g = -9.8
		G1 = (self.m1 + self.m2)*g*self.l1*C1 + self.m2*g*self.l2*C12
		G2 = self.m2*g*self.l2*C12

		return np.array([G1, G2])

	def gen_Fq(self, forceEE=None):
		"""
		Generates the forces at each joint from a force at the end-effector.
		Fq = jacEE.T(q)*Fx
		"""
		Fx = self.fEE if forceEE is None else forceEE

		JEE = self.gen_jacEE()
		Fq = np.dot(JEE.T,Fx)

		return Fq

	def fwd_kinematics(self, q=None):
		"""
		Compute (x,y) position of the hand given a configuration

		q np.array: a set of angles to return positions for

		Pos_Elbow = [ l1*cos(q0); <--(x1)
		l1*sin(q0); <--(y1)
		0 ]

		Pos_EE = [ l1*cos(q0) + l2*cos(q0+q1); <--(x2)
		l1*sin(q0) + l2*sin(q0+q1); <--(y2)
		0 ]
		"""
		q = self.q if q is None else q

		x = np.cumsum([0,
			self.l1 * np.cos(q[0]),
			self.l2 * np.cos(q[0]+q[1])])
		y = np.cumsum([0,
			self.l1 * np.sin(q[0]),
			self.l2 * np.sin(q[0]+q[1])])
		return np.array([x, y])

	def inv_kinematics(self, xy):
		"""
		Calculate the joint angles for a given (x,y) hand position
		"""
		import scipy.optimize
		# function to optimize
		def distance_to_target(q, xy, L):
			x = L[0] * np.cos(q[0]) + L[1] * np.cos(q[0] + q[1])
			y = L[0] * np.sin(q[0]) + L[1] * np.sin(q[0] + q[1])
			return np.sqrt((x - xy[0])**2 + (y - xy[1])**2)

		return scipy.optimize.minimize(fun=distance_to_target, x0=self.q,
			args=([xy[0], xy[1]], self.L))['x']

	def apply_torque(self, u, dt=None):
		"""
		Takes in a torque and time step and updates the
		arm simulation accordingly.

		Solves for ddq, dq, q given u: 
		u = M(q)*ddq + V(q,dq) + G(q)

		u np.array: the control signal to apply
		dt float: the time step
		"""
		if dt is None:
			dt = self.dt

		# equations solved for angles
		C1 = np.cos(self.q[0])
		C2 = np.cos(self.q[1])
		S2 = np.sin(self.q[1])
		C12 = np.cos(self.q[0] + self.q[1])

		# generate entries of M(q) mass matrix
		M11 = (self.K1 + self.K2*C2)
		M12 = (self.K3 + self.K4*C2)
		M21 = M12
		M22 = self.K3

		# generate coriolis forces matrix V(q,dq)
		#print "dq: " + str(self.dq)
		V1 = -self.K4*S2*(2.0*self.dq[0]*self.dq[1] +self.dq[1]**2.0) 
		V2 = self.K4*S2*self.dq[0]**2.0

		# generate gravity forces matrix G(q)
		g = -9.81
		G1 = (self.m1 + self.m2)*g*self.l1*C1 + self.m2*g*self.l2*C12
		G2 = self.m2*g*self.l2*C12

		# computation without gravity terms
		#ddq1 = (M11*u[1] - M21*u[0] + M21*V1 - M11*V2) / (M11*M22 - M12**2.0)
		#ddq0 = (u[0] - M12*ddq1 - V1 ) / M11

		# computation WITH gravity terms
		#ddq1 = (M11*u[1] - M21*u[0] + M21*V1 + M21*G1 - M11*V2 - M11*G2) / (M11*M22 - M12**2.0)
		ddq1 = (M12*u[0] - M12*V1 - M12*G1 + M11*(V2 + G2 - u[1])) / (M12**2.0 - M22*M11)        
		ddq0 = (u[0] - M12*ddq1 - V1 - G1) / M11

		self.dq += np.array([ddq0, ddq1]) * dt
		self.q += self.dq * dt

		# transfer to next time step
		self.t += dt

	def joint_callback(self, msg):
		"""
		Save joint state message from ROS.
		"""
		self.joints_msg = msg
		print self.joints_msg.position[0], self.joints_msg.position[1]

	def to_JointState(self):
		"""
		Converts current configuration into a joint state message for ROS.
		"""
		msg = sensor_msgs.msg.JointState()
		msg.header.stamp = rospy.Time.now()
		msg.header.frame_id = ''
		msg.name = ['joint0', 'joint1']
		msg.position = [self.q[0], self.q[1]]
		msg.velocity = [self.dq[0],self.dq[1]]
		return msg

if __name__ == '__main__':
	# setup 2 link arm
	init_q = [2.04279802,  1.38453601]
	init_dq = [0.0, 0.0]
	l1 = 0.4
	l2 = 0.4
	dt = 1e-2
	arm = TwoLinkArm(init_q, init_dq, l1, l2, dt=dt)

	# create a ros node
	rospy.init_node('two_link_arm')
	rate = rospy.Rate(100)

	# set up subscribers to the joint state and publishers for markers
	#joint_state_sub = rospy.Subscriber('/joint_states', JointState, arm.joint_callback)
	joint_state_pub = rospy.Publisher('/joint_states', JointState, queue_size=10)
	marker_pub = rospy.Publisher('/visualization_marker', Marker, queue_size=10)

	print "\nEnter <x-dir force> <y-dir force> to apply a force at the end-effector:"
	while not rospy.is_shutdown():
		if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
			line = raw_input()
			splt = line.split()
			if len(splt) == 2:
				fin = map(float, splt)
				arm.fEE = np.array(fin)
				print "Force applied: ", arm.fEE
				f = arm.gen_Fq(forceEE=arm.fEE)
				q = arm.q
				normf = np.linalg.norm(f)
				normq = np.linalg.norm(q)
				direction  = np.dot(f, q)
				print "Direction: ", direction
				print "Config (before): ", arm.q

		        # if user entered a force at EE, apply it to torque
				arm.apply_torque(f, dt=arm.dt)
				print "Config (after): ", arm.q

		# update the visualization with the current joint state
		joint_state_pub.publish(arm.to_JointState())

		"""
		if arm.joints_msg is not None:
			q = [arm.joints_msg.position[0], arm.joints_msg.position[1]]
			xy = arm.fwd_kinematics(q)
			#print "xy: ", xy
			#print "ee xy:", (xy[0][2], xy[1][2]) 

			marker = Marker()
			marker.header.frame_id = 'root'
			marker.ns = 'end_effector_frame'
			marker.id = 0
			marker.type = Marker.SPHERE
			marker.pose.position.x = 0
			marker.pose.position.y = -xy[1][2] 
			marker.pose.position.z = xy[0][2]
			
			marker.pose.orientation.x = 0
			marker.pose.orientation.y = 0
			marker.pose.orientation.z = 0
			marker.pose.orientation.w = 1 
			marker.scale.x = 0.2
			marker.scale.y = 0.2
			marker.scale.z = 0.2
			marker.color.r = 0.4
			marker.color.g = 0.0
			marker.color.b = 0.8
			marker.color.a = 0.5
			marker_pub.publish(marker)
		"""

	rate.sleep()