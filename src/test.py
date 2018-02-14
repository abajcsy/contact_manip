#! /usr/bin/env python
import numpy as np
import rospy
import geometry_msgs.msg
import std_msgs.msg
import sensor_msgs.msg

if __name__ == '__main__':
	rospy.init_node("manipulation")

	# create pose publisher
	pose_pub = rospy.Publisher('/joint_states', sensor_msgs.msg.JointState, queue_size=1)

	msg = sensor_msgs.msg.JointState()
	msg.header.stamp = rospy.Time.now()
	msg.header.seq = 0
	msg.header.frame_id = ''
	msg.name = ['base_joint', 'joint1', 'joint2']
	msg.position = [0,0.5,0]
	msg.velocity = [0,0,0]
	msg.effort = [0,0,0]

	# publish to ROS at 100hz
	r = rospy.Rate(100) 
	while not rospy.is_shutdown():

		pose_pub.publish(msg)
		r.sleep()
