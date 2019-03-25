from bash_command import bash_command as bc
import rebound
import numpy as np

import end_aleksey_config_b


def test_rotate_1():
	ang=np.pi/2.0
	vec=np.array([1,0,0])
	axis=np.array([0,0,1])
	vec=np.array(end_aleksey_config_b.rotate_vec(ang,axis,vec))
	print vec, np.array([np.cos(ang),np.sin(ang),0])
	assert np.allclose(vec, np.array([np.cos(ang),np.sin(ang),0]))

def test_rotate_2():
	ang=np.pi/3.0
	vec=np.array([1,0,0])
	axis=np.array([0,0,1])
	vec=np.array(end_aleksey_config_b.rotate_vec(ang,axis,vec))
	print vec, np.array([np.cos(ang),-np.sin(ang),0])
	assert np.allclose(vec, np.array([np.cos(ang),np.sin(ang),0]))
