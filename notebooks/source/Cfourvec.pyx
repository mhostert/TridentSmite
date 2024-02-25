#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#cython: boundscheck=False
#cython: language_level=3
#cython: wraparound=False
#cython: nonecheck=False
#cython: define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#cython: cdivision=True

import numpy as np
cimport numpy as np
from numpy cimport ndarray

#######################################
# C functions to be used
from libc.math cimport sqrt, abs, log, cos, sin, acos

#######################################
# C implementation of RANDOM
from libc.stdlib cimport rand, RAND_MAX

cdef double UniformRand():
	return ( rand() + 1. )/(RAND_MAX + 1. )



#######################################
# FOURVECTOR FUNCTIONS
#######################################

#******************************
def build_fourvec(ndarray[double, ndim=1] E, ndarray[double, ndim=1] p, ndarray[double, ndim=1] cost, ndarray[double, ndim=1] phi):

	cdef int i,m
	m = phi.shape[0]
	cdef ndarray[double,ndim=2] s = np.empty((m,4))

	with nogil:
		for i in range(m):
			s[i,0] = E[i]
			s[i,1] = p[i]*cos(phi[i])*sqrt(1.0-cost[i]*cost[i])
			s[i,2] = p[i]*sin(phi[i])*sqrt(1.0-cost[i]*cost[i])
			s[i,3] = p[i]*cost[i]

	return s

#******************************
def momentum_scalar(ndarray[double] E, double mass):

	cdef int i,m
	m = E.shape[0]
	cdef ndarray[double,ndim=1] s = np.empty((m))

	with nogil:
		for i in range(m):
			s[i] = sqrt(E[i]*E[i] - mass*mass)
	return s

#******************************
def get_theta_3vec(ndarray[double, ndim=2] r):

	cdef int i,m
	m = r.shape[0]
	cdef ndarray[double, ndim=1] s = np.empty((m))

	with nogil:
		for i in range(m):
			s[i] = acos(r[i,3]/sqrt(r[i,1]*r[i,1]+r[i,2]*r[i,2]+r[i,3]*r[i,3]))
	return s

#******************************
def dot4(ndarray[double, ndim=2] x, ndarray[double, ndim=2] y):
	cdef int i,m
	m= x.shape[0]
	cdef ndarray[double,ndim=1] s = np.empty((m))
	with nogil:
		for i in range(m):
			s[i] = x[i,0]*y[i,0] - x[i,1]*y[i,1] - x[i,2]*y[i,2] - x[i,3]*y[i,3]
	return  s

#******************************
def dot3(ndarray[double, ndim=2] x, ndarray[double, ndim=2] y):
	cdef int i,m
	m= x.shape[0]
	cdef ndarray[double,ndim=1] s = np.empty((m))
	for i in range(m):
		s[i] = x[i,1]*y[i,1] + x[i,2]*y[i,2] + x[i,3]*y[i,3]
	return  s

#******************************
def get_cosTheta(ndarray[double, ndim=2] x):
	cdef int i,m
	m= x.shape[0]
	cdef ndarray[double, ndim=1] s = np.empty((m))
	with nogil:
		for i in range(m):
			s[i] = x[i,3]/sqrt(x[i,1]*x[i,1] + x[i,2]*x[i,2] + x[i,3]*x[i,3])
	return s

#******************************
def rotationx(ndarray[double, ndim=2] v4, ndarray[double] theta):

	cdef int i, m;

	m = v4.shape[0]

	cdef ndarray[double, ndim=2] res = np.empty((m,4))
	cdef ndarray[double, ndim=3] R = np.empty((m,4,4))

	with nogil:
		for i in range(m):
			R[i,0,0] = 1.0
			R[i,0,1] = 0.0
			R[i,0,2] = 0.0
			R[i,0,3] = 0.0

			R[i,1,0] = 0.0
			R[i,1,1] = 1.0
			R[i,1,2] = 0.0
			R[i,1,3] = 0.0

			R[i,2,0] = 0.0
			R[i,2,1] = 0.0
			R[i,2,2] = cos(theta[i])
			R[i,2,3] = -sin(theta[i])

			R[i,3,0] = 0.0
			R[i,3,1] = 0.0
			R[i,3,2] = sin(theta[i])
			R[i,3,3] = cos(theta[i])

			res[i,0] = R[i,0,0]*v4[i,0] + R[i,0,1]*v4[i,1] + R[i,0,2]*v4[i,2] + R[i,0,3]*v4[i,3]
			res[i,1] = R[i,1,0]*v4[i,0] + R[i,1,1]*v4[i,1] + R[i,1,2]*v4[i,2] + R[i,1,3]*v4[i,3]
			res[i,2] = R[i,2,0]*v4[i,0] + R[i,2,1]*v4[i,1] + R[i,2,2]*v4[i,2] + R[i,2,3]*v4[i,3]
			res[i,3] = R[i,3,0]*v4[i,0] + R[i,3,1]*v4[i,1] + R[i,3,2]*v4[i,2] + R[i,3,3]*v4[i,3]
		    
	return res

# #******************************
# def rotationy(ndarray[double, ndim=2] v4, ndarray[double] theta):

# 	cdef int i, m;

# 	m = v4.shape[0]

# 	cdef ndarray[double, ndim=2] res = np.empty((m,4))
# 	cdef ndarray[double, ndim=3] R = np.empty((m,4,4))

# 	with nogil:
# 		for i in range(m):
# 			R[i,0,0] = 1.0
# 			R[i,0,1] = 0.0
# 			R[i,0,2] = 0.0
# 			R[i,0,3] = 0.0

# 			R[i,1,0] = 0.0
# 			R[i,1,1] = cos(theta[i])
# 			R[i,1,2] = 0.0
# 			R[i,1,3] = -sin(theta[i])

# 			R[i,2,0] = 0.0
# 			R[i,2,1] = 0.0
# 			R[i,2,2] = 1.0
# 			R[i,2,3] = 0.0

# 			R[i,3,0] = 0.0
# 			R[i,3,1] = sin(theta[i])
# 			R[i,3,2] = 0.0
# 			R[i,3,3] = cos(theta[i])

# 			res[i,0] = R[i,0,0]*v4[i,0] + R[i,0,1]*v4[i,1] + R[i,0,2]*v4[i,2] + R[i,0,3]*v4[i,3]
# 			res[i,1] = R[i,1,0]*v4[i,0] + R[i,1,1]*v4[i,1] + R[i,1,2]*v4[i,2] + R[i,1,3]*v4[i,3]
# 			res[i,2] = R[i,2,0]*v4[i,0] + R[i,2,1]*v4[i,1] + R[i,2,2]*v4[i,2] + R[i,2,3]*v4[i,3]
# 			res[i,3] = R[i,3,0]*v4[i,0] + R[i,3,1]*v4[i,1] + R[i,3,2]*v4[i,2] + R[i,3,3]*v4[i,3]
		    
# 	return res
#******************************
def rotationy(ndarray[double, ndim=2] v4, ndarray[double] theta):

	cdef int i, m;
	m = v4.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))

	with nogil:
		for i in range(m):
			res[i,0] = v4[i,0]
			res[i,1] = cos(theta[i])*v4[i,1] - sin(theta[i])*v4[i,3]
			res[i,2] = v4[i,2]
			res[i,3] = sin(theta[i])*v4[i,1] + cos(theta[i])*v4[i,3]
		    
	return res
#******************************
def rotationy_sin(ndarray[double, ndim=2] v4, ndarray[double] stheta):

	cdef int i, m;
	m = v4.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))

	with nogil:
		for i in range(m):
			res[i,0] = v4[i,0]
			res[i,1] = sqrt(1.0-stheta[i]*stheta[i])*v4[i,1] - stheta[i]*v4[i,3]
			res[i,2] = v4[i,2]
			res[i,3] = stheta[i]*v4[i,1] + sqrt(1.0-stheta[i]*stheta[i])*v4[i,3]
		    
	return res
#******************************
def rotationy_cos(ndarray[double, ndim=2] v4, ndarray[double] ctheta, int sign=1):

	cdef int i, m;
	m = v4.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))

	with nogil:
		for i in range(m):
			res[i,0] = v4[i,0]
			res[i,1] = ctheta[i]*v4[i,1] - sign*sqrt(1.0-ctheta[i]*ctheta[i])*v4[i,3]
			res[i,2] = v4[i,2]
			res[i,3] = sign*sqrt(1.0-ctheta[i]*ctheta[i])*v4[i,1] + ctheta[i]*v4[i,3]
		    
	return res

#******************************
def rotationz(ndarray[double, ndim=2] v4, ndarray[double] theta):

	cdef int i, m;
	m = v4.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))
	with nogil:
		for i in range(m):
			res[i,0] = v4[i,0]
			res[i,1] = cos(theta[i])*v4[i,1] - sin(theta[i])*v4[i,2] 
			res[i,2] = sin(theta[i])*v4[i,1] + cos(theta[i])*v4[i,2] 
			res[i,3] = v4[i,3]
		    
	return res
#******************************
def rotationz_sin(ndarray[double, ndim=2] v4, ndarray[double] stheta):

	cdef int i, m;
	m = v4.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))
	with nogil:
		for i in range(m):
			res[i,0] = v4[i,0]
			res[i,1] = sqrt(1.0-stheta[i]*stheta[i])*v4[i,1] - stheta[i]*v4[i,2] 
			res[i,2] = stheta[i]*v4[i,1] + sqrt(1.0-stheta[i]*stheta[i])*v4[i,2] 
			res[i,3] = v4[i,3]
		    
	return res
#******************************
def rotationz_cos(ndarray[double, ndim=2] v4, ndarray[double] ctheta, int sign = 1):

	cdef int i, m;
	m = v4.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))
	with nogil:
		for i in range(m):
			res[i,0] = v4[i,0]
			res[i,1] = ctheta[i]*v4[i,1] - sign*sqrt(1.0-ctheta[i]*ctheta[i])*v4[i,2] 
			res[i,2] = sign*sqrt(1.0-ctheta[i]*ctheta[i])*v4[i,1] + ctheta[i]*v4[i,2] 
			res[i,3] = v4[i,3]
		    
	return res
#******************************
# def rotationz(ndarray[double, ndim=2] v4, ndarray[double] theta):

# 	cdef int i, m;

# 	m = v4.shape[0]

# 	cdef ndarray[double, ndim=2] res = np.empty((m,4))
# 	cdef ndarray[double, ndim=3] R = np.empty((m,4,4))

# 	with nogil:
# 		for i in range(m):
# 			R[i,0,0] = 1.0
# 			R[i,0,1] = 0.0
# 			R[i,0,2] = 0.0
# 			R[i,0,3] = 0.0

# 			R[i,1,0] = 0.0
# 			R[i,1,1] = cos(theta[i])
# 			R[i,1,2] = -sin(theta[i])
# 			R[i,1,3] = 0.0

# 			R[i,2,0] = 0.0
# 			R[i,2,1] = sin(theta[i])
# 			R[i,2,2] = cos(theta[i])
# 			R[i,2,3] = 0.0

# 			R[i,3,0] = 0.0
# 			R[i,3,1] = 0.0
# 			R[i,3,2] = 0.0
# 			R[i,3,3] = 1.0

# 			res[i,0] = R[i,0,0]*v4[i,0] + R[i,0,1]*v4[i,1] + R[i,0,2]*v4[i,2] + R[i,0,3]*v4[i,3]
# 			res[i,1] = R[i,1,0]*v4[i,0] + R[i,1,1]*v4[i,1] + R[i,1,2]*v4[i,2] + R[i,1,3]*v4[i,3]
# 			res[i,2] = R[i,2,0]*v4[i,0] + R[i,2,1]*v4[i,1] + R[i,2,2]*v4[i,2] + R[i,2,3]*v4[i,3]
# 			res[i,3] = R[i,3,0]*v4[i,0] + R[i,3,1]*v4[i,1] + R[i,3,2]*v4[i,2] + R[i,3,3]*v4[i,3]
		    
# 	return res

#******************************
def L(ndarray[double, ndim=2] v4, ndarray[double] beta):
	cdef int i, m;
	m = beta.shape[0]
	cdef ndarray[double, ndim=2] res = np.empty((m,4))
	with nogil:
		for i in range(m):
			res[i,0] = 1.0/sqrt(1.0 - beta[i]*beta[i])*v4[i,0] - beta[i]/sqrt(1.0 - beta[i]*beta[i])*v4[i,3]
			res[i,1] = v4[i,1] 
			res[i,2] = v4[i,2] 
			res[i,3] = -beta[i]/sqrt(1.0 - beta[i]*beta[i])*v4[i,0] + 1.0/sqrt(1.0 - beta[i]*beta[i])*v4[i,3]
		    
	return res

#******************************
def T(ndarray[double, ndim=2] v4, ndarray[double] beta, ndarray[double] theta, ndarray[double] phi):
	return L( rotationy( rotationz(v4,-phi), theta), -beta)

#******************************
def Tinv(ndarray[double, ndim=2] v4, ndarray[double] beta, ndarray[double] ctheta, ndarray[double] phi):
	return rotationz( rotationy_cos( L(v4, beta), ctheta, sign=-1), phi)
