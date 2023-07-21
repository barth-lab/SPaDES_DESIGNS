#! /usr/bin/env python
import numpy as np
import math 

#''' transform matrix format ''
# [1,0,0]
# [0,1,0]
# [0,0,1] ---- rotation matrix
# [v,p,n] ---- translation vector

def is_near_zero(a):
    return a < 1E-6

def Heron_Formulate(a, b, c):
    a, b, c = map(float, [a,b,c])
    CS = (a**2 + b**2 - c**2)/(2*a*b)
    C = math.acos( CS )
    return C

def Coordinate_Center(Coordinate_set):
    center_coordinate = sum( Coordinate_set )/len(Coordinate_set)
    return center_coordinate

def Vector_generate(Coord1,Coord2):
    v = np.array(Coord2) -  np.array(Coord1)
    return v

def DistanceMeasure(TupleA,TupleB):
    d = math.sqrt((TupleA[0]-TupleB[0])**2+(TupleA[1]-TupleB[1])**2+(TupleA[2]-TupleB[2])**2)
    return(d)

def Distance_Search( TupleA, Matrix_select, treshold ):
    I,J = Matrix_select.shape
    d_list = []
    for i in range(I) :
        d = DistanceMeasure( Matrix_select[i], TupleA )
        d_list.append( d )
    d_list = np.array(d_list)
    J = d_list <= treshold
    return J

def Matrix_Select( TupleA, MatrixB, treshold):
    '''find all the row in MatrixB closer to TupleA than the treshold'''
    #Poplace = TupleA + np.array( [treshold, treshold, treshold] )
    #Neplace = TupleA - np.array( [treshold, treshold, treshold] )
    J_list  = []
    for i in [0,1,2]:
        J_list.append( np.absolute( MatrixB[:,i] - TupleA[i] ) <= treshold  )
        #Criteria = [ Poplace[0,i], Neplace[0:i] ].sort()
        #J = dstack( (MatrixB[i] > Criteria[0], MaitrxB[i] < Criteria[1]) ) 
        #J_list.append( J.np.array( [any(j) for j in D[0] ] ) )
    J = np.dstack( tuple( J_list) )
    J = np.array( [any(j) for j in J[0] ] )
    return MatrixB[J], J

def Squared_DistanceMeasure( TupleA, TupleB):
    d = (TupleA[0]-TupleB[0])**2+(TupleA[1]-TupleB[1])**2+(TupleA[2]-TupleB[2])**2
    return d

def Angle_Between_Vector(u,v,lu,lv):
    # '''  angle is 0->pi
    # u and v are vectors, lu and lv are their length respectively'''
    a1 = math.fsum(list(u*v))
    a2 = lu*lv
    a = math.acos( a1/a2 )
    return a

def Scale_Vector(u):
    '''find the unit vector to the same direction '''
    scale = math.sqrt( u[0]**2+ u[1]**2 + u[2]**2 )
    return u/scale

def Transform_Vector(m, v):
    '''m is a matrix that could transform vector v'''
    x = sum( m[:3,0]*v ) + m[3,0]
    y = sum( m[:3,1]*v ) + m[3,1]
    z = sum( m[:3,2]*v ) + m[3,2]
    d = np.array([x,y,z])
    return d

def RotationAtOrigin(axis, angle):
    '''it generate a matrix to rotate a coordinate with angle and center at the axis  '''
    v = Scale_Vector( axis )

    c = math.cos(float(angle))
    s = math.sin(float(angle))
    t = 1.0 - c
    
    m = np.zeros(12)
    m.shape = 4,3
    m[0,0] = t*v[0]*v[0] + c
    m[0,1] = t*v[0]*v[1] + v[2]*s
    m[0,2] = t*v[0]*v[2] - v[1]*s
  
    m[1,0] = t*v[1]*v[0] - v[2]*s
    m[1,1] = t*v[1]*v[1] + c
    m[1,2] = t*v[1]*v[2] + v[0]*s

    m[2,0] = t*v[2]*v[0] + v[1]*s
    m[2,1] = t*v[2]*v[1] - v[0]*s
    m[2,2] = t*v[2]*v[2] + c
 
    return m

def PlaneFitTLSTSQ( coordinate ):
    ''' using total least-squares (TLSTSQ) solution to give a plane non-curve to fit a set of points n >= 3
        At least three points and not in a line
        assume: ax + by + cz +d = 0 
        we need S is minimal 
        method: SVD singular value decomposition
        coordinate is an array of atom coordinates
    '''
    A, center = Set_Center_Zero( coordinate )
    u,s,vh = np.linalg.svd(A)
    v = vh.conj().transpose()
    return v[:,-1], center

def PointPlaneDis( n, c, v ):
    '''  a(x-c1) + b(y-c2) + c(z-c3) = 0 
         v = np.array([a,b,c])
         '''
    dv = n - c
    D = math.fabs( v[0]*dv[0] + v[1]*dv[1] + v[2]*dv[2] )/np.sqrt( (v[0]**2) + (v[1]**2) +(v[2]**2) )
    #D =( v[0]*n[0] + v[1]*n[1] - n[2] + v[2])/np.sqrt( ((v[0]**2) + (v[1]**2) +1))
    return D

def Set_Center_Zero(v):
    center = np.array( [Coordinate_Center(v[:,0]), Coordinate_Center(v[:,1]), Coordinate_Center(v[:,2])] )
    v = v-center
    return v, center

