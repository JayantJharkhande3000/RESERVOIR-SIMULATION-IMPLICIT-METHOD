# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:40:39 2020

@author: Dell
"""


import numpy as np
import matplotlib.pyplot as plt 

print(   "\t ########################################################################### \t"  )
print(   "\t ############# Reservoir Simulation implicit solution to 1D  ############### \t"  )
print(   "\t ######## For dirichlet(constant pressure) boundary consition   ############ \t"  )
print(   "\t ############ left side constant pressure boundary consition ############### \t"  )
print(   "\t ################## Right side no flow  boundary consition   ############### \t"  )
print(   "\t ########################################################################### \t"  )
    



L = 10000
print("\n\tThe lenght of the reservoir is ", str(L))
n = 4 
print("\n\tBlock node in the reservoir is ", str(n))
P0 = 1000
print("\n\tThe intial pressure of the reservoir is ", str(P0))
P_left = 2000
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_left))
dx = L/n
print("\n\tThe lenght of Blocks in the reservoir is ", str(dx))

porosity = 0.2
print("\n\tthe porosity value of the reservoir is ", str(porosity))

permeability  = 50
print("\n\tthe permeability(md) value of the reservoir is ", str(permeability))

viscosity = 1
print("\n\tthe viscosity value is ", str(viscosity))

area = 200000
print("\n\tCross sectional area of the reservoir ", str(area))

compressibility = 1*10**(-6)
print("\n\tcompressibility of the reservoir is ", str(compressibility))

Bw = 1
print("\n\t water formation volume factor is " +str(Bw)+ " rb/stb" )


#################### final time for simulation is  ##############################
t_final = 3
print("\n\t the reservoir simulation time in days is  ", str(t_final))

#################### time step  ###############################################
dt = 1
print("\n\t the reservoir simulation incremental time step in days is "+ str(dt)+ "day")

########  x represent the lenght distribution array in the L direction with n nodes ########
x = np.linspace(dx/2, L-dx/2, n)


########  P represent the pressure distribution array in the reservoir with n nodes ########
P = np.ones((n,n), dtype=float)

###### dPdt represent the pressure_gradient update distribution array in the reservoir with n nodes ########
dPdt = np.empty(n)

############  T represent the time distribution array in with the dt increment #############
t =  np.arange(0 , t_final, dt)



alpha = float(permeability /(porosity*viscosity*compressibility))

print("\n\t the neta constant value is  ", str(alpha))

neta = float((alpha*dt)/(dx**2))
print("\n\t the neta constant value is  ", str(neta))



neta = 0.2532
print("\n\t the neta constant value is  ", str(neta))


print("\n############## pressure distribution ################\n")
print("pressure distribution at day 0 is\n", str(P))



"""
here M is the time step for large value of M the solution will converges

"""
M = 8
N = 4
pressure_previous = np.ones([N,1])*P0

boundary_condition_array =  np.zeros([N,1])




for k in range(0, M, 1):  # k only reachs M - 1, coz need to stop at t = T which is at index M
    # initilise the trigonal matrix*
    mat_dig = np.zeros([N,N])
    print("\n time step value",k+1)

    for i in range(1,N,1):
        #print("i value" , i)
        mat_dig[i][i] = 1 + 2 * neta
        mat_dig[i][i-1] = - neta
        mat_dig[i-1][i] = - neta
    mat_dig[0][0]= 1 + 3 * neta
    mat_dig[N-1][N-1]= 1 + neta
    boundary_condition_array[0][0] = 2*neta*P_left
    pressure_previous = np.dot(np.linalg.inv(mat_dig), (pressure_previous + boundary_condition_array))
    print("\n",pressure_previous)

    
    
print(boundary_condition_array.shape) 
   
print(boundary_condition_array)    
print(mat_dig)

print(pressure_previous.shape)

print("\n",pressure_previous)















