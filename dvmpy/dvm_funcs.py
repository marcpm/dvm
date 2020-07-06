import numpy as np
import argparse



def arg_parser():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--naca", required=True, help="NACA 4digit")
    ap.add_argument("-a", "--alpha", required=True, type=int, help="Angle of attack (degrees)")
    ap.add_argument("-g", "--geometry", required=True, help="uniform or cosine")
    ap.add_argument("-p", "--panels", required=True, help="Number of Panels")
    ap.add_argument("-fp", "--flappos", required=False, type =int, default=False, help="flap position cord")
    ap.add_argument("-fa", "--flapalpha", required=False, type =int, default=False, help="flap deflection angle")
    args = vars(ap.parse_args())
    return args


def airfoil_discretization(args):
    # Takes CLI parsed args or dictionary as input
    
    naca = str(args["naca"]) # convert to string to separate digits easier
    alpha = args["alpha"]
    geo = args["geometry"]
    npanels = args["panels"]
    flap_pos = args["flappos"]
    flap_alpha = args["flapalpha"]
    
    max_camber = int(naca[0])/100 
    max_camber_pos = int(naca[1])/10
    thick = int(naca[2:])
    npoints = npanels + 1
    alpha = np.deg2rad(alpha)
        
    if geo == "uniform":
        X = np.linspace(0,1, num=npoints, endpoint=True)
    elif geo == "fullcosine":
        _a = np.arange(1,npoints+1) # arange func stops at n-1
        X = 0.5 * (1 - np.cos(((_a-1) / npanels) * np.pi))
        
    if flap_pos != False:
        flap_alpha = np.deg2rad(flap_alpha)
            
    Z = np.zeros(npoints)    
    for i in range(npoints):
        if X[i] <= max_camber_pos:
            Z[i] = (max_camber/(max_camber_pos**2)) * (2 * max_camber_pos*(X[i])- (X[i]**2))
            
        elif X[i] > max_camber_pos:
             Z[i] = (max_camber/(1-max_camber_pos)**2)*(1-2*max_camber_pos+2*max_camber_pos*(X[i])-(X[i]**2)) 
                
        # Flap
        if X[i] >= flap_pos:
            Z[i] = Z[i] - np.tan(flap_alpha) * (X[i] - flap_pos)
    
    XZ_space = np.vstack([X, Z])    
    return XZ_space  


def dvm_solver(args):
    alpha = args["alpha"]
    npanels = args["panels"]
    npoints = npanels + 1
    alpha = np.deg2rad(alpha)
    
    XZ_space = airfoil_discretization(args)
        
    vortex_points = np.zeros((2,npanels))
    ctrl_points = np.zeros((2,npanels))  
    for i in range(npanels): 
        # panel length
        cordlen = np.sqrt( (XZ_space[0][i+1] - XZ_space[0][i])**2 +
                           (XZ_space[1][i+1] - XZ_space[1][i])**2 )
        
        angle = np.arctan( (XZ_space[1][i+1] - XZ_space[1][i]) / 
                            (XZ_space[0][i+1] - XZ_space[0][i])   )
        
        vortex_points[0][i] = XZ_space[0][i] + 1/4 * cordlen * np.cos(angle)
        vortex_points[1][i] = XZ_space[1][i] + 1/4 * cordlen * np.sin(angle)
        
        ctrl_points[0][i] = XZ_space[0][i] + 3/4 * cordlen * np.cos(angle)
        ctrl_points[1][i] = XZ_space[1][i] + 3/4 * cordlen * np.sin(angle)
        
        
    inf_coeff = np.zeros((npanels,npanels))    
    RHS = np.zeros((npanels,1))
    for i in range(npanels): #add npanels -1 ??? last i+1 gonna crash
        XZ_deltas = np.zeros((2,1))
        #XZ_deltas = XZ_space[:,i+1] - XZ_space[:,i]
        XZ_deltas[0] = XZ_space[0][i+1] - XZ_space[0][i]
        XZ_deltas[1] = XZ_space[1][i+1] - XZ_space[1][i]
        
        tan_vector = np.zeros((2,1))
        #tan_vector = XZ_deltas / np.sqrt(XZ_deltas[]**2)
        tan_vector[0] = XZ_deltas[0] / np.sqrt( XZ_deltas[0]**2 + XZ_deltas[1]**2 )
        tan_vector[1] = XZ_deltas[1] / np.sqrt( XZ_deltas[0]**2 + XZ_deltas[1]**2 )
        
        norm_vector = np.zeros((2,1))
        norm_vector[0] = - tan_vector[1]
        norm_vector[1] = tan_vector[0]
        
        for j in range(npanels):
            XZ_deltas[0] = ctrl_points[0][i] - vortex_points[0][j]
            XZ_deltas[1] = ctrl_points[1][i] -vortex_points[1][j]
            
            r = np.sqrt(XZ_deltas[0]**2 + XZ_deltas[1]**2 )
            
            induced_vel = np.zeros((2,1))
            induced_vel[0] = 1/(2*np.pi)* XZ_deltas[1]/(r)**2
            induced_vel[1] = -1/(2*np.pi)* XZ_deltas[0]/(r)**2
            
            inf_coeff[i][j] = induced_vel[0] * norm_vector[0] + induced_vel[1]*norm_vector[1]
            
        RHS[i] = -( np.cos(alpha)*norm_vector[0] + np.sin(alpha)* norm_vector[1] )
        
    gamma = np.linalg.solve( inf_coeff, RHS)
    # flat plate test
    #gamba = gamma/np.sin(alpha)*5/np.pi

    cl = 2*np.sum(gamma)
    
    _sm = np.zeros((npanels, 1))
    for i in range(npanels):
        _sm[i] = gamma[i]*vortex_points[0][i]*np.cos(alpha)
    
    
    cmle = -2*np.sum(_sm)
    
    return np.array([[cl,cmle]])