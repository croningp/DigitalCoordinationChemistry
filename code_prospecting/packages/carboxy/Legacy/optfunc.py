###########################################################################################
# 1| NECESSARY LIBRARIES                                                                  #
###########################################################################################

# maths and data handling
import numpy as np
import pandas as pd

# visualisation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import seaborn as sns
import PyQt5

# random number generation
import random

# efficiency and timing
import itertools
import timeit



###########################################################################################
# 2| GENERAL FUNCTIONS                                                                    #
###########################################################################################

def grid_1d (llim:float=-5, ulim:float=5, num:int=200):
    '''
    Returns an incremented variable list. Equivalent to numpy.linspace()
    '''
    x_array = np.linspace(llim, ulim, num)  
    return x_array

"""
def grid_2d_sq (llim=-5, ulim=5, num=50):
    # returns two square matrices (of dimension num x num)
    x1 = np.linspace(llim, ulim, num)
    x1_elongate = np.array([])
    for i in range(num):
        x1_elongate = np.concatenate((x1_elongate, x1))
    x1_array = x1_elongate.reshape((num,num))
    
    x2 = np.linspace(llim, ulim, num)
    x2_elongate = np.array([])
    for i in range(num):
        for j in range(num):
            x2_elongate = np.append(x2_elongate, x2[i])
    x2_array = x2_elongate.reshape((num,num))
    
    return (x1_array, x2_array)


def grid_2d_irr (llim1=-5, ulim1=5, num1=50, llim2=-5, ulim2=5, num2=50):
    # returns two rectangular matrices (of dimension num1 x num2)
    x1 = np.linspace(llim1, ulim1, num1)
    x1_elongate = np.array([])
    for i in range(num2):
        x1_elongate = np.concatenate((x1_elongate, x1))
    x1_array = x1_elongate.reshape((num1,num2))
    
    x2 = np.linspace(llim2, ulim2, num2)
    x2_elongate = np.array([])
    for i in range(num2):
        for j in range(num1):
            x2_elongate = np.append(x2_elongate, x2[i])
    x2_array = x2_elongate.reshape((num1,num2))
    
    return (x1_array, x2_array)
"""

def generate_grid(variables:list, file_name:str="variable"):
    '''
    Transforms lists of incremented variables to individual arrays corresponding to the values
    of each axis in a grid composed of points at every combination of the variable increments.
    Taken together these arrays define each point in the 'n'-dimensional space.
    Generates as many arrays as there are variables provided.
    Each array is returned as an individual python variable.
    '''
    num_list = []
    variable_list = []
    n=0
    
    for i in variables:
        num = int(len(i))
        num_list.append(len(i))
    
    for i in range(len(variables)):
        array = np.array([])
        for t in itertools.product(*variables):
            array = np.append(array, t[-i])   # needs to be -ve to return the grids in the right
                                              # order for the display2d class to interpret them
        new_array = array.reshape(num_list)
        variable_list.append(new_array)

        n += 1
        print(n, " array(s) of ", len(num_list), " completed.")

        # save each "new_array" variable generated as a .txt file
        with open((file_name+"_x"+str(n)+".txt"), 'w') as outfile:
            outfile.write('# Array shape: {0}\n'.format(new_array.shape))

            for data_slice in new_array:
                np.savetxt(outfile, data_slice, fmt='%-7.2f')
                outfile.write('# New slice\n')

    return variable_list


def find_nearest(value:float, array):
    '''
    Returns nearest value present in an array. 
    Requires array to be values from a plotted mathematical function as opposed to the product 
    of the 'generate_grid' function.
    '''
    # returns nearest value present in an array to that entered as an argument 
    array = np.asarray(array)
    idx = (abs(array-value)).argmin()
    return array.flat[idx]


def func_min(function):
    '''
    Undefined.
    '''
    True

def normalise_func(function):
    '''
    Undefined.
    '''
    True



###########################################################################################
# 3| DISPLAY OPTIONS                                                                      #
###########################################################################################

class displaynd:

    def matrix(z, xi:list):
        True

class display2d:
    '''
    Provides options to plot arrays of points in 3-dimensions, where 2 given dimensions are 
    input variables, and the third an output.
    '''
    def scatter(z, x1, x2):
        '''
        Produces a scatter graph to demonstrate the exact position of each sampled point for 
        a function. Does not accurately represent relative scales of x1 and x2.
        '''
        fig = plt.figure(figsize=(9,6))
        ax = fig.gca(projection='3d')
        ax.view_init(elev=40, azim=-60)
        ax.set_zlabel('f(variable1, variable2)')
        ax.set_xlabel('variable1')
        ax.set_ylabel('variable2')
        ax.scatter(x1, x2, z, s=1, c='k', marker='o', alpha=1)
    
    def surface(z, x1, x2, colormap="jet"):
        '''
        Produces a surface built from the sampling of a function at each point in the provided
        grid. Where too few points are given, surface returned may not be smooth. Does not 
        accurately represent relative scales of x1 and x2.
        '''
        fig = plt.figure(figsize=(9,6))
        ax = fig.gca(projection='3d')
        ax.view_init(elev=40, azim=-60)
        ax.set_zlabel('f(variable1, variable2)')
        ax.set_xlabel('variable1')
        ax.set_ylabel('variable2')
        ax.plot_trisurf(x1.flatten(), x2.flatten(), z.flatten(), cmap=colormap)
    
    def contour(z, x1, x2, colormap="jet", labels=True, scale=True, levels:int=9):
        '''
        Produces a contour plot built from the sampling of a function at each point in the provided
        grid. Accurately represents the relative scales of x1 and x2. May not always show the
        true minima - depending on the density of points given as input.
        '''
        # create table containing function values
        df = pd.DataFrame(z, index=x1[:,0], columns=x2[0])
        
        x1_range = max(x1[:,0]) - min(x1[:,0])
        x2_range = max(x2[0]) - min(x2[0])
        factor = max([x1_range, x2_range])/6
        
        if scale==True:
            fig, ax = plt.subplots(figsize=(x1_range/factor, x2_range/factor))
        elif scale==False:
            fig, ax = plt.subplots(figsize=(6,6))
        CS = ax.contour(x1, x2, df, levels=levels, cmap=colormap)
        ax.set_xlabel('variable1')
        ax.set_ylabel('variable2')
        if labels == True:
            ax.clabel(CS, inline=True, fontsize=10, fmt='%1.1f')

    def heatmap(z, x1, x2, colormap="jet"):
        '''
        Produces a heatmap of provided samples of a function. Does not accurately represent 
        relative scales of x1 and x2.
        '''
        df = pd.DataFrame(z, index=np.around(x1[0],2), columns=np.around(x2[:,0],2))
        fig = plt.figure()
        sns.heatmap(df)
    
class display1d:
    '''
    Provides options to plot arrays of points in 2-dimensions, where 1 given dimension is an 
    input variable, and the other an output.
    '''    
    def line(y, x):
        '''
        Plots a line graph from the sampling of a function. Where there are too few points 
        present, may be some warping.
        '''
        fig = plt.figure(figsize=(9,6))
        #ax = add_axes(ax)
        #ax.set_ylabel('f(x)')
        #ax.set_xlabel('x')
        #ax.
        plt.plot(x.flatten(), y.flatten(), c='k')

class explore:
    '''
    To be added.
    '''
    def map_contour(z, x1, x2, y1, y2, colormap="jet", labels=False, scale=True):
        '''
        To be added.
        '''
        display2d.contour(z, x1, x2, colormap=colormap, labels=labels, scale=scale)
        annotation = range(len(y1))
        plt.plot(y1, y2, 'ko', ms=12)
        for i, txt in enumerate(annotation):
            plt.annotate(txt, (y1[i], y2[i]), xytext=(-3,-4), textcoords='offset pixels', color='white')



###########################################################################################
# 4| DISPLAY OPTIONS                                                                      #
###########################################################################################

def ackley(xi:list, a:float=20, b:float=0.2, c:float=2*np.pi):
    n = len(xi)

    f2_inputs = 0
    f5_inputs = 0

    for x in xi:
        f1 = (x**2)
        f2_inputs += f1
        f3 = (np.cos(c * x))
        f5_inputs += f3

    f2 = (1/n) * f2_inputs
    f4 = -b * np.sqrt(f2)
    f5 = (1/n) * f5_inputs
    f6 = -a * np.exp(f4) - np.exp(f5) + a + np.exp(1)
    return f6


def alpine1(xi:list):
    n = len(xi)
    f2 = 0

    for x in xi:
        f1 = abs(x * np.sin(x) + (0.1 * x))
        f2 += f1
    return f2


def alpine2_r(xi:list):
    n = len(xi)
    f2 = 1
    
    for x in xi:
        f1 = np.sqrt(abs(x)) * np.sin(x)
        f2 *= f1
    return -f2


def brown(xi:list):
    n = len(xi)
    f6 = 0
    
    for x in xi[0:(n-1)]:
        for xj in xi[1:n]:
            f1 = x**2
            f2 = xj**2
            f3_inputs = f2+1
            f3 = f1**(f3_inputs)
            f4_inputs = f1+1
            f4 = f2**(f4_inputs)
            f5 = f3+f4
            f6 += f5
    return f6
    

def exponential(xi:list):
    n = len(xi)
    f2_inputs = 0
    
    for x in xi:
        f1 = x**2
        f2_inputs += f1

    f2 = -0.5 * f2_inputs
    f3 = - np.exp(f2)
    return f3


def griewank(xi:list, a=1, b=4000):
    n = len(xi)
    i = 1
    f3a_inputs = 0
    f3b_inputs = 1

    for x in xi:
        f1 = (x**2)/b
        f3a_inputs += f1
        f2 = np.cos(x/(np.sqrt(i)))
        f3b_inputs *= f2
        i += 1

    f3 = f3a_inputs - f3b_inputs + a
    return f3


def happycat(xi:list, alpha:float=0.125):
    n = len(xi)
    i = 1
    sum_x = 0
    sum_x2 = 0
    
    for x in xi:
        sum_x += x
        sum_x2 += x**2
    
    norm_x = np.sqrt(sum_x2)
    f1 = ((norm_x**2) - n)**2
    f2 = f1**alpha
    f3 = (0.5*(norm_x**2) + sum_x)/n 
    f4 = f2 + f3 + 0.5
    return f4   


def periodic(xi:list, a=1, b=0.1):
    n = len(xi)
    f3_inputs = 0
    f4_inputs = 0
    
    for x in xi:
        f1 = 0.5 *(1-(np.cos(2*x)))   # should be equivalent to sin^2(x)
        f2 = x**2
        f4_inputs += f1
        f3_inputs += f2
    
    f3 = np.exp(f3_inputs)
    f4 = a + f4_inputs - (0.1 * f3)
    return f4


def powellsum(xi:list):
    n = len(xi)
    i = 1
    f2 = 0
    
    for x in xi:
        f1 = abs(x)**(i+1)
        f2 += f1
        i += 1
    return f2


def qing(xi:list):
    n = len(xi)
    i = 1
    f3 = 0
    
    for x in xi:
        f1 = (x**2) - i
        f2 = f1**2
        f3 += f2
        i += 1
    return f3


def quartic(xi:list):
    n = len(xi)
    i = 1
    f3 = 0
    
    for x in xi:
        rand = random.uniform(0,1)
        f1 = x**4
        f2 = (i * f1) + rand
        i += 1
        f3 += f2
    return f3


def rastrigin(xi:list, a:float=10, b:float=2):
    n = len(xi)
    sum_3 = 0

    for item in xi:
        sum_1 = a * np.cos(b * np.pi * item)
        sum_2 = item**2 - sum_1
        sum_3 += sum_2

    f_1 = 10 * n + sum_3
    return f_1


def salomon(xi:list, a:float=1, b:float=0.1, c:float=2*np.pi):
    n = len(xi)
    f1 = 0
    
    for x in xi:
        f1 += x**2
    
    f2 = np.sqrt(f1)
    f3 = a - np.cos(c * f2) + (0.1 * f2)
    return f3


def schwefel(xi:list, a:float=418.9829):
    n = len(xi)
    sum_2 = 0
    for item in xi:
        sum_1 = item * np.sin(np.sqrt(abs(item)))
        sum_2 += sum_1
    f_1 = a * n - sum_2
    return f_1
    
    
def schwefel220(xi:list):
    n = len(xi)
    f = 0
    
    for x in xi:
        f += abs(x)
    return f


def schwefel222(xi:list):
    n = len(xi)
    f1a_inputs = 0
    f1b_inputs = 1
    
    for x in xi:
        f1a_inputs += abs(x)
        f1b_inputs *= abs(x)
        
    f1 = f1a_inputs + f1b_inputs
    return f1


def schwefel223(xi:list):
    n = len(xi)
    f = 0    
    
    for x in xi:
        f += x**10
    return f


def shubert(xi:list):
    n = len(xi)
    f2 = 0
    f3 = 1
    
    for x in xi:
        for j in range(1,6):
            f1 = (j + 1) * x + j
            f2 += np.cos(f1)
        f3 *= f2
        f2 = 0
    return f3


def shubert3(xi:list):
    n = len(xi)
    f3_inputs = 0
    f3 = 0
    
    for x in xi:
        for j in range(1,6):
            f1 = (j+1)*x
            f2 = j * np.sin(f1+j)
            f3_inputs += f2
        f3 += f3_inputs
    return f3


def shubert4(xi:list):
    n = len(xi)
    f3_inputs = 0
    f3 = 0
    
    for x in xi:
        for j in range(1,6):
            f1 = (j+1)*x
            f2 = j * np.cos(f1+j)
            f3_inputs += f2
        f3 += f3_inputs
    return f3


def sphere(xi:list):
    n = len(xi)
    f = 0
    
    for x in xi:
        f += x**2
    return f


def styblinskitank(xi:list):
    n = len(xi)
    f4_inputs = 0
    
    for x in xi:
        f1 = x**4
        f2 = 16 * (x**2)
        f3 = 5 * x
        f4_inputs += f1 - f2 + f3
    
    f4 = 0.5 * f4_inputs
    return f4


def sumsquares(xi:list):
    n = len(xi)
    f = 0
    i = 1
    
    for x in xi:
        f += i * (x**2)
        i += 1
    return f


def xinsheyang(xi:list):
    n = len(xi)
    i = 1
    f2 = 0

    for x in xi:
        epsilon = random.uniform(0,1)
        f1 = epsilon * (abs(x)**i)
        f2 += f1
        i += 1
    return f2


def xinsheyang2(xi:list):
    n = len(xi)
    f3_inputs = 0
    f4_inputs = 0
    
    for x in xi:
        f1 = abs(x)
        f2 = np.sin(x**2)
        f4_inputs += f1 
        f3_inputs += f2
      
    f3 = np.exp(-f3_inputs)
    f4 = f4_inputs * f3
    return f4


def xinsheyang3(xi:list, m:float=5, beta:float=15):
    n = len(xi)
    f4_inputs = 0
    f5_inputs = 0
    f6_inputs = 1    
    
    for x in xi:
        f1 = (x / beta)**(2*m)
        f2 = x**2
        f3 = 0.5 *(1+(np.cos(2*x)))   # should be equivalent to cos^2(x)
        f4_inputs += f1
        f5_inputs += f2
        f6_inputs *= f3
        
    f4 = np.exp(-f4_inputs)
    f5 = 2 * np.exp(-f5_inputs)
    f6 = f4 - (f5 * f6_inputs)
    return f6


def xinsheyang4(xi:list):
    n = len(xi)
    f4a_inputs = 0
    f4b_inputs = 0
    f5_inputs = 0
    
    for x in xi:
        f1 = 0.5 *(1-(np.cos(2*x)))   # should be equivalent to sin^2(x) 
        f2 = x**2
        f4a_inputs += f1
        f4b_inputs += f2
        f3 =  0.5 *(1-(np.cos(2*(np.sqrt(abs(x))))))   # should be equivalent to sin^2(sqrt(x))
        f5_inputs += f3
                
    f4 = f4a_inputs - np.exp(-f4b_inputs)
    f5 = f4 * np.exp(-f5_inputs)
    return f5

    
def zakharov(xi:list):
    n = len(xi)
    i = 1
    f1 = 0
    f34_inputs = 0
    
    for x in xi:
        f1 += x**2
        f2 = 0.5 * i * x
        f34_inputs += f2
        
        i += 1
        
    f3 = f34_inputs**2
    f4 = f34_inputs**4
    f5 = f1 + f3 + f4
    return f5
  


###########################################################################################
# 5| SEARCH AND HELP FUNCTIONALITY                                                        #
###########################################################################################

functions = {
    "ackley" : {
        "name": "Ackley Function",
        "class": "localminima",
        "global minima": [0,0],
        "input domain": [-10,10],
        "continuous": True,
        "convex": False,
        "differentiable": True,
        "separable": False,
        "scalable": True,
        "modality": "multimodal",
        "dimensionality": "n"
    },
    "alpine1" : {
        "name" : "Alpine Function N.1",
        "class" : "localminima",
        "global minima" : [0],
        "input domain" : [0,10],
        "continuous" : np.nan,
        "convex" : False,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : np.nan,
        "dimensionality" : "n"
    },
    "alpine2_r" : {
        "name" : "Alpine Function N.2",
        "class" : "localminima",
        "global minima" : [7.917],
        "input domain" : [0,10],
        "continuous" : np.nan,
        "convex" : False,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : np.nan,
        "dimensionality" : "n"
    },
    "exponential" : {
        "name" : "Exponential Function",
        "class" : "highD",
        "global minima" : [0],
        "input domain" : [-1,1],
        "continuous" : True,
        "convex" : True,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "griewank" : {
        "name" : "Griewank Function",
        "class" : "localminima",
        "global minima" : [0],
        "input domain" : [-600,600],
        "continuous" : True,
        "convex" : False,
        "differentiable" : np.nan,
        "separable" : np.nan,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "periodic" : {
        "name" : "Periodic Function",
        "class" : "localminima",
        "global minima" : [0],
        "input domain" : [-10,10],
        "continuous" : True,
        "convex" : False,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "powellsum" : {
        "name" : "Powell Sum Function",
        "class" : "highD",
        "global minima" : [0],
        "input domain" : [-1,1],
        "continuous" : True,
        "convex" : True,
        "differentiable" : False,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "qing" : {
        "name" : "Qing Function",
        "class" : "highD",
        "global minima" : "+/- sqrt(i)",
        "input domain" : [-500,500],
        "continuous" : True,
        "convex" : False,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "quartic" : {
        "name" : "Quartic Function",
        "class" : "highD",
        "global minima" : [0],
        "input domain" : [-1.28,1.28],
        "continuous" : True,
        "convex" : False,
        "differentiable" : True,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "rastrigin" : {
        "name" : "Rastrigin Function",
        "class" : "localminima",
        "global minima" : "f=0 at (0,0)",
        "input domain" : "x[i] = [-5.12,5.12]",
        "continuous" : True,
        "convex" : True,
        "differentiable" : True,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "salomon" : {
        "name" : "Salomon Function",
        "class" : "highD",
        "global minima" : [0],
        "input domain" : [-100,100],
        "continuous" : True,
        "convex" : False,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "schwefel" : {
        "name" : "Schwefel Function",
        "class" : "localminima",
        "global minima" : [420.9687],
        "input domain" : [-500,500],
        "continuous" : True,
        "convex" : False,
        "differentiable" : np.nan,
        "separable" : np.nan,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "schwefel220" : {
        "name" : "Schwefel Function 2.20",
        "class" : "highD",
        "global minima" : [0],
        "input domain" : [-100,100],
        "continuous" : True,
        "convex" : True,
        "differentiable" : False,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "schwefel222" : {
        "name" : "Schwefel Function 2.22",
        "class" : "highD",
        "global minima" : [0],
        "input domain" : [-100,100],
        "continuous" : True,
        "convex" : True,
        "differentiable" : False,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "schwefel223" : {
        "name" : "Schwefel Function 2.23",
        "class" : "bowl-shaped",
        "global minima" : [0],
        "input domain" : [-10,10],
        "continuous" : True,
        "convex" : True,
        "differentiable" : False,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "shubert" : {
        "name" : "Shubert Function",
        "class" : "bowl-shaped",
        "global minima" : [0],
        "input domain" : [-5.12,5.12],
        "continuous" : True,
        "convex" : True,
        "differentiable" : True,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "sphere" : {
        "name" : "Sphere Function",
        "class" : "localminima",
        "global minima" : "18 global minima of approx. f=-186.7309",
        "input domain" : [-10,10],
        "continuous" : True,
        "convex" : False,
        "differentiable" : True,
        "separable" : False,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "styblinskitank" : {
        "name" : "Styblinski-Tank Function",
        "class" : "highD",
        "global minima" : [-2.903534],
        "input domain" : [-5,5],
        "continuous" : True,
        "convex" : False,
        "differentiable" : np.nan,
        "separable" : np.nan,
        "scalable" : np.nan,
        "modality" : "multimodal",
        "dimensionality" : "n"
    },
    "sumsquares" : {
        "name" : "Sum Squares Function",
        "class" : "bowl-shaped",
        "global minima" : [0],
        "input domain" : [-10,10],
        "continuous" : True,
        "convex" : True,
        "differentiable" : True,
        "separable" : True,
        "scalable" : np.nan,
        "modality" : "unimodal",
        "dimensionality" : "n"
    },
    "xinsheyang" : {
        "name": "Xin-She Yang Function",
        "class": "highD",
        "global minima": [0],
        "input domain": [-5,5],
        "continuous": np.nan,
        "convex": False,
        "differentiable": False,
        "separable": True,
        "scalable": np.nan,
        "modality": np.nan,
        "dimensionality": "n"    
    },
    "xinsheyang2" : {
        "name": "Xin-She Yang Function N.2",
        "class": "highD",
        "global minima": [0],
        "input domain": [-2*np.pi,2*np.pi],
        "continuous": np.nan,
        "convex": False,
        "differentiable": False,
        "separable": False,
        "scalable": np.nan,
        "modality": np.nan,
        "dimensionality": "n" 
    },
    "xinsheyang3" : {
        "name": "Xin-She Yang Function N.3",
        "class": "highD",
        "global minima": [0],
        "input domain": [-2*np.pi,2*np.pi],
        "continuous": np.nan,
        "convex": True,
        "differentiable": True,
        "separable": False,
        "scalable": np.nan,
        "modality": np.nan,
        "dimensionality": "n" 
    },
    "xinsheyang4" : {
        "name": "Xin-She Yang Function N.4",
        "class": "highD",
        "global minima": [0],
        "input domain": [-10,10],
        "continuous": np.nan,
        "convex": False,
        "differentiable": False,
        "separable": False,
        "scalable": np.nan,
        "modality": np.nan,
        "dimensionality": "n" 
    },
    "zakharov" : {
        "name": "Zakharov Function",
        "class": "plate-shaped",
        "global minima": [0],
        "input domain": [-5,10],
        "continuous": True,
        "convex": True,
        "differentiable": np.nan,
        "separable": np.nan,
        "scalable": np.nan,
        "modality": "unimodal",
        "dimensionality": "n" 
    },
}

def list_func():
    print("AVAILABLE FUNCTIONS: ", [*functions])
    print("\nNumber of Available Functions:", len([*functions]), "\n")
    
    """
    list_of_functions = []
    for key in [*functions]:
        list_of_functions.append(key)
    print(list_of_functions)
    """

def lookup_func(search_term):
    found = 0
    for key in [*functions]:
        if key == search_term:
            print("LOOKUP_FUNCTION:", search_term, functions[key])
            found += 1
    if found < 1:
        print("Error: Search term not found.")
    elif found > 1:
        print("Error: Multiple results found.")