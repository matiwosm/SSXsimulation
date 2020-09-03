import numpy as np
import matplotlib.pyplot as plt

def makeTrajectoryPuncture(trajectoryArray, dimension, value, precision, min1,
max1, min2, max2):
    '''
    Takes in an array and creates a trajectory puncture plot
    Arguments:
        trajectoryArray: (3xN) numpy array where N is len of trajectory.
            dimensions are ordered x,y,z
        dimension: int (0:x, 1:y, 2:z) along which to slice to make plot
        value: number at which to slice for plot
        precision: number indicating how far from value to find points
        min1, max1, min2, max2: number; extents of puncture plot on the two
        non-"dimension" axes
    Return: void (makes plot)
    '''
    dict = {0:'x', 1:'y', 2:'z'}

    dim1Num = (dimension+1)%3
    dim2Num = (dimension+2)%3

    if(trajectoryArray.shape[0]==3):
        dim = trajectoryArray[dimension]
        dim1 = trajectoryArray[dim1Num]
        dim2 = trajectoryArray[dim2Num]
    else:
        st = str(trajectoryArray.shape)
        raise ValueError('array must have shape 3xN. array has shape %s' % st)

    dim1Array = np.array([])
    dim2Array = np.array([])
    for i in range(len(dim)):
        if np.abs(dim[i]-value) < precision:
            new = np.array([dim1[i]])
            dim1Array = np.append(dim1Array, new)
            new = np.array([dim2[i]])
            dim2Array = np.append(dim2Array, new)

    plt.plot(dim1Array, dim2Array, 'bo')
    plt.ylabel(dict[dim2Num])
    plt.xlabel(dict[dim1Num])
    plt.ylim((min2, max2))
    plt.xlim((min1, max1))
    plt.title("Dipole Trajectory Puncture Plot: %s vs. %s" % (dict[dim1Num], dict[dim2Num]))
    plt.axis('equal')
    plt.show()

data = np.load('path_4d.npy')
array= np.swapaxes(data, 1, 0)
array = array[:3]

makeTrajectoryPuncture(array, 2, 0, .001, -2, 2, -2, 2)
