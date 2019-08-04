'''
This code contains the following functions:
    theta: convert chromosome accumulative coords to theta
    complex: convert chromosome coords to complex array list
    arc: return a list of complex array which can be used to draw an arc based on coordinate intervals
    heatmap_val: convert input values to heatmap values, used for plotting heatmap
'''

import numpy as np
import colorlover as cl

def to_theta(accum_coord, SUM, degreerange=[0,360]):
    # SUM=sum(ideogram_chr_config['chr_size'])+sum(ideogram_chr_config['chr_spacing'])
    # accum_coord could be np.ndarray, float or a list
    # returns theta value (same data type as accum_coord)
    PI = np.pi
    a = PI*(degreerange[1]-degreerange[0])/(SUM*180)
    
    if not isinstance(accum_coord, list):
        # the constant to convert accumulated coordinates to theta
        theta = accum_coord*a
        theta += degreerange[0]*(PI/180)
    else: 
        theta = [*map(lambda x: x*a + degreerange[0]*(PI/180), accum_coord)]
    return theta


def to_complex(theta, radius):
    '''this function considers all possible data type of theta and radius and returns the complex needed for path and trace'''
    # always output a flatten ndarray or a list of flatten ndarray
    # verticalize the theta and radius when horizontal expansion is needed, but flatten at the end! In the case of bezier curves, we might need to reshape the array

    if isinstance(radius, (float, int)):
        if isinstance(theta, np.ndarray):
            
            Complex = np.zeros(theta.shape, dtype='complex')
            Complex.real = np.sin(theta.astype('float'))*radius
            Complex.imag = np.cos(theta.astype('float'))*radius
            return Complex

        elif isinstance(theta, list):
            Complex = [*map(lambda x: to_complex(x, radius), theta)]
            return Complex

    elif isinstance(radius, np.ndarray):
        # for ideogram, heatmap, cytoband etc where radius.shape[1] > theta.shape[1]
        
        if isinstance(theta, np.ndarray):
            
            if theta.ndim == 1:
                theta = theta.reshape((len(theta), 1))
            if radius.ndim == 1:
                radius = radius.reshape((len(radius), 1))

            if theta.shape[1] == radius.shape[1]:
                Complex = np.zeros(theta.shape, dtype='complex')
                Complex.real = np.sin(theta.astype('float'))*radius
                Complex.imag = np.cos(theta.astype('float'))*radius

                return Complex.flatten()

            elif theta.shape[1] == 1 and theta.shape[1] < radius.shape[1]:
                theta = np.repeat(theta, radius.shape[1]).reshape((len(theta),radius.shape[1]))
                Complex = to_complex(theta, radius)
                return Complex
            elif radius.shape[1] == 1 and theta.shape[1] > radius.shape[1]:
                radius = np.repeat(radius, theta.shape[1]).reshape((len(radius),theta.shape[1]))
                Complex = to_complex(theta, radius)
                return Complex
        elif isinstance(theta, list):
            Complex = [*map(lambda x: to_complex(x, radius), theta)]
            return Complex

    elif isinstance(radius, list):
        if isinstance(theta, list):
            Complex = [*map(lambda x, y: to_complex(x, y), theta, radius)]
            return Complex

        elif isinstance(theta, np.ndarray):
            # this should not happen, but I'll convert radius list into radius ndarray
            radius = np.array(radius).reshape((1,len(radius)))
            Complex = to_complex(theta, radius)
            return Complex


def to_arc(interval_theta_array, arc_radius, ideogram_theta_list):
    # make sure the two ndarray are of the same length
    # the interval is inclusive on both ends
    # returns a list of arc complex array(each element of the list being a np.array of complex coordinates within the interval), 1X coordinates, if used for drawing rings, you need to repeat it
    # the arc is inclusive on both ends (start and end coordinates are always included, plus whatever you have in the ideogram)
    # ideogram_theta_list is a list of theta values (2X), they need to be split in half first

    try:
        assert np.all(np.logical_and(interval_theta_array>=0, interval_theta_array<=2*np.pi))
        assert np.all(interval_theta_array[:,0]-interval_theta_array[:,1]<=0)       # allow if only one data point on a given chr
    except AssertionError:
        print ('interval_theta_array_list: \t')
        print (interval_theta_array.tolist())
        print ('\t')

        print ('assertionError, printing interval_theta_array:\t')
        #print (np.where(interval_theta_array[:,0]-interval_theta_array[:,1]>0 ))
        print (interval_theta_array[np.where(interval_theta_array >= 2*np.pi)])
        print ('the problem index starts from:')
        print (np.where(interval_theta_array >= 2*np.pi))
        print ('interval_theta_array shape:')
        print (interval_theta_array.shape)

        #print (interval_theta_array)
        # DEBUG
        print ('arc_radius is:')
        print (arc_radius)

    if not isinstance(arc_radius, (float, int)):
        assert len(interval_theta_array) == len(arc_radius)
    
    ideogram_theta_list = [*map(lambda x: np.split(x, 2)[0], ideogram_theta_list)]          
    # now ideogram_theta_list is 1X
    ideogram_theta_array = np.concatenate(ideogram_theta_list)


    arc_theta_list, arc_complex_list = [], []
    for i in range(len(interval_theta_array)):
        arc_theta_list.append(ideogram_theta_array[np.where(np.logical_and(ideogram_theta_array>interval_theta_array[i,0], ideogram_theta_array<interval_theta_array[i,1]))])
        arc_complex_list.append(np.zeros(len(arc_theta_list[i])+2, dtype='complex'))
        try:
            arc_complex_list[i].real = np.sin(np.insert(interval_theta_array[i,:], 1, arc_theta_list[i]).astype(float))*arc_radius[i]
            arc_complex_list[i].imag = np.cos(np.insert(interval_theta_array[i,:], 1, arc_theta_list[i]).astype(float))*arc_radius[i]
        except Exception:
            arc_complex_list[i].real = np.sin(np.insert(interval_theta_array[i,:], 1, arc_theta_list[i]).astype(float))*arc_radius
            arc_complex_list[i].imag = np.cos(np.insert(interval_theta_array[i,:], 1, arc_theta_list[i]).astype(float))*arc_radius

    return arc_complex_list


def val2heatmap(input_val_array, palatte_dict={'palatte': 'RdBu',
                                               'scale': 'div',
                                               'reverse': True,
                                               'ncolor': 11}):


    # palatte, heatmap_scale='div', heatmap_ncol=11):
    # convert input values to heatmap values
    # heatmap values must be integer from range(9) for sequential scale or range(11) for divergent scale
    

    if not isinstance(input_val_array, np.ndarray):
        input_val_array = np.array(input_val_array, dtype='float')

    heatmap_val = np.zeros(input_val_array.shape, dtype='float')
    maximum = max(input_val_array)
    minimum = min(input_val_array)

    if palatte_dict['scale'] == 'div':
        heatmap_val[np.where(input_val_array[:]<0)] = (-5)*input_val_array[np.where(input_val_array[:]<0)]/minimum+5
        heatmap_val[np.where(input_val_array[:]>=0)] = (5)*input_val_array[np.where(input_val_array[:]>=0)]/maximum+5 
    elif palatte_dict['scale'] == 'seq':
        if maximum > 0 and minimum >= 0:
            heatmap_val = (9)*(input_val_array/maximum)
        elif maximum <= 0 and minimum < 0:
            heatmap_val = (-9)*(input_val_array/minimum)
        else:
            raise ValueError('Please make sure all data values either positive or negative if you choose sequential palatte')
    
    try:
        #print('debugging...')
        #print(palatte_dict)
        palatte_cl = cl.scales[str(palatte_dict['ncolor'])][palatte_dict['scale']][palatte_dict['palatte']]
    
    except Exception:
        try:
            assert palatte_dict['ncolor'] in range(4,12)
            if palatte_dict['ncolor'] > 9 and palatte_dict['scale'] == 'seq':
                print('Sequential palatte {} does not have more than 9 colors, setting ncolor to 9'.format(palatte_dict['palatte']))
                palatte_cl = cl.scales['9'][palatte_dict['scale']][palatte_dict['palatte']]
            else:
                if palatte_dict['scale'] == 'seq':
                    try:
                        palatte_cl = cl.scales[str(palatte_dict['ncolor'])]['div'][palatte_dict['palatte']]
                    except KeyError:
                        raise KeyError('palatte {} cannot be found in colorlover package!'.format(palatte_dict['palatte']))
                elif palatte_dict['scale'] == 'div':
                    try:
                        palatte_cl = cl.scales[str(palatte_dict['ncolor'])]['seq'][palatte_dict['palatte']]
                    except KeyError:
                        raise KeyError('palatte {} cannot be found in colorlover package!'.format(palatte_dict['palatte']))        

        except Exception:
            palatte_cl = palatte_dict['palatte']
 

    int_heatmap_val = [*map(lambda x: int(round(x)), heatmap_val)]
    if palatte_dict['reverse'] == False:
        int_heatmap_val = [*map(lambda x: palatte_dict['ncolor']-x-1, int_heatmap_val)]
    
    rgb_heatmap_val = [*map(lambda x: palatte_cl[x], int_heatmap_val)]

    return np.array(rgb_heatmap_val)

def val2radius(data_val, R0, R1):
    radius = np.zeros(data_val.shape, dtype='float')
    if type(R0) is not list:
        d = R1-R0
    else:
        d = [*map(lambda x, y: x-y, R1, R0)]
    maximum = max(data_val)

    ## ONGOING, NEED TO ADD possibility for negative values?

    try: 
        minimum = min(data_val)
        radius[:][np.where(data_val[:]>=0)] = data_val[np.where(data_val[:]>0)]*(d/maximum)+R0
        #radius[:][np.where(data_val[:]==0)] = R0
        radius[:][np.where(data_val[:]<0)] = data_val[np.where(data_val[:]<0)]*(d/minimum)*(-1)+R0
    except ZeroDivisionError:
        print (data_val)
        print (R0)
        print (R1)

    return radius

def bezier_complex(data_accum_coord, ends_radius, bezier_radius, SUM, type='link', degreerange=[0,360]):
    # input accum_accord data and output complex array!
    # need to calculate
    assert isinstance(data_accum_coord, np.ndarray)

    ends_theta = to_theta(data_accum_coord, SUM, degreerange)
    ctrlpoint_theta = np.zeros((len(ends_theta),2), dtype='float')

    if type == 'link' or type == 'twistedribbon':
        ctrlpoint_theta[:,0] = (ends_theta[:,0]+ends_theta[:,2])/2.0
        ctrlpoint_theta[:,1] = (ends_theta[:,1]+ends_theta[:,3])/2.0
    elif type == 'ribbon':
        ctrlpoint_theta[:,0] = (ends_theta[:,0]+ends_theta[:,3])/2.0
        ctrlpoint_theta[:,1] = (ends_theta[:,1]+ends_theta[:,2])/2.0
    else:
        raise ValueError('Please select a valid bezier curve type (link,ribbon,twistedribbon)')
    theta_array = np.column_stack((ends_theta[:,:2], ctrlpoint_theta, ends_theta[:,2:]))
    radius_array = np.ones(theta_array.shape)*ends_radius
    radius_array[:,2:4] = bezier_radius
    return to_complex(theta_array, radius_array)
