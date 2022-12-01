import numpy as np
from numpy.core.fromnumeric import size
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt


class SR:
    '''
    Parameters
    ----------
    x : float
        Average distance to 'n' nearest neighbours for a given point.
    n_pts : int
        Number of datapoints recorded for the space.
    n_maxima : int
        Number of maxima under the threshold set as a 'valid' result.
    d : float, optional
        (default = 1.0) Maximum diameter of the maxima wells.
    '''
    
    '''def __init__(
        self, av_dist:float, n_pts:int, n_maxima:int, d_well:float=1.0
    ) -> None:
        
        self.av_dist = av_dist
        self.n_pts = n_pts
        self.n_maxima = n_maxima
        self.d = d_well
        
        self.solve()
    '''    
   
    def bind_to_df(x_iters, y_vals):
        df = pd.DataFrame(x_iters, columns=['x1','x2','x3'])
        df['y'] = y_vals
        return df


    def successes_dst(df, threshold:float, d_well):
        arr = None
        arr_x = None
        arr_y = None
        container = None
        n = 0   # used to overcome an issue
                # where an arr of length 1 reports length as 3 with len()
                # as it has a lower dimensionality than arrays with multiple sets of coords

        # unpack from df any datapoints below threshold
        for index, row in df.iterrows():
            if row['y'] <= threshold:
                entry = list([row[df.columns[0]],row[df.columns[1]],row[df.columns[2]]])
                n+=1

                # add to arr
                if np.all(arr == None):
                    arr = np.array(entry)
                else:
                    arr = np.vstack([arr,entry])
        
        # if no datapoints found below threshold, report nil values
        if np.all(arr == None):
                arr_new = np.array([])
                n_successes = 0.0

        # else
        else:
            # create first array for point-by-point comparison
            for _ in range(n):
                if np.all(arr_x == None):
                    arr_x = arr
                else:
                    arr_x = np.vstack([arr_x, arr])

            # create second array for point-by-point comparison
            for i in range(n):
                if arr.shape == (3,):
                    item = arr
                else:
                    item = arr[i]
                container = item

                for _ in range(n-1):
                    container = np.vstack([container,item])

                if np.all(arr_y == None):
                    arr_y = container
                else:
                    arr_y = np.vstack([arr_y, container])

            # reshape
            arr_x = arr_x.reshape(n,n,3)
            arr_y = arr_y.reshape(n,n,3)

            # take euclidean distance between all sets of points in arr
            arr_new = (arr_x-arr_y)**2

            # reshape and sum up all values for each datapoint below threshold
            arr_new = arr_new.reshape(n**2,3)
            arr_new = np.sum(arr_new, axis=1)
            arr_new = arr_new.reshape(n,n)
            
            # replace any values above distance, d_well with distance, d_well itself 
            arr_new = np.where(arr_new>d_well, d_well, arr_new)
            
            # calculate number of successes
            n_successes = len(arr_new)

        return arr_new, n_successes


    def average_distance(arr_row, nearest_neighbours:int=2):
        arr_row = sorted(arr_row)
        arr_sum = arr_row[1:nearest_neighbours+1]
        av_dist = np.sum(arr_sum)/nearest_neighbours
        return av_dist


    def success_rate(av_dist, n_successes, d_well, toggle_abs=True):
        '''
        Estimates the percentage and absolute number of successes in a search based on well size.

        Parameters
        ----------
        toggle_abs : bool, optional
            (default = True) Toggles whether absolute success rate is returned.
            Relative success rate is always returned.
        '''
        k = 4.61512051684
        rel_success_rate = (np.exp(k*av_dist/d_well)-1)

        if toggle_abs != True:
            abs_success_rate = n_successes * rel_success_rate/100
            return rel_success_rate, abs_success_rate
        else:
            return rel_success_rate


    def solve(x_iters, y_vals, threshold, nearest_neighbours, d_well):
        data = SR.bind_to_df(x_iters, y_vals)
        distances, n_successes = SR.successes_dst(data, threshold, d_well)
        
        averages=[]
        for i in range(len(distances)):
            pt_av_dist = SR.average_distance(distances[i], nearest_neighbours)
            averages.append(pt_av_dist)
        
        av_dist = np.sum(averages)/len(averages)

        relative, absolute = SR.success_rate(av_dist, n_successes, d_well, False)

        return relative, absolute


class matrix():

    def plot_6D_mechanics(
        ax, data, colors=['indigo','darkviolet','violet','blue','DarkCyan','g']
    ):
        cols = data.columns

        color1 = colors[0]       # cluster + cluster
        color2 = colors[1]       # cluster + template
        color3 = colors[2]       # cluster + solvent
        color4 = colors[3]       # template + template
        color5 = colors[4]       # solvent + template
        color6 = colors[5]       # solvent + solvent

        # Co3O(OH)
        rec1 = ax[5,0].add_patch(matplotlib.patches.Rectangle((-0.05,-0.05), 1.1, 1.1, color='lightgrey'))
        rec1.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[0], color=color1, s=2, ax=ax[5,0], zorder=2)
        data.plot.scatter(x=cols[1], y=cols[0], color=color1, s=2, ax=ax[5,1])
        data.plot.scatter(x=cols[2], y=cols[0], color=color1, s=2, ax=ax[5,2])
        data.plot.scatter(x=cols[4], y=cols[0], color=color2, s=2, ax=ax[5,3])
        data.plot.scatter(x=cols[6], y=cols[0], color=color2, s=2, ax=ax[5,4])
        data.plot.scatter(x=cols[7], y=cols[0], color=color3, s=2, ax=ax[5,5])

        # Co3O
        rec2 = ax[4,1].add_patch(matplotlib.patches.Rectangle((-0.05,-0.05), 1.1, 1.1, color='lightgrey'))
        rec2.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[1], color=color1, s=2, ax=ax[4,0])
        data.plot.scatter(x=cols[1], y=cols[1], color=color1, s=2, ax=ax[4,1], zorder=2)
        data.plot.scatter(x=cols[2], y=cols[1], color=color1, s=2, ax=ax[4,2])
        data.plot.scatter(x=cols[4], y=cols[1], color=color2, s=2, ax=ax[4,3])
        data.plot.scatter(x=cols[6], y=cols[1], color=color2, s=2, ax=ax[4,4])
        data.plot.scatter(x=cols[7], y=cols[1], color=color3, s=2, ax=ax[4,5])

        # Co4O4
        rec3 = ax[3,2].add_patch(matplotlib.patches.Rectangle((-0.05,-0.05), 1.1, 1.1, color='lightgrey'))
        rec3.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[2], color=color1, s=2, ax=ax[3,0])
        data.plot.scatter(x=cols[1], y=cols[2], color=color1, s=2, ax=ax[3,1])
        data.plot.scatter(x=cols[2], y=cols[2], color=color1, s=2, ax=ax[3,2], zorder=2)
        data.plot.scatter(x=cols[4], y=cols[2], color=color2, s=2, ax=ax[3,3])
        data.plot.scatter(x=cols[6], y=cols[2], color=color2, s=2, ax=ax[3,4])
        data.plot.scatter(x=cols[7], y=cols[2], color=color3, s=2, ax=ax[3,5])

        # rT1
        rec4 = ax[2,3].add_patch(matplotlib.patches.Rectangle((-0.05,-0.05), 0.6, 0.6, color='lightgrey'))
        rec4.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[4], color=color2, s=2, ax=ax[2,0])
        data.plot.scatter(x=cols[1], y=cols[4], color=color2, s=2, ax=ax[2,1])
        data.plot.scatter(x=cols[2], y=cols[4], color=color2, s=2, ax=ax[2,2])
        data.plot.scatter(x=cols[4], y=cols[4], color=color4, s=2, ax=ax[2,3], zorder=2)
        data.plot.scatter(x=cols[6], y=cols[4], color=color4, s=2, ax=ax[2,4])
        data.plot.scatter(x=cols[7], y=cols[4], color=color5, s=2, ax=ax[2,5])

        # rT2
        rec5 = ax[1,4].add_patch(matplotlib.patches.Rectangle((-0.05,-0.05), 0.6, 0.6, color='lightgrey'))
        rec5.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[6], color=color2, s=2, ax=ax[1,0])
        data.plot.scatter(x=cols[1], y=cols[6], color=color2, s=2, ax=ax[1,1])
        data.plot.scatter(x=cols[2], y=cols[6], color=color2, s=2, ax=ax[1,2])
        data.plot.scatter(x=cols[4], y=cols[6], color=color4, s=2, ax=ax[1,3])
        data.plot.scatter(x=cols[6], y=cols[6], color=color4, s=2, ax=ax[1,4], zorder=2)
        data.plot.scatter(x=cols[7], y=cols[6], color=color5, s=2, ax=ax[1,5])

        # V_TOT
        rec6 = ax[0,5].add_patch(matplotlib.patches.Rectangle((4.9,4.9), 5.2, 5.2, color='lightgrey'))
        rec6.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[7], color=color3, s=2, ax=ax[0,0])
        data.plot.scatter(x=cols[1], y=cols[7], color=color3, s=2, ax=ax[0,1])
        data.plot.scatter(x=cols[2], y=cols[7], color=color3, s=2, ax=ax[0,2])
        data.plot.scatter(x=cols[4], y=cols[7], color=color5, s=2, ax=ax[0,3])
        data.plot.scatter(x=cols[6], y=cols[7], color=color5, s=2, ax=ax[0,4])
        data.plot.scatter(x=cols[7], y=cols[7], color=color6, s=2, ax=ax[0,5], zorder=2)


    def plot_6D(data):
        fig, ax = plt.subplots(6,6, figsize=(11,11), sharex='col', sharey='row')
        
        matrix.plot_6D_mechanics(ax, data)


    def over_6D(data1, data2):
        fig, ax = plt.subplots(6,6, figsize=(11,11), sharex='col', sharey='row')
        
        matrix.plot_6D_mechanics(ax, data1, colors=['b','b','b','b','b','b'])
        matrix.plot_6D_mechanics(ax, data2, colors=['r','r','r','r','r','r'])


    def plot_3D_mechanics(
        ax, data, color='k'#, x_min=0, x_max=1
    ):
        cols = data.columns

        # x1
        #rec1 = ax[2,0].add_patch(matplotlib.patches.Rectangle((x_min-0.05,x_min-0.05), x_max+0.1, x_max+0.1, color='lightgrey'))
        #rec1.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[0], color=color, s=2, ax=ax[2,0], zorder=2)
        data.plot.scatter(x=cols[1], y=cols[0], color=color, s=2, ax=ax[2,1])
        data.plot.scatter(x=cols[2], y=cols[0], color=color, s=2, ax=ax[2,2])

        # x2
        #rec2 = ax[1,1].add_patch(matplotlib.patches.Rectangle((x_min-0.05,x_min-0.05), x_max+0.1, x_max+0.1, color='lightgrey'))
        #rec2.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[1], color=color, s=2, ax=ax[1,0])
        data.plot.scatter(x=cols[1], y=cols[1], color=color, s=2, ax=ax[1,1], zorder=2)
        data.plot.scatter(x=cols[2], y=cols[1], color=color, s=2, ax=ax[1,2])

        # x3
        #rec3 = ax[0,2].add_patch(matplotlib.patches.Rectangle((x_min-0.05,x_min-0.05), x_max+0.1, x_max+0.1, color='lightgrey'))
        #rec3.set_zorder=1
        data.plot.scatter(x=cols[0], y=cols[2], color=color, s=2, ax=ax[0,0])
        data.plot.scatter(x=cols[1], y=cols[2], color=color, s=2, ax=ax[0,1])
        data.plot.scatter(x=cols[2], y=cols[2], color=color, s=2, ax=ax[0,2], zorder=2)


    def over_3D(
        data1, data2#, x_min=0, x_max=1
    ):
        fig, ax = plt.subplots(3,3, figsize=(7,7), sharex='col', sharey='row')
        
        matrix.plot_3D_mechanics(ax, data1, color='b')#, x_min=x_min, x_max=x_max)
        matrix.plot_3D_mechanics(ax, data2, color='r')#, x_min=x_min, x_max=x_max)


class plot3D():
    def format_data(x_iters, y_vals, threshold):
        df = SR.bind_to_df(x_iters, y_vals)

        df_successes = pd.DataFrame(columns=df.columns)
        df_failures = pd.DataFrame(columns=df.columns)

        for index, row in df.iterrows():
            if row['y'] <= threshold:
                DICT = {df.columns[0]:row[df.columns[0]],
                        df.columns[1]:row[df.columns[1]],
                        df.columns[2]:row[df.columns[2]],
                        df.columns[3]:row[df.columns[3]]}
                
                df_successes = df_successes.append(DICT, True)
                
            else:
                DICT = {df.columns[0]:row[df.columns[0]],
                        df.columns[1]:row[df.columns[1]],
                        df.columns[2]:row[df.columns[2]],
                        df.columns[3]:row[df.columns[3]]}
                
                df_failures = df_failures.append(DICT, True)

        print(
            '# OF SUCCESSES:', len(df_successes), 
            '\n# OF FAILURES:', len(df_failures), 
            '\n TOTAL # POINTS:', len(df_successes + df_failures)
        )

        return df_successes, df_failures