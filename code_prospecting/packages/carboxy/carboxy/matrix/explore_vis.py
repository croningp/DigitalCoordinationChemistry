### CODE: DANIEL KOWALSKI

#from re import I
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from sklearn import metrics
import spectra_visualiser as vis


def CCE_matrix(
    df_volumes, prev_res=None, 
    range_cluster:float=4.7, range_template:float=2.6, range_solv:float=6.0
):
    '''
    '''
    data = df_volumes
    cols = data.columns

    fig, ax = plt.subplots(10,10, figsize=(11,11), sharex='col', sharey='row')

    color1 = 'indigo'       # cluster + cluster
    color2 = 'darkviolet'   # cluster + template
    color3 = 'violet'       # cluster + solvent
    color4 = 'blue'         # template + template
    color5= 'DarkCyan'      # solvent + template
    color6 = 'g'            # solvent + solvent

    ####################################################################################
    
    if prev_res != None:
        previous_results_legible = pd.read_csv(prev_res, index_col=0)
        data_prev = previous_results_legible
        data_prev = data_prev.drop(["prime1","prime2"],axis=0)

        data_prev.replace(0.00, np.nan, inplace=True)
        cols_prev = data_prev.columns

        color_prev = 'grey'

        # Co3O(OH)
        rec1 = ax[9,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec1.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,0], zorder=2)
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[0], color=color_prev, s=1, ax=ax[9,9])

        # Co3O
        rec2 = ax[8,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec2.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,1], zorder=2)
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[1], color=color_prev, s=1, ax=ax[8,9])

        # Co4O4
        rec3 = ax[7,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec3.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,2], zorder=2)
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[2], color=color_prev, s=1, ax=ax[7,9])


        # Ce
        rec4 = ax[6,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec4.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,3], zorder=2)
        ax[6,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[6,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[3], color=color_prev, s=1, ax=ax[6,9])

        # Dy
        rec5 = ax[5,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec5.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,2])
        ax[5,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,4], zorder=2)
        ax[5,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[4], color=color_prev, s=2, ax=ax[5,9])

        # Yb
        rec6 = ax[4,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec6.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,2])
        ax[4,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[4,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,5], zorder=2)
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[5], color=color_prev, s=1, ax=ax[4,9])

        # OA
        rec7 = ax[3,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec7.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,6], zorder=2)
        ax[3,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[3,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[6], color=color_prev, s=1, ax=ax[3,9])

        # SA
        rec8 = ax[2,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec8.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,5])
        ax[2,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,7], zorder=2)
        ax[2,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[7], color=color_prev, s=1, ax=ax[2,9])

        # TMTACN
        rec9 = ax[1,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec9.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,5])
        ax[1,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[1,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,8], zorder=2)
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[8], color=color_prev, s=1, ax=ax[1,9])

        # MeOH
        rec10 = ax[0,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_solv, range_solv, color='lightgrey'))
        rec10.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[9], color=color_prev, s=1, ax=ax[0,9], zorder=2)

    ####################################################################################

    # Co3O(OH)
    rec1 = ax[9,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec1.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[0], color=color1, s=2, ax=ax[9,0], zorder=2)
    data.plot.scatter(x=cols[1], y=cols[0], color=color1, s=2, ax=ax[9,1])
    data.plot.scatter(x=cols[2], y=cols[0], color=color1, s=2, ax=ax[9,2])
    data.plot.scatter(x=cols[3], y=cols[0], color=color2, s=2, ax=ax[9,3])
    data.plot.scatter(x=cols[4], y=cols[0], color=color2, s=2, ax=ax[9,4])
    data.plot.scatter(x=cols[5], y=cols[0], color=color2, s=2, ax=ax[9,5])
    data.plot.scatter(x=cols[6], y=cols[0], color=color2, s=2, ax=ax[9,6])
    data.plot.scatter(x=cols[7], y=cols[0], color=color2, s=2, ax=ax[9,7])
    data.plot.scatter(x=cols[8], y=cols[0], color=color2, s=2, ax=ax[9,8])
    data.plot.scatter(x=cols[9], y=cols[0], color=color3, s=2, ax=ax[9,9])

    # Co3O
    rec2 = ax[8,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec2.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[1], color=color1, s=2, ax=ax[8,0])
    data.plot.scatter(x=cols[1], y=cols[1], color=color1, s=2, ax=ax[8,1], zorder=2)
    data.plot.scatter(x=cols[2], y=cols[1], color=color1, s=2, ax=ax[8,2])
    data.plot.scatter(x=cols[3], y=cols[1], color=color2, s=2, ax=ax[8,3])
    data.plot.scatter(x=cols[4], y=cols[1], color=color2, s=2, ax=ax[8,4])
    data.plot.scatter(x=cols[5], y=cols[1], color=color2, s=2, ax=ax[8,5])
    data.plot.scatter(x=cols[6], y=cols[1], color=color2, s=2, ax=ax[8,6])
    data.plot.scatter(x=cols[7], y=cols[1], color=color2, s=2, ax=ax[8,7])
    data.plot.scatter(x=cols[8], y=cols[1], color=color2, s=2, ax=ax[8,8])
    data.plot.scatter(x=cols[9], y=cols[1], color=color3, s=2, ax=ax[8,9])

    # Co4O4
    rec3 = ax[7,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec3.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[2], color=color1, s=2, ax=ax[7,0])
    data.plot.scatter(x=cols[1], y=cols[2], color=color1, s=2, ax=ax[7,1])
    data.plot.scatter(x=cols[2], y=cols[2], color=color1, s=2, ax=ax[7,2], zorder=2)
    data.plot.scatter(x=cols[3], y=cols[2], color=color2, s=2, ax=ax[7,3])
    data.plot.scatter(x=cols[4], y=cols[2], color=color2, s=2, ax=ax[7,4])
    data.plot.scatter(x=cols[5], y=cols[2], color=color2, s=2, ax=ax[7,5])
    data.plot.scatter(x=cols[6], y=cols[2], color=color2, s=2, ax=ax[7,6])
    data.plot.scatter(x=cols[7], y=cols[2], color=color2, s=2, ax=ax[7,7])
    data.plot.scatter(x=cols[8], y=cols[2], color=color2, s=2, ax=ax[7,8])
    data.plot.scatter(x=cols[9], y=cols[2], color=color3, s=2, ax=ax[7,9])

    # Ce
    rec4 = ax[6,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec4.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[3], color=color2, s=2, ax=ax[6,0])
    data.plot.scatter(x=cols[1], y=cols[3], color=color2, s=2, ax=ax[6,1])
    data.plot.scatter(x=cols[2], y=cols[3], color=color2, s=2, ax=ax[6,2])
    data.plot.scatter(x=cols[3], y=cols[3], color=color4, s=2, ax=ax[6,3], zorder=2)
    ax[6,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[6,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[6], y=cols[3], color=color4, s=2, ax=ax[6,6])
    data.plot.scatter(x=cols[7], y=cols[3], color=color4, s=2, ax=ax[6,7])
    data.plot.scatter(x=cols[8], y=cols[3], color=color4, s=2, ax=ax[6,8])
    data.plot.scatter(x=cols[9], y=cols[3], color=color5, s=2, ax=ax[6,9])

    # Dy
    rec5 = ax[5,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec5.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[4], color=color2, s=2, ax=ax[5,0])
    data.plot.scatter(x=cols[1], y=cols[4], color=color2, s=2, ax=ax[5,1])
    data.plot.scatter(x=cols[2], y=cols[4], color=color2, s=2, ax=ax[5,2])
    ax[5,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[4], y=cols[4], color=color4, s=2, ax=ax[5,4], zorder=2)
    ax[5,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[6], y=cols[4], color=color4, s=2, ax=ax[5,6])
    data.plot.scatter(x=cols[7], y=cols[4], color=color4, s=2, ax=ax[5,7])
    data.plot.scatter(x=cols[8], y=cols[4], color=color4, s=2, ax=ax[5,8])
    data.plot.scatter(x=cols[9], y=cols[4], color=color5, s=2, ax=ax[5,9])

    # Yb
    rec6 = ax[4,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec6.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[5], color=color2, s=2, ax=ax[4,0])
    data.plot.scatter(x=cols[1], y=cols[5], color=color2, s=2, ax=ax[4,1])
    data.plot.scatter(x=cols[2], y=cols[5], color=color2, s=2, ax=ax[4,2])
    ax[4,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[4,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[5], y=cols[5], color=color4, s=2, ax=ax[4,5], zorder=2)
    data.plot.scatter(x=cols[6], y=cols[5], color=color4, s=2, ax=ax[4,6])
    data.plot.scatter(x=cols[7], y=cols[5], color=color4, s=2, ax=ax[4,7])
    data.plot.scatter(x=cols[8], y=cols[5], color=color4, s=2, ax=ax[4,8])
    data.plot.scatter(x=cols[9], y=cols[5], color=color5, s=2, ax=ax[4,9])

    # OA
    rec7 = ax[3,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec7.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[6], color=color2, s=2, ax=ax[3,0])
    data.plot.scatter(x=cols[1], y=cols[6], color=color2, s=2, ax=ax[3,1])
    data.plot.scatter(x=cols[2], y=cols[6], color=color2, s=2, ax=ax[3,2])
    data.plot.scatter(x=cols[3], y=cols[6], color=color4, s=2, ax=ax[3,3])
    data.plot.scatter(x=cols[4], y=cols[6], color=color4, s=2, ax=ax[3,4])
    data.plot.scatter(x=cols[5], y=cols[6], color=color4, s=2, ax=ax[3,5])
    data.plot.scatter(x=cols[6], y=cols[6], color=color4, s=2, ax=ax[3,6], zorder=2)
    ax[3,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[3,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[9], y=cols[6], color=color5, s=2, ax=ax[3,9])

    # SA
    rec8 = ax[2,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec8.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[7], color=color2, s=2, ax=ax[2,0])
    data.plot.scatter(x=cols[1], y=cols[7], color=color2, s=2, ax=ax[2,1])
    data.plot.scatter(x=cols[2], y=cols[7], color=color2, s=2, ax=ax[2,2])
    data.plot.scatter(x=cols[3], y=cols[7], color=color4, s=2, ax=ax[2,3])
    data.plot.scatter(x=cols[4], y=cols[7], color=color4, s=2, ax=ax[2,4])
    data.plot.scatter(x=cols[5], y=cols[7], color=color4, s=2, ax=ax[2,5])
    ax[2,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[7], y=cols[7], color=color4, s=2, ax=ax[2,7], zorder=2)
    ax[2,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[9], y=cols[7], color=color5, s=2, ax=ax[2,9])

    # TMTACN
    rec9 = ax[1,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec9.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[8], color=color2, s=2, ax=ax[1,0])
    data.plot.scatter(x=cols[1], y=cols[8], color=color2, s=2, ax=ax[1,1])
    data.plot.scatter(x=cols[2], y=cols[8], color=color2, s=2, ax=ax[1,2])
    data.plot.scatter(x=cols[3], y=cols[8], color=color4, s=2, ax=ax[1,3])
    data.plot.scatter(x=cols[4], y=cols[8], color=color4, s=2, ax=ax[1,4])
    data.plot.scatter(x=cols[5], y=cols[8], color=color4, s=2, ax=ax[1,5])
    ax[1,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[1,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[8], y=cols[8], color=color4, s=2, ax=ax[1,8], zorder=2)
    data.plot.scatter(x=cols[9], y=cols[8], color=color5, s=2, ax=ax[1,9])

    # MeOH
    rec10 = ax[0,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_solv, range_solv, color='lightgrey'))
    rec10.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[9], color=color3, s=2, ax=ax[0,0])
    data.plot.scatter(x=cols[1], y=cols[9], color=color3, s=2, ax=ax[0,1])
    data.plot.scatter(x=cols[2], y=cols[9], color=color3, s=2, ax=ax[0,2])
    data.plot.scatter(x=cols[3], y=cols[9], color=color5, s=2, ax=ax[0,3])
    data.plot.scatter(x=cols[4], y=cols[9], color=color5, s=2, ax=ax[0,4])
    data.plot.scatter(x=cols[5], y=cols[9], color=color5, s=2, ax=ax[0,5])
    data.plot.scatter(x=cols[6], y=cols[9], color=color5, s=2, ax=ax[0,6])
    data.plot.scatter(x=cols[7], y=cols[9], color=color5, s=2, ax=ax[0,7])
    data.plot.scatter(x=cols[8], y=cols[9], color=color5, s=2, ax=ax[0,8])
    data.plot.scatter(x=cols[9], y=cols[9], color=color6, s=2, ax=ax[0,9], zorder=2)

    fig.savefig('output_exploration\\CCE_n_visualisation.png')



def CCE_matrix_concise(
    df_ratios, prev_res=None, drop_metadata:bool=True
):

    cols_cluster = ['Co3O(OH)','Co3O','Co4O4','V_TOT']
    data = [df_ratios['Co3O(OH)'],df_ratios['Co3O'],df_ratios['Co4O4'],df_ratios['V_TOT']]
    data = pd.DataFrame(data)
    data = data.T
    data.columns=['Co3O(OH)','Co3O','Co4O4','V_TOT']

    fig, ax = plt.subplots(4,4, figsize=(6.6,6.6), sharex='col', sharey='row')

    color1 = 'indigo'       # cluster + cluster
    color2 = 'darkviolet'   # cluster + template
    color3 = 'violet'       # cluster + solvent
    color4 = 'blue'         # template + template
    color5= 'DarkCyan'      # solvent + template
    color6 = 'g'            # solvent + solvent

    ####################################################################################

    if prev_res != None:
        previous_results_ratio = pd.read_csv(prev_res, index_col=0)
        if drop_metadata == True:
            data_prev = previous_results_ratio.drop(["TYPE"],axis=0)
        else:
            data_prev = previous_results_ratio

        column_floats = ['Co3O(OH)','Co3O','Co4O4','rT1','rT2','V_TOT']

        for i in column_floats:
            data_prev[i] = data_prev[i].astype(float)

        data_prev.replace(0.00, np.nan, inplace=True)
        cols_prev = data_prev.columns

        color_prev = 'grey'

        # Co3O(OH)
        rec1 = ax[3,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec1.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[0], color=color_prev, s=1, ax=ax[3,0], zorder=2)
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[0], color=color_prev, s=1, ax=ax[3,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[0], color=color_prev, s=1, ax=ax[3,2])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[0], color=color_prev, s=1, ax=ax[3,3])

        # Co3O
        rec2 = ax[2,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec2.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[1], color=color_prev, s=1, ax=ax[2,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[1], color=color_prev, s=1, ax=ax[2,1], zorder=2)
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[1], color=color_prev, s=1, ax=ax[2,2])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[1], color=color_prev, s=1, ax=ax[2,3])

        # Co4O4
        rec3 = ax[1,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec3.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[2], color=color_prev, s=1, ax=ax[1,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[2], color=color_prev, s=1, ax=ax[1,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[2], color=color_prev, s=1, ax=ax[1,2], zorder=2)
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[2], color=color_prev, s=1, ax=ax[1,3])

        # V_TOT
        rec4 = ax[0,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 7.6, 7.6, color='lightgrey'))
        rec4.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[7], color=color_prev, s=1, ax=ax[0,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[7], color=color_prev, s=1, ax=ax[0,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[7], color=color_prev, s=1, ax=ax[0,2])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[7], color=color_prev, s=1, ax=ax[0,3], zorder=2)

    ####################################################################################

    # Co3O(OH)
    rec1 = ax[3,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec1.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[0], color=color1, s=2, ax=ax[3,0], zorder=2)
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[0], color=color1, s=2, ax=ax[3,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[0], color=color1, s=2, ax=ax[3,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[0], color=color3, s=2, ax=ax[3,3])

    # Co3O
    rec2 = ax[2,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec2.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[1], color=color1, s=2, ax=ax[2,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[1], color=color1, s=2, ax=ax[2,1], zorder=2)
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[1], color=color1, s=2, ax=ax[2,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[1], color=color3, s=2, ax=ax[2,3])

    # Co4O4
    rec3 = ax[1,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec3.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[2], color=color1, s=2, ax=ax[1,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[2], color=color1, s=2, ax=ax[1,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[2], color=color1, s=2, ax=ax[1,2], zorder=2)
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[2], color=color3, s=2, ax=ax[1,3])

    # V_TOT
    rec4 = ax[0,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 7.6, 7.6, color='lightgrey'))
    rec4.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[3], color=color3, s=2, ax=ax[0,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[3], color=color3, s=2, ax=ax[0,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[3], color=color3, s=2, ax=ax[0,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[3], color=color3, s=2, ax=ax[0,3], zorder=2)

    fig.savefig('output_exploration\\CCE_n_visualisation_ratio.png')



def ISE_matrix(
    df_volumes, prev_res=None, 
    range_cluster:float=4.7, range_template:float=2.6, range_solv:float=6.0
):
    '''
    '''
    data = df_volumes
    cols = data.columns

    fig, ax = plt.subplots(11,11, figsize=(12.1,12.1), sharex='col', sharey='row')

    color1 = 'indigo'       # cluster + cluster
    color2 = 'darkviolet'   # cluster + template
    color3 = 'violet'       # cluster + solvent
    color4 = 'blue'         # template + template
    color5= 'DarkCyan'      # solvent + template
    color6 = 'g'            # solvent + solvent

    ####################################################################################
    
    if prev_res != None:
        previous_results_legible = pd.read_csv(prev_res, index_col=0)
        data_prev = previous_results_legible
        data_prev = data_prev.drop(["prime1","prime2"],axis=0)

        data_prev.replace(0.00, np.nan, inplace=True)
        cols_prev = data_prev.columns

        color_prev = 'grey'
        
        # Cr3O
        rec1 = ax[10,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec1.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,0], zorder=2)
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[0], color=color_prev, s=1, ax=ax[10,10])

        # Mn3O
        rec2 = ax[9,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec2.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,1], zorder=2)
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[1], color=color_prev, s=1, ax=ax[9,10])

        # Fe3O
        rec3 = ax[8,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec3.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,2], zorder=2)
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[2], color=color_prev, s=1, ax=ax[8,10])

        # Co3O
        rec3 = ax[7,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
        rec3.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,3], zorder=2)
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[3], color=color_prev, s=1, ax=ax[7,10])

        # Ce
        rec4 = ax[6,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec4.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,3])
        ax[4,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[4,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,4], zorder=2)
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[4], color=color_prev, s=1, ax=ax[6,10])

        # Dy
        rec5 = ax[5,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec5.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,2])
        ax[5,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,4])
        ax[5,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,5], zorder=2)
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[5], color=color_prev, s=1, ax=ax[5,10])

        # Yb
        rec6 = ax[4,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec6.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,2])
        ax[6,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[6,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,6], zorder=2)
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[6], color=color_prev, s=1, ax=ax[4,10])

        # OA
        rec7 = ax[3,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec7.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,6])
        ax[3,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[3,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,7], zorder=2)
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[7], color=color_prev, s=1, ax=ax[3,10])

        # SA
        rec8 = ax[2,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec8.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,5])
        ax[2,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,7])
        ax[2,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,8], zorder=2)
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[8], color=color_prev, s=1, ax=ax[2,10])

        # TMTACN
        rec9 = ax[1,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
        rec9.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,5])
        ax[1,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        ax[1,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,9], zorder=2)
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[9], color=color_prev, s=1, ax=ax[1,10])

        # MeOH
        rec10 = ax[0,10].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_solv, range_solv, color='lightgrey'))
        rec10.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,3])
        data_prev.plot.scatter(x=cols_prev[4], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,4])
        data_prev.plot.scatter(x=cols_prev[5], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,5])
        data_prev.plot.scatter(x=cols_prev[6], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,6])
        data_prev.plot.scatter(x=cols_prev[7], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,7])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,8])
        data_prev.plot.scatter(x=cols_prev[9], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,9])
        data_prev.plot.scatter(x=cols_prev[10], y=cols_prev[10], color=color_prev, s=1, ax=ax[0,10], zorder=2)


    ####################################################################################

    # Cr3O
    rec1 = ax[10,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec1.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[0], color=color1, s=2, ax=ax[10,0], zorder=2)
    data.plot.scatter(x=cols[1], y=cols[0], color=color1, s=2, ax=ax[10,1])
    data.plot.scatter(x=cols[2], y=cols[0], color=color1, s=2, ax=ax[10,2])
    data.plot.scatter(x=cols[3], y=cols[0], color=color1, s=2, ax=ax[10,3])
    data.plot.scatter(x=cols[4], y=cols[0], color=color2, s=2, ax=ax[10,4])
    data.plot.scatter(x=cols[5], y=cols[0], color=color2, s=2, ax=ax[10,5])
    data.plot.scatter(x=cols[6], y=cols[0], color=color2, s=2, ax=ax[10,6])
    data.plot.scatter(x=cols[7], y=cols[0], color=color2, s=2, ax=ax[10,7])
    data.plot.scatter(x=cols[8], y=cols[0], color=color2, s=2, ax=ax[10,8])
    data.plot.scatter(x=cols[9], y=cols[0], color=color2, s=2, ax=ax[10,9])
    data.plot.scatter(x=cols[10], y=cols[0], color=color3, s=2, ax=ax[10,10])

    # Mn3O
    rec2 = ax[9,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec2.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[1], color=color1, s=2, ax=ax[9,0])
    data.plot.scatter(x=cols[1], y=cols[1], color=color1, s=2, ax=ax[9,1], zorder=2)
    data.plot.scatter(x=cols[2], y=cols[1], color=color1, s=2, ax=ax[9,2])
    data.plot.scatter(x=cols[3], y=cols[1], color=color1, s=2, ax=ax[9,3])
    data.plot.scatter(x=cols[4], y=cols[1], color=color2, s=2, ax=ax[9,4])
    data.plot.scatter(x=cols[5], y=cols[1], color=color2, s=2, ax=ax[9,5])
    data.plot.scatter(x=cols[6], y=cols[1], color=color2, s=2, ax=ax[9,6])
    data.plot.scatter(x=cols[7], y=cols[1], color=color2, s=2, ax=ax[9,7])
    data.plot.scatter(x=cols[8], y=cols[1], color=color2, s=2, ax=ax[9,8])
    data.plot.scatter(x=cols[9], y=cols[1], color=color2, s=2, ax=ax[9,9])
    data.plot.scatter(x=cols[10], y=cols[1], color=color3, s=2, ax=ax[9,10])

    # Fe3O
    rec3 = ax[8,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec3.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[2], color=color1, s=2, ax=ax[8,0])
    data.plot.scatter(x=cols[1], y=cols[2], color=color1, s=2, ax=ax[8,1])
    data.plot.scatter(x=cols[2], y=cols[2], color=color1, s=2, ax=ax[8,2], zorder=2)
    data.plot.scatter(x=cols[3], y=cols[2], color=color1, s=2, ax=ax[8,3])
    data.plot.scatter(x=cols[4], y=cols[2], color=color2, s=2, ax=ax[8,4])
    data.plot.scatter(x=cols[5], y=cols[2], color=color2, s=2, ax=ax[8,5])
    data.plot.scatter(x=cols[6], y=cols[2], color=color2, s=2, ax=ax[8,6])
    data.plot.scatter(x=cols[7], y=cols[2], color=color2, s=2, ax=ax[8,7])
    data.plot.scatter(x=cols[8], y=cols[2], color=color2, s=2, ax=ax[8,8])
    data.plot.scatter(x=cols[9], y=cols[2], color=color2, s=2, ax=ax[8,9])
    data.plot.scatter(x=cols[10], y=cols[2], color=color3, s=2, ax=ax[8,10])

    # Co3O
    rec3 = ax[7,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_cluster, range_cluster, color='lightgrey'))
    rec3.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[3], color=color1, s=2, ax=ax[7,0])
    data.plot.scatter(x=cols[1], y=cols[3], color=color1, s=2, ax=ax[7,1])
    data.plot.scatter(x=cols[2], y=cols[3], color=color1, s=2, ax=ax[7,2])
    data.plot.scatter(x=cols[3], y=cols[3], color=color1, s=2, ax=ax[7,3], zorder=2)
    data.plot.scatter(x=cols[4], y=cols[3], color=color2, s=2, ax=ax[7,4])
    data.plot.scatter(x=cols[5], y=cols[3], color=color2, s=2, ax=ax[7,5])
    data.plot.scatter(x=cols[6], y=cols[3], color=color2, s=2, ax=ax[7,6])
    data.plot.scatter(x=cols[7], y=cols[3], color=color2, s=2, ax=ax[7,7])
    data.plot.scatter(x=cols[8], y=cols[3], color=color2, s=2, ax=ax[7,8])
    data.plot.scatter(x=cols[9], y=cols[3], color=color2, s=2, ax=ax[7,9])
    data.plot.scatter(x=cols[10], y=cols[3], color=color3, s=2, ax=ax[7,10])

    # Ce
    rec4 = ax[6,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec4.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[4], color=color2, s=2, ax=ax[6,0])
    data.plot.scatter(x=cols[1], y=cols[4], color=color2, s=2, ax=ax[6,1])
    data.plot.scatter(x=cols[2], y=cols[4], color=color2, s=2, ax=ax[6,2])
    data.plot.scatter(x=cols[3], y=cols[4], color=color2, s=2, ax=ax[6,3])
    ax[4,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[4,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[4], y=cols[4], color=color4, s=2, ax=ax[6,4], zorder=2)
    data.plot.scatter(x=cols[5], y=cols[4], color=color4, s=2, ax=ax[6,5])
    data.plot.scatter(x=cols[6], y=cols[4], color=color4, s=2, ax=ax[6,6])
    data.plot.scatter(x=cols[7], y=cols[4], color=color4, s=2, ax=ax[6,7])
    data.plot.scatter(x=cols[8], y=cols[4], color=color4, s=2, ax=ax[6,8])
    data.plot.scatter(x=cols[9], y=cols[4], color=color4, s=2, ax=ax[6,9])
    data.plot.scatter(x=cols[10], y=cols[4], color=color5, s=2, ax=ax[6,10])

    # Dy
    rec5 = ax[5,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec5.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[5], color=color2, s=2, ax=ax[5,0])
    data.plot.scatter(x=cols[1], y=cols[5], color=color2, s=2, ax=ax[5,1])
    data.plot.scatter(x=cols[2], y=cols[5], color=color2, s=2, ax=ax[5,2])
    ax[5,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[3], y=cols[5], color=color2, s=2, ax=ax[5,3])
    data.plot.scatter(x=cols[4], y=cols[5], color=color4, s=2, ax=ax[5,4])
    ax[5,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[5], y=cols[5], color=color4, s=2, ax=ax[5,5], zorder=2)
    data.plot.scatter(x=cols[6], y=cols[5], color=color4, s=2, ax=ax[5,6])
    data.plot.scatter(x=cols[7], y=cols[5], color=color4, s=2, ax=ax[5,7])
    data.plot.scatter(x=cols[8], y=cols[5], color=color4, s=2, ax=ax[5,8])
    data.plot.scatter(x=cols[9], y=cols[5], color=color4, s=2, ax=ax[5,9])
    data.plot.scatter(x=cols[10], y=cols[5], color=color5, s=2, ax=ax[5,10])

    # Yb
    rec6 = ax[4,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec6.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[6], color=color2, s=2, ax=ax[4,0])
    data.plot.scatter(x=cols[1], y=cols[6], color=color2, s=2, ax=ax[4,1])
    data.plot.scatter(x=cols[2], y=cols[6], color=color2, s=2, ax=ax[4,2])
    ax[6,5].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[6,6].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[3], y=cols[6], color=color2, s=2, ax=ax[4,3])
    data.plot.scatter(x=cols[4], y=cols[6], color=color4, s=2, ax=ax[4,4])
    data.plot.scatter(x=cols[5], y=cols[6], color=color4, s=2, ax=ax[4,5])
    data.plot.scatter(x=cols[6], y=cols[6], color=color4, s=2, ax=ax[4,6], zorder=2)
    data.plot.scatter(x=cols[7], y=cols[6], color=color4, s=2, ax=ax[4,7])
    data.plot.scatter(x=cols[8], y=cols[6], color=color4, s=2, ax=ax[4,8])
    data.plot.scatter(x=cols[9], y=cols[6], color=color4, s=2, ax=ax[4,9])
    data.plot.scatter(x=cols[10], y=cols[6], color=color5, s=2, ax=ax[4,10])

    # OA
    rec7 = ax[3,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec7.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[7], color=color2, s=2, ax=ax[3,0])
    data.plot.scatter(x=cols[1], y=cols[7], color=color2, s=2, ax=ax[3,1])
    data.plot.scatter(x=cols[2], y=cols[7], color=color2, s=2, ax=ax[3,2])
    data.plot.scatter(x=cols[3], y=cols[7], color=color2, s=2, ax=ax[3,3])
    data.plot.scatter(x=cols[4], y=cols[7], color=color4, s=2, ax=ax[3,4])
    data.plot.scatter(x=cols[5], y=cols[7], color=color4, s=2, ax=ax[3,5])
    data.plot.scatter(x=cols[6], y=cols[7], color=color4, s=2, ax=ax[3,6])
    ax[3,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[3,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[7], y=cols[7], color=color4, s=2, ax=ax[3,7], zorder=2)
    data.plot.scatter(x=cols[8], y=cols[7], color=color4, s=2, ax=ax[3,8])
    data.plot.scatter(x=cols[9], y=cols[7], color=color4, s=2, ax=ax[3,9])
    data.plot.scatter(x=cols[10], y=cols[7], color=color5, s=2, ax=ax[3,10])

    # SA
    rec8 = ax[2,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec8.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[8], color=color2, s=2, ax=ax[2,0])
    data.plot.scatter(x=cols[1], y=cols[8], color=color2, s=2, ax=ax[2,1])
    data.plot.scatter(x=cols[2], y=cols[8], color=color2, s=2, ax=ax[2,2])
    data.plot.scatter(x=cols[3], y=cols[8], color=color2, s=2, ax=ax[2,3])
    data.plot.scatter(x=cols[4], y=cols[8], color=color4, s=2, ax=ax[2,4])
    data.plot.scatter(x=cols[5], y=cols[8], color=color4, s=2, ax=ax[2,5])
    ax[2,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[6], y=cols[8], color=color4, s=2, ax=ax[2,6])
    data.plot.scatter(x=cols[7], y=cols[8], color=color4, s=2, ax=ax[2,7])
    ax[2,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[8], y=cols[8], color=color4, s=2, ax=ax[2,8], zorder=2)
    data.plot.scatter(x=cols[9], y=cols[8], color=color4, s=2, ax=ax[2,9])
    data.plot.scatter(x=cols[10], y=cols[8], color=color5, s=2, ax=ax[2,10])

    # TMTACN
    rec9 = ax[1,9].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='lightgrey'))
    rec9.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[9], color=color2, s=2, ax=ax[1,0])
    data.plot.scatter(x=cols[1], y=cols[9], color=color2, s=2, ax=ax[1,1])
    data.plot.scatter(x=cols[2], y=cols[9], color=color2, s=2, ax=ax[1,2])
    data.plot.scatter(x=cols[3], y=cols[9], color=color2, s=2, ax=ax[1,3])
    data.plot.scatter(x=cols[4], y=cols[9], color=color4, s=2, ax=ax[1,4])
    data.plot.scatter(x=cols[5], y=cols[9], color=color4, s=2, ax=ax[1,5])
    ax[1,7].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    ax[1,8].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_template, range_template, color='grey'))
    data.plot.scatter(x=cols[6], y=cols[9], color=color4, s=2, ax=ax[1,6])
    data.plot.scatter(x=cols[7], y=cols[9], color=color4, s=2, ax=ax[1,7])
    data.plot.scatter(x=cols[8], y=cols[9], color=color4, s=2, ax=ax[1,8])
    data.plot.scatter(x=cols[9], y=cols[9], color=color4, s=2, ax=ax[1,9], zorder=2)
    data.plot.scatter(x=cols[10], y=cols[9], color=color5, s=2, ax=ax[1,10])

    # MeOH
    rec10 = ax[0,10].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), range_solv, range_solv, color='lightgrey'))
    rec10.set_zorder=1
    data.plot.scatter(x=cols[0], y=cols[10], color=color3, s=2, ax=ax[0,0])
    data.plot.scatter(x=cols[1], y=cols[10], color=color3, s=2, ax=ax[0,1])
    data.plot.scatter(x=cols[2], y=cols[10], color=color3, s=2, ax=ax[0,2])
    data.plot.scatter(x=cols[3], y=cols[10], color=color3, s=2, ax=ax[0,3])
    data.plot.scatter(x=cols[4], y=cols[10], color=color5, s=2, ax=ax[0,4])
    data.plot.scatter(x=cols[5], y=cols[10], color=color5, s=2, ax=ax[0,5])
    data.plot.scatter(x=cols[6], y=cols[10], color=color5, s=2, ax=ax[0,6])
    data.plot.scatter(x=cols[7], y=cols[10], color=color5, s=2, ax=ax[0,7])
    data.plot.scatter(x=cols[8], y=cols[10], color=color5, s=2, ax=ax[0,8])
    data.plot.scatter(x=cols[9], y=cols[10], color=color5, s=2, ax=ax[0,9])
    data.plot.scatter(x=cols[10], y=cols[10], color=color6, s=2, ax=ax[0,10], zorder=2)

    fig.savefig('output_exploration\\ISE_n_visualisation.png')



def ISE_matrix_concise(
    df_ratios, prev_res_ratio=None, drop_metadata:bool=True
):
    cols_cluster = ['Cr3O','Mn3O','Fe3O','Co3O','V_TOT']
    data = [df_ratios['Cr3O'],df_ratios['Mn3O'],df_ratios['Fe3O'],df_ratios['Co3O'],df_ratios['V_TOT']]
    data = pd.DataFrame(data)
    data = data.T
    data.columns=['Cr3O','Mn3O','Fe3O','Co3O','V_TOT']

    fig, ax = plt.subplots(5,5, figsize=(8.25,8.25), sharex='col', sharey='row')

    color1 = 'indigo'       # cluster + cluster
    color2 = 'darkviolet'   # cluster + template
    color3 = 'violet'       # cluster + solvent
    color4 = 'blue'         # template + template
    color5= 'DarkCyan'      # solvent + template
    color6 = 'g'            # solvent + solvent

    ####################################################################################

    if prev_res_ratio != None:
        previous_results_ratio = pd.read_csv(prev_res_ratio, index_col=0)
        if drop_metadata == True:
            data_prev = previous_results_ratio.drop(["TYPE"],axis=0)
        else:
            data_prev = previous_results_ratio

        column_floats = ['Cr3O','Mn3O','Fe3O','Co3O','rT1','rT2','V_TOT']

        for i in column_floats:
            data_prev[i] = data_prev[i].astype(float)

        data_prev.replace(0.00, np.nan, inplace=True)
        cols_prev = data_prev.columns

        color_prev = 'grey'

        # Cr3O
        rec1 = ax[4,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec1.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[0], color=color_prev, s=1, ax=ax[4,0], zorder=2)
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[0], color=color_prev, s=1, ax=ax[4,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[0], color=color_prev, s=1, ax=ax[4,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[0], color=color_prev, s=1, ax=ax[4,3])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[0], color=color_prev, s=1, ax=ax[4,4])

        # Mn3O
        rec2 = ax[3,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec2.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[1], color=color_prev, s=1, ax=ax[3,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[1], color=color_prev, s=1, ax=ax[3,1], zorder=2)
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[1], color=color_prev, s=1, ax=ax[3,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[1], color=color_prev, s=1, ax=ax[3,3])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[1], color=color_prev, s=1, ax=ax[3,4])

        # Fe3O
        rec3 = ax[2,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec3.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[2], color=color_prev, s=1, ax=ax[2,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[2], color=color_prev, s=1, ax=ax[2,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[2], color=color_prev, s=1, ax=ax[2,2], zorder=2)
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[2], color=color_prev, s=1, ax=ax[2,3])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[2], color=color_prev, s=1, ax=ax[2,4])

        # Co3O
        rec4 = ax[1,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
        rec4.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[3], color=color_prev, s=1, ax=ax[1,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[3], color=color_prev, s=1, ax=ax[1,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[3], color=color_prev, s=1, ax=ax[1,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[3], color=color_prev, s=1, ax=ax[1,3], zorder=2)
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[3], color=color_prev, s=1, ax=ax[1,4])

        # V_TOT
        rec5 = ax[0,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 7.6, 7.6, color='lightgrey'))
        rec5.set_zorder=1
        data_prev.plot.scatter(x=cols_prev[0], y=cols_prev[8], color=color_prev, s=1, ax=ax[0,0])
        data_prev.plot.scatter(x=cols_prev[1], y=cols_prev[8], color=color_prev, s=1, ax=ax[0,1])
        data_prev.plot.scatter(x=cols_prev[2], y=cols_prev[8], color=color_prev, s=1, ax=ax[0,2])
        data_prev.plot.scatter(x=cols_prev[3], y=cols_prev[8], color=color_prev, s=1, ax=ax[0,3])
        data_prev.plot.scatter(x=cols_prev[8], y=cols_prev[8], color=color_prev, s=1, ax=ax[0,4], zorder=2)

    ####################################################################################

    # Cr3O
    rec1 = ax[4,0].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec1.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[0], color=color1, s=2, ax=ax[4,0], zorder=2)
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[0], color=color1, s=2, ax=ax[4,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[0], color=color1, s=2, ax=ax[4,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[0], color=color1, s=2, ax=ax[4,3])
    data.plot.scatter(x=cols_cluster[4], y=cols_cluster[0], color=color3, s=1, ax=ax[4,4])

    # Mn3O
    rec2 = ax[3,1].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec2.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[1], color=color1, s=2, ax=ax[3,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[1], color=color1, s=2, ax=ax[3,1], zorder=2)
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[1], color=color1, s=2, ax=ax[3,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[1], color=color1, s=2, ax=ax[3,3])
    data.plot.scatter(x=cols_cluster[4], y=cols_cluster[1], color=color3, s=1, ax=ax[3,4])

    # Fe3O
    rec3 = ax[2,2].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec3.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[2], color=color1, s=2, ax=ax[2,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[2], color=color1, s=2, ax=ax[2,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[2], color=color1, s=2, ax=ax[2,2], zorder=2)
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[2], color=color1, s=2, ax=ax[2,3])
    data.plot.scatter(x=cols_cluster[4], y=cols_cluster[2], color=color3, s=1, ax=ax[2,4])

    # Co3O
    rec4 = ax[1,3].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 1.2, 1.2, color='lightgrey'))
    rec4.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[3], color=color1, s=2, ax=ax[1,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[3], color=color1, s=2, ax=ax[1,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[3], color=color1, s=2, ax=ax[1,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[3], color=color1, s=2, ax=ax[1,3], zorder=2)
    data.plot.scatter(x=cols_cluster[4], y=cols_cluster[3], color=color3, s=1, ax=ax[1,4])

    # V_TOT
    rec5 = ax[0,4].add_patch(matplotlib.patches.Rectangle((-0.1,-0.1), 7.6, 7.6, color='lightgrey'))
    rec5.set_zorder=1
    data.plot.scatter(x=cols_cluster[0], y=cols_cluster[4], color=color3, s=2, ax=ax[0,0])
    data.plot.scatter(x=cols_cluster[1], y=cols_cluster[4], color=color3, s=2, ax=ax[0,1])
    data.plot.scatter(x=cols_cluster[2], y=cols_cluster[4], color=color3, s=2, ax=ax[0,2])
    data.plot.scatter(x=cols_cluster[3], y=cols_cluster[4], color=color3, s=2, ax=ax[0,3])
    data.plot.scatter(x=cols_cluster[4], y=cols_cluster[4], color=color3, s=1, ax=ax[0,4], zorder=2)

    fig.savefig('output_exploration\\ISE_n_visualisation_ratio.png')



def exploration_3D(
    df, cols:list, prev_res=None, elev:float=20, azim:float=30
):
    if len(cols) > 3:
        raise ValueError('Only three columns can be plotted.')

    prev_res = pd.read_csv(prev_res, index_col=0)
    prev_res = prev_res.drop(["TYPE"],axis=0)

    fig = plt.figure(figsize=(10,10))

    ax = fig.add_subplot(projection='3d')

    ax.scatter(df[cols[0]], df[cols[1]], df[cols[2]], marker='o', color='r')
    ax.scatter(prev_res[cols[0]], prev_res[cols[1]], prev_res[cols[2]], marker='o', color='r')

    ax.set_xlabel(cols[0])
    ax.set_ylabel(cols[1])
    ax.set_zlabel(cols[2])

    ax.view_init(elev=elev, azim=azim)



def AMF_MS_correlation_plot(df, labels_main, labels_aux, save:str=None):
    fig = plt.figure(figsize=(10,7)) 
    grid = plt.GridSpec(2, 3, wspace=0.4, hspace=0.3)
    ax1 = plt.subplot(grid[:2,:2])
    ax2 = plt.subplot(grid[0,2])
    ax3 = plt.subplot(grid[1,2])

    # Av % Diff vs. Uniqueness
    ax1.set_xlim(-0.05,1.05)
    ax1.set_ylim(-0.05,1.05)
    ax1.set_facecolor('#f2f2f2')
    ax1.set_ylabel('Uniqueness', fontsize=14)
    ax1.set_xlabel('% Diff to Stds', fontsize=14)
    ax1.scatter(df['% Diff Scaled'], df['Uniqueness Scaled'], s=40, c=df['Score'], cmap='cool')

    for i, txt in enumerate(labels_main):
        ax1.annotate(txt, (df['% Diff Scaled'][i]+0.01, df['Uniqueness Scaled'][i]+0.01), fontsize=12)
        

    # Score vs. Av % Diff
    R2_Diff = metrics.r2_score(df['% Diff Scaled'], df['Score'])

    ax2.set_ylabel('Score', fontsize=14)
    ax2.set_xlabel('% Diff to Stds', fontsize=14)
    ax2.scatter(df['% Diff Scaled'], df['Score'], s=25, color='#1085c6')
    ax2.plot([0,1],[0,1], ':', color='#de32ad')

    for i, txt in enumerate(labels_aux):
        ax2.annotate(txt, (df['% Diff Scaled'][i]+0.015, df['Score'][i]+0.01), fontsize=10)

    string = ("R2 = " + str(round(R2_Diff,2)))
    ax2.annotate(string, (0.55, 0.015), fontsize=12, color='#de32ad', 
            bbox=dict(facecolor='none', edgecolor='#de32ad'))


    # Score vs. Uniqueness
    R2_Unique = metrics.r2_score(df['Uniqueness Scaled'], df['Score'])

    ax3.set_ylabel('Score', fontsize=14)
    ax3.set_xlabel('Uniqueness', fontsize=14)
    ax3.scatter(df['Uniqueness Scaled'], df['Score'], s=25, color='#1085c6')
    ax3.plot([0,1],[0,1], ':', color='#de32ad')

    for i, txt in enumerate(labels_aux):
        ax3.annotate(txt, (df['Uniqueness Scaled'][i]+0.015, df['Score'][i]+0.01), fontsize=10)

    string = ("R2 = " + str(round(R2_Unique,2)))
    ax3.annotate(string, (0.55, 0.015), fontsize=12, color='#de32ad', 
            bbox=dict(facecolor='none', edgecolor='#de32ad'))

    plt.tight_layout()

    if save != None:
        plt.savefig(save+'.png')



def AMF_MS_uniqueness_plot(df, title:str="", rule:list=[], color:str='darkgray', xlim:tuple=(200,1000), save:str=None):
    df_to_plot = pd.DataFrame(columns=['m/z','Intensity'])
    df_to_plot['m/z'] = df['index']
    df_to_plot['Intensity'] = 1/df['count']
    unique_obj1 = vis.spec_obj(df_to_plot, 'ms', "Unique Peaks")
    unique_obj1 = vis.core.feature_scale_y(unique_obj1)

    fig,ax = plt.subplots(figsize=(12,4))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0,1.05)
    ax.set_xlim(xlim)
    ax.set_xlabel("m/z", fontsize='x-large')
    ax.set_ylabel("1/Count", fontsize='x-large')
    ax.set_title(
        title, color=color,
        fontsize='x-large', fontweight='bold', loc='right'
    )

    for index,row in unique_obj1.data.iterrows():
        ax.plot([row['m/z'],row['m/z']],[0,row['Intensity']], color) #ed31fd #01bae5 #9062ff

    for i in rule:
        ax.plot([i,i],[0,1], 'k:', linewidth=2.5)
        ax.plot([i],[1], 'ko', linewidth=2.5)

    plt.tight_layout()

    if save != None:
        plt.savefig(save+'.png')



def CCE_tracking_plot(result_data, csv:bool=True, labels:list=[], title:str="", max_n_iters:int=4, pop_size:int=48):
    # import data
    if csv == True:
        df = pd.read_csv(result_data, index_col=0)
    else:
        df = result_data

    # define plot features
    x_range=range(0,len(df))

    plt.rc('axes', labelsize=14)

    fig,ax = plt.subplots(figsize=(12,5))   
    ax.set_xlabel('Experiment #')
    ax.set_ylabel('Combined AMF_Score')
    ax.set_title(title, loc='right', fontsize=14, weight='bold')
    ax.set_ylim(-0.05,1.05)  

    # plot scatter data
    ax.scatter(x=x_range,y=df['AMF'], s=12, c='black')

    # define parameters for coloured bars
    col_height_max = []
    col_height_min = []
    col_average = []

    for n in range(max_n_iters):
        df_slice = df.iloc[n*pop_size:(n+1)*pop_size]
        try: 
            col_height_max.append(max(df_slice['AMF']))
        except ValueError:
            col_height_max.append(0.07)
        try: 
            col_height_min.append(min(df_slice['AMF']))
        except ValueError:
            col_height_min.append(0.0)
        try: 
            col_average.append(sum(df_slice['AMF'])/len(df_slice))
        except (ValueError, ZeroDivisionError):
            col_average.append(0.0) 

    shade = [
        [-0.5, 47.5, 'cyan', 0.1, col_height_max[0], col_height_min[0], col_average[0]],
        [47.5, 95.5, 'darkviolet', 0.1, col_height_max[1], col_height_min[1], col_average[1]],
        [95.5, 143.5, 'dodgerblue', 0.1, col_height_max[2], col_height_min[2], col_average[2]],
        [143.5, 191.5, 'darkviolet', 0.1, col_height_max[3],col_height_min[3],  col_average[3]],
    #    [191.5, 239.5, 'darkviolet', 0.1, col_heights[4]],
    #    [239.5, 287.5, 'dodgerblue', 0.1, col_heights[5]]
    ]

    # plot coloured bars
    y1,y2=0.00,1.05
    for i in shade:
        x = [i[0],i[1]]
        y1s = [y1,y1]
        y2s = [y2,y2]

        if len(i) == 7:
            facecolor=str(i[2])
            alpha=float(i[3])  
            y2s=([float(i[4]),float(i[4])])
            y1s=([float(i[5]),float(i[5])])
            y_av = [i[6],i[6]]
            ax.plot(x, y_av, 'k:')
        if len(i) == 6:
            facecolor=str(i[2])
            alpha=float(i[3])  
            y2s=([float(i[4]),float(i[4])])
            y1s=([float(i[5]),float(i[5])])
        if len(i) == 5:
            facecolor=str(i[2])
            alpha=float(i[3])  
            y2s=([float(i[4]),float(i[4])])
        if len(i) == 4:
            facecolor=str(i[2])
            alpha=float(i[3])
        elif len(i) == 3:
            facecolor=str(i[2])
            alpha=0.15                    
        elif len(i) == 2:
            facecolor='deepskyblue'
            alpha=0.15

        ax.fill_between(x, y1s, y2s, facecolor=facecolor, edgecolor=None, alpha=alpha)
    
    # annotate any points for which labels are requested
    for i,txt in enumerate(labels):
        ax.annotate(txt, (x_range[i], df['AMF'][i]), color='gray')
        
    # annotate each bar with the conditions for that exploration
    ax.annotate('CCE_1 (LHS)', (2,0.02), weight='bold', fontsize=10, color='darkcyan')
    ax.annotate('CCE_2 (GPBO)', (50,0.02), weight='bold', fontsize=10, color='darkviolet')
    ax.annotate('CCE_3 (rand)', (98,0.02), weight='bold', fontsize=10, color='dodgerblue')
    ax.annotate('CCE_4 (GPBO)', (146,0.02), weight='bold', fontsize=10, color='darkviolet')
    #ax.annotate('CCE_5 (rand)', (194,0.02), weight='bold', fontsize=10, color='darkviolet')
    #ax.annotate('CCE_6 (GPBO)', (242,0.02), weight='bold', fontsize=10, color='dodgerblue')



def CCE_boxplot(result_data, csv:bool=True, label_list:list=['1','2','3','4'], title:str="", n_iters:int=4, pop_size:int=48):
    # import data
    if csv == True:
        df = pd.read_csv(result_data, index_col=0)
    else:
        df = result_data

    # set up plot area
    plt.rc('axes', labelsize=14)

    fig,ax = plt.subplots(figsize=(12,5))   
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Combined AMF_Score')
    ax.set_title(title, loc='right', fontsize=14, weight='bold')
    ax.set_ylim(-0.05,1.05) 

    # split data by iteration 
    boxplot_data = []

    for n in range(n_iters):
        df_slice = df.iloc[n*pop_size:(n+1)*pop_size]
        boxplot_data.append(df_slice['AMF'].array)
    
    # plot boxplot
    bp = ax.boxplot(boxplot_data, showfliers=False, patch_artist=True)

    # boxplot formatting
    ax.set_xticklabels(label_list)

    patch_colors = ['#e5ffff', '#f4e5fa', '#e8f3ff', '#f4e5fa']
    median_colors = ['darkcyan', 'darkviolet', 'dodgerblue', 'darkviolet']

    for patch, color in zip(bp['boxes'], patch_colors):
        patch.set_facecolor(color)

    for median, color in zip(bp['medians'], median_colors):
        median.set(color=color, linewidth=3)
        
    for whisker in bp['whiskers']:
        whisker.set(color='gray', linewidth=1.5, linestyle=":")

    for cap in bp['caps']:
        cap.set(color='gray', linewidth=2)



def ISE_tracking_plot(result_data, csv:bool=True, labels:list=[], title:str="", max_n_iters:int=2, pop_size:int=48):
    # import data
    if csv == True:
        df = pd.read_csv(result_data, index_col=0)
    else:
        df = result_data

    # define plot features
    x_range=range(0,len(df))

    plt.rc('axes', labelsize=14)

    fig,ax = plt.subplots(figsize=(12,5))   
    ax.set_xlabel('Experiment #')
    ax.set_ylabel('Combined AMF_Score')
    ax.set_title(title, loc='right', fontsize=14, weight='bold')
    ax.set_ylim(-0.05,1.05)  

    # plot scatter data
    ax.scatter(x=x_range,y=df['AMF'], s=12, c='black')

    # define parameters for coloured bars
    col_height_max = []
    col_height_min = []
    col_average = []

    for n in range(max_n_iters):
        df_slice = df.iloc[n*pop_size:(n+1)*pop_size]
        try: 
            col_height_max.append(max(df_slice['AMF']))
        except ValueError:
            col_height_max.append(0.07)
        try: 
            col_height_min.append(min(df_slice['AMF']))
        except ValueError:
            col_height_min.append(0.0)
        try: 
            col_average.append(sum(df_slice['AMF'])/len(df_slice))
        except (ValueError, ZeroDivisionError):
            col_average.append(0.0) 

    shade = [
        [-0.5, 47.5, 'cyan', 0.1, col_height_max[0], col_height_min[0], col_average[0]],
        [47.5, 95.5, 'darkviolet', 0.1, col_height_max[1], col_height_min[1], col_average[1]],
    #    [95.5, 143.5, 'dodgerblue', 0.1, col_height_max[2], col_height_min[2], col_average[2]],
    #    [143.5, 191.5, 'darkviolet', 0.1, col_height_max[3],col_height_min[3],  col_average[3]],
    #    [191.5, 239.5, 'darkviolet', 0.1, col_heights[4]],
    #    [239.5, 287.5, 'dodgerblue', 0.1, col_heights[5]]
    ]

    # plot coloured bars
    y1,y2=0.00,1.05
    for i in shade:
        x = [i[0],i[1]]
        y1s = [y1,y1]
        y2s = [y2,y2]

        if len(i) == 7:
            facecolor=str(i[2])
            alpha=float(i[3])  
            y2s=([float(i[4]),float(i[4])])
            y1s=([float(i[5]),float(i[5])])
            y_av = [i[6],i[6]]
            ax.plot(x, y_av, 'k:')
        if len(i) == 6:
            facecolor=str(i[2])
            alpha=float(i[3])  
            y2s=([float(i[4]),float(i[4])])
            y1s=([float(i[5]),float(i[5])])
        if len(i) == 5:
            facecolor=str(i[2])
            alpha=float(i[3])  
            y2s=([float(i[4]),float(i[4])])
        if len(i) == 4:
            facecolor=str(i[2])
            alpha=float(i[3])
        elif len(i) == 3:
            facecolor=str(i[2])
            alpha=0.15                    
        elif len(i) == 2:
            facecolor='deepskyblue'
            alpha=0.15

        ax.fill_between(x, y1s, y2s, facecolor=facecolor, edgecolor=None, alpha=alpha)
    
    # annotate any points for which labels are requested
    for i,txt in enumerate(labels):
        ax.annotate(txt, (x_range[i], df['AMF'][i]), color='gray')
        
    # annotate each bar with the conditions for that exploration
    ax.annotate('ISE_1 (LHS)', (2,0.02), weight='bold', fontsize=10, color='darkcyan')
    ax.annotate('ISE_2 (GPBO)', (50,0.02), weight='bold', fontsize=10, color='darkviolet')
    # ax.annotate('CCE_3 (rand)', (98,0.02), weight='bold', fontsize=10, color='dodgerblue')
    # ax.annotate('CCE_4 (GPBO)', (146,0.02), weight='bold', fontsize=10, color='darkviolet')
    #ax.annotate('CCE_5 (rand)', (194,0.02), weight='bold', fontsize=10, color='darkviolet')
    #ax.annotate('CCE_6 (GPBO)', (242,0.02), weight='bold', fontsize=10, color='dodgerblue')



def ISE_boxplot(result_data, csv:bool=True, label_list:list=['1','2'], title:str="", n_iters:int=2, pop_size:int=48):
    # import data
    if csv == True:
        df = pd.read_csv(result_data, index_col=0)
    else:
        df = result_data

    # set up plot area
    plt.rc('axes', labelsize=14)

    fig,ax = plt.subplots(figsize=(12,5))   
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Combined AMF_Score')
    ax.set_title(title, loc='right', fontsize=14, weight='bold')
    ax.set_ylim(-0.05,1.05) 

    # split data by iteration 
    boxplot_data = []

    for n in range(n_iters):
        df_slice = df.iloc[n*pop_size:(n+1)*pop_size]
        boxplot_data.append(df_slice['AMF'].array)
    
    # plot boxplot
    bp = ax.boxplot(boxplot_data, showfliers=False, patch_artist=True)

    # boxplot formatting
    ax.set_xticklabels(label_list)

    patch_colors = ['#e5ffff', '#f4e5fa', '#e8f3ff', '#f4e5fa']
    median_colors = ['darkcyan', 'darkviolet', 'dodgerblue', 'darkviolet']

    for patch, color in zip(bp['boxes'], patch_colors):
        patch.set_facecolor(color)

    for median, color in zip(bp['medians'], median_colors):
        median.set(color=color, linewidth=3)
        
    for whisker in bp['whiskers']:
        whisker.set(color='gray', linewidth=1.5, linestyle=":")

    for cap in bp['caps']:
        cap.set(color='gray', linewidth=2)