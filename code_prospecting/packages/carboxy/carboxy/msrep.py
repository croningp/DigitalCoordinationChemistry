### CODE: CATRIONA MACGREGOR & DANIEL KOWALSKI

from copy import deepcopy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

"""
def gen_table(spectrum1, spectrum2):

    #create deepcopies to avoid editting original spectrum object
    spec1 = deepcopy(spectrum1)
    spec2 = deepcopy(spectrum2)

    #define dataframe as spec1.data, not list as list will be harder to work with later
    #run this and will save having it in each cell
    df = pd.DataFrame(spec1.data)

    #add column
    df["Relative Intensity 2"]=""

    #delete column
    del df['Intensity']

    #rename column
    df.rename(columns ={"Relative Intensity":"Relative Intensity 1"}, inplace=True)

    #round to 0 decimal places
    #every time change dataframe can make df = original line of code
    df = df.round(0)

    #to check if columns were added/deleted/renamed - not necessary
    #can say df or spec1.data but spec1.data won't include decimal places, only changes to columns... why?
    #spec1.data

    #define dataframe2 as spec2.data to save having to do in future cells
    df2 = pd.DataFrame(spec2.data)

    #delete column
    del df2['Intensity']

    #rename column to  Relative Intensity 2
    df2.rename(columns = {"Relative Intensity":"Relative Intensity 2"}, inplace=True)

    #round to 0 decimal places
    df2 = df2.round(0)

    #to check - again with spec2.data it won't show decimals
    #print(spec2.data)

    #make a functioning loop using iterrows where it goes through each m/z and prints
    for index, row in df2.iterrows():
        v = row['m/z']
        x = row['Relative Intensity 2']
        i = None
        
        for index, row in df.iterrows():
            if row['m/z'] == v:
                i = index

        if i != None:
            df.loc[i, 'Relative Intensity 2'] = x
        else:
            values_to_add = {'m/z': v, 'Relative Intensity 2': x}
            row_to_add = pd.Series(values_to_add)
            df = df.append(row_to_add, ignore_index=True)

    # replace field that's entirely space (or empty) with NaN
    df = df.replace(r'^\s*$', np.nan, regex=True)

    # replace 'NaN's with '0.000'
    df = df.replace(np.nan, 0)

    return df
"""



def gen_table(spectrum1, spectrum2):

    #create deepcopies to avoid editting original spectrum object
    spec1 = deepcopy(spectrum1)
    spec2 = deepcopy(spectrum2)

    #define dataframe as spec1.data, not list as list will be harder to work with later
    #run this and will save having it in each cell
    df = pd.DataFrame(spec1.data)
    #add column
    df["Relative Intensity 2"]=""
    #delete column
    del df['Intensity']
    #rename column
    df.rename(columns ={"Relative Intensity":"Relative Intensity 1"}, inplace=True)
    #round to 0 decimal places
    #every time change dataframe can make df = original line of code
    df = df.round(0)
    # create copy of df as df1
    df1 = deepcopy(df)

    # remove duplicates from df1
    # setup empty variables/lists
    prev_mz1 = None
    prev_ri1 = None
    mz1 = []
    ri1 = []

    # search through DataFrame for duplicate entries
    for index, row in df1.iterrows():
        if prev_mz1 == row['m/z']:
            # where a duplicate has a higher value, overwrite with higher value
            if prev_ri1 < row['Relative Intensity 1']:
                ri1[-1] = row['Relative Intensity 1']
        # otherwise just append 
        # (should also append first instance of any duplicate)
        else:
            mz1.append(row['m/z'])
            ri1.append(row['Relative Intensity 1'])

        prev_mz1 = row['m/z']
        prev_ri1 = row['Relative Intensity 1']
        
    
    df1_dict = {'m/z': mz1, 'Relative Intensity 1': ri1}
    df1 = pd.DataFrame(df1_dict)
    
    
    #define dataframe2 as spec2.data to save having to do in future cells
    df2 = pd.DataFrame(spec2.data)
    #delete column
    del df2['Intensity']
    #rename column to  Relative Intensity 2
    df2.rename(columns = {"Relative Intensity":"Relative Intensity 2"}, inplace=True)
    #round to 0 decimal places
    df2 = df2.round(0)

    # remove duplicates from df2
    # setup empty variables/lists
    prev_mz2 = None
    prev_ri2 = None

    mz2 = []
    ri2 = []

    # search through DataFrame for duplicate entries
    for index, row in df2.iterrows():
        if prev_mz2 == row['m/z']:
            # where a duplicate has a higher value, overwrite with higher value
            if prev_ri2 < row['Relative Intensity 2']:
                ri2[-1] = row['Relative Intensity 2']
        # otherwise just append 
        # (should also append first instance of any duplicate)
        else:
            mz2.append(row['m/z'])
            ri2.append(row['Relative Intensity 2'])

        prev_mz2 = row['m/z']
        prev_ri2 = row['Relative Intensity 2']
        
    
    df2_dict = {'m/z': mz2, 'Relative Intensity 2': ri2}
    df2 = pd.DataFrame(df2_dict)
    
    #combination of the two datasets
    for index, row in df2.iterrows():
        # where a value of m/z is present in both dataset add to Relative Intensity 2 column
        v = row['m/z']
        x = row['Relative Intensity 2']
        i = None
        
        for index, row in df1.iterrows():
            if row['m/z'] == v:
                i = index
        
        # where a value of m/z is present only in the second dataset add to bottom of DataFrame
        if i != None:
            df1.loc[i, 'Relative Intensity 2'] = x
        else:
            values_to_add = {'m/z': v, 'Relative Intensity 2': x}
            row_to_add = pd.Series(values_to_add)
            df1 = df1.append(row_to_add, ignore_index=True)

    # replace field that's entirely space (or empty) with NaN
    df = df1.replace(r'^\s*$', np.nan, regex=True)

    # replace 'NaN's with '0.000'
    df = df.replace(np.nan, 0)

    return df



def make_graph(df, label1, label2):
    
    # formatting
    SMALL_SIZE = 10
    MEDIUM_SIZE = 13
    BIGGER_SIZE = 15

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    # sets up plot area
    fig, ax = plt.subplots(figsize=(5,5))
    
    # plots data
    ax.scatter(
        df['Relative Intensity 1'], df['Relative Intensity 2'],
        s=20, color='darkorchid', marker='x'
    )

    tit = str(label1 + " vs. " + label2)
    ax.set(xlabel='Relative Intensity 1', ylabel='Relative Intensity 2', title=tit)
    
    # compare to the optimum result (straight line diagonally)
    ax.plot([0,100],[0,100], ':', color='#DE32AD')

    R2_Diff = metrics.r2_score(df['Relative Intensity 2'], df['Relative Intensity 1'])
    string = ("R2 = " + str(round(R2_Diff,2)))

    ax.annotate(string, (75, 40), fontsize=12, color='#DE32AD',
            bbox=dict(facecolor='none', edgecolor='#DE32AD'))
    
    # saves and displays
    file_name = str("output_reproducibility\\"+label1+"+"+label2+".png")
    fig.savefig(file_name)
    fig.tight_layout()
    plt.show()



def compare_spectra(data_list):
    # create empty list to store results
    df_list = []

    # iterate through data_list, comparing every pair of spectra listed
    for n in range(len(data_list)):
        dataset = data_list.pop(0)

        for i in data_list:
            df_repr = gen_table(dataset, i)
            make_graph(df_repr, dataset.label, i.label)
            df_list.append([dataset.label, i.label, df_repr])
    
    return df_list