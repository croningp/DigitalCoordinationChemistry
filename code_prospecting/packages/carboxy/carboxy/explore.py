### CODE: DANIEL KOWALSKI

from summit.domain import *
from summit.utils.dataset import DataSet
from summit.strategies.random import LHS
from summit.strategies.random import Random
from summit.strategies import SOBO
import pandas as pd
import numpy as np
import skopt
from copy import deepcopy


class CCE():
    def __init__(
        self, expl_name:str="CCE_n", prev_res:str=None, package:str='skopt',
        max_volume:float=7.5, cluster_max:float=1.0, template_max:float=0.5,
        n_initial_points:int=48, kappa:float=3, double_weighting:bool=True
    ) -> None:
        
        self.expl_name = expl_name
        self.prev_res = prev_res
        self.package = package   # must be either 'skopt' or 'summit'
        self.max_volume = max_volume
        self.cluster_max = cluster_max
        self.template_max = template_max
        self.n_initial_points = n_initial_points
        self.kappa = kappa
        self.double_weighting = double_weighting


    def summit_domain(self):
        self.domain = Domain()

        # define cluster space
        self.domain += ContinuousVariable(
            name='Co3O(OH)', 
            description='ratio of {Co3O(OH)} intermediate to other components', 
            bounds=[0,self.cluster_max]
        )

        self.domain += ContinuousVariable(
            name='Co3O', 
            description='ratio of {Co3O} intermediate to other components', 
            bounds=[0,self.cluster_max]
        )

        self.domain += ContinuousVariable(
            name='Co4O4', 
            description='ratio of {Co4} cubane to other components', 
            bounds=[0,self.cluster_max]
        )

        # define template space
        self.domain += CategoricalVariable(
            name='T1', 
            description='choice of template from set 1: {Ce, Dy, Yb, None}', 
            levels=['Ce','Dy','Yb','MeOH']
        )

        self.domain += ContinuousVariable(
            name='rT1', 
            description='ratio of template choice 1 to other components', 
            bounds=[0,self.template_max]
        )

        self.domain += CategoricalVariable(
            name='T2', 
            description='choice of template from set 2: {oxalic acid, succinic acid, TMTACN, None}', 
            levels=['OA','SA','TMTACN','MeOH']
        )

        self.domain += ContinuousVariable(
            name='rT2', 
            description='ratio of template choice 2 to other components', 
            bounds=[0,self.template_max]
        )

        # define total volume
        volume_UB = self.max_volume
        volume_LB = self.max_volume/2

        self.domain += ContinuousVariable(
            name='V_TOT', 
            description='total volume of reagents before dilution in MeOH', 
            bounds=[volume_LB,volume_UB]
        )

        # define response variable
        self.domain += ContinuousVariable(
            name='AMF', 
            description='analytical mapping function', 
            bounds=[0,1], 
            is_objective=True, 
            maximise=False
        )

        return self.domain


    def skopt_optimisation(self) -> None:
        self.dimensions = [
            (0.0, self.cluster_max),      # Co3O(OH)
            (0.0, self.cluster_max),      # Co3O
            (0.0, self.cluster_max),      # Co4O4
            ('Ce','Dy','Yb','MeOH'),      # T1
            (0.0, self.template_max),     # rT1
            ('OA','SA','TMTACN','MeOH'),  # T2
            (0.0, self.template_max),     # rT2
            (3.75, 7.5),                  # V_TOT
        ]

        # define optimisation strategy
        self.optimizer = skopt.Optimizer(
            dimensions=self.dimensions, 
            base_estimator="gp",
            n_initial_points=self.n_initial_points,
            initial_point_generator="LHS",
            acq_func="LCB",
            acq_func_kwargs = {"kappa":self.kappa},
            acq_optimizer="lbfgs",
            random_state=None
        )

        # update algorithm with previous results
        self.previous_results = pd.read_csv(self.prev_res, index_col=0)

        x_list = self.previous_results[self.previous_results.columns[:8]].values.tolist()
        y_list = self.previous_results['AMF'].values.tolist()

        # correction as skopt minimizes
        for n in range(len(y_list)):
            y_list[n] = 1 - y_list[n]

        # add double_weighting to the initial screen if required
        if self.double_weighting == True:
            first_x_iter = x_list[:self.n_initial_points]
            x_list = x_list + first_x_iter

            first_y_iter = y_list[:self.n_initial_points]
            y_list = y_list + first_y_iter

            self.tell = self.optimizer.tell(
                x = x_list,
                y = y_list,
                fit = True
            )

        elif self.double_weighting == False:
            self.tell = self.optimizer.tell(
                x = x_list,
                y = y_list,
                fit = True
            )

        # ask for points
        self.ratio = self.optimizer.ask(self.n_initial_points, strategy='cl_mean')
        self.ratio = pd.DataFrame(self.ratio, columns=['Co3O(OH)','Co3O','Co4O4','T1','rT1','T2','rT2','V_TOT'])

        # print data
        self.ratio.to_csv('output_exploration\\CCE_n_ratio.csv')

        return self.ratio


    def summit_optimisation(self) -> None:
        self.summit_domain()
        # define any datapoints already known
        c1 = []
        c2 = []
        c3 = []
        T1 = []
        rT1 = []
        T2 = []
        rT2 = []
        V_TOT = []
        AMF = []

        # create Summit DataSet 
        values = {
            'Co3O(OH)': c1, 
            'Co3O': c2, 
            'Co4O4': c3,
            'T1': T1,
            'rT1': rT1,
            'T2': T2,
            'rT2': rT2,
            'V_TOT': V_TOT,
            'AMF': AMF
        }

        if self.prev_res == None:
            # define optimisation strategy
            strategy = LHS(self.domain, random_state=np.random.RandomState(3))
            
            # ask for points
            self.ratio = strategy.suggest_experiments(self.n_initial_points)

        else:
            # define optimisation strategy
            strategy = SOBO(self.domain, acquisition_type='LCB')

            # update algorithm with previous results
            previous_results_df = pd.read_csv(self.prev_res, index_col=0)
            previous_results = DataSet.from_df(previous_results_df)

            # ask for points
            self.ratio = strategy.suggest_experiments(self.n_initial_points, previous_results)

        # print data
        self.ratio.to_csv('output_exploration\\CCE_n_ratio.csv')

        return self.ratio
    

    def random_optimisation(self):
        self.summit_domain()

        # define optimisation strategy
        strategy = Random(self.domain, random_state=None)
            
        # ask for points
        self.ratio = strategy.suggest_experiments(self.n_initial_points)

        # print data
        self.ratio.to_csv('output_exploration\\CCE_n_ratio.csv')

        return self.ratio


    def generate_new_values(self):
        if self.package == 'summit':
            df_to_return = self.summit_optimisation()
        elif self.package == 'skopt':
            df_to_return = self.skopt_optimisation()
        elif self.package == 'random':
            df_to_return = self.random_optimisation()

        return df_to_return


    def calculate_volumes(self):
        # unpack one line at a time into dataframe - calculating actual volumes as we go
        self.cols = ['Co3O(OH)','Co3O','Co4O4','Ce','Dy','Yb','OA','SA','TMTACN','MeOH']
        self.df = pd.DataFrame(columns=self.cols)

        for index, row in self.ratio.iterrows():
            if (self.package == 'summit') or (self.package == 'random'):
                # formatting datatypes
                c1_val = float(row['Co3O(OH)'].DATA)
                c2_val = float(row['Co3O'].DATA)
                c3_val = float(row['Co4O4'].DATA)
                T1_val = str(row['T1'].DATA)
                rT1_val = float(row['rT1'].DATA)
                T2_val = str(row['T2'].DATA)
                rT2_val = float(row['rT2'].DATA)
                V_TOT_val = float(row['V_TOT'].DATA)

            if self.package == 'skopt':
                # formatting datatypes
                c1_val = float(row['Co3O(OH)'])
                c2_val = float(row['Co3O'])
                c3_val = float(row['Co4O4'])
                T1_val = str(row['T1'])
                rT1_val = float(row['rT1'])
                T2_val = str(row['T2'])
                rT2_val = float(row['rT2'])
                V_TOT_val = float(row['V_TOT'])
            
            # adjust volumes to match V_TOT
            adjusted_volumes=[]
            ratio = [c1_val,c2_val,c3_val,rT1_val,rT2_val]
            sum_volume = sum(ratio)
            
            for i in ratio:
                adjusted_individual = i * V_TOT_val/sum_volume 
                adjusted_volumes.append(adjusted_individual)
            
            if round(sum(adjusted_volumes)) != round(V_TOT_val):
                print(r"'adjusted_volumes' does not match precisely with 'V_TOT_val' for "+str(index)+". adjusted_volumes=%s, V_TOT_val=%s" %(round(sum(adjusted_volumes)), round(V_TOT_val)))
        
            # reformat into dictionary for addition to DataFrame
            entry={}
            entry['Co3O(OH)'] = adjusted_volumes[0]
            entry['Co3O'] = adjusted_volumes[1]
            entry['Co4O4'] = adjusted_volumes[2]
            extra_MeOH = self.max_volume - V_TOT_val
            if T1_val == T2_val:
                entry['MeOH'] = extra_MeOH + adjusted_volumes[3] + adjusted_volumes[4]
            elif T1_val == 'MeOH':
                entry['MeOH'] = extra_MeOH + adjusted_volumes[3]
                entry[T2_val]=adjusted_volumes[4]
            elif T2_val == 'MeOH':
                entry[T1_val]=adjusted_volumes[3]
                entry['MeOH'] = extra_MeOH + adjusted_volumes[4]
            else:
                entry[T1_val]=adjusted_volumes[3]
                entry[T2_val]=adjusted_volumes[4]
                entry['MeOH'] = extra_MeOH

            # append to df
            self.df = self.df.append(entry, ignore_index=True)
            self.df = self.df.round(2)

        return self.df


    def write_legible_df(self):
        # write .csv suitable for ClusterBot
        df_copy = deepcopy(self.df)
        df_copy = df_copy.fillna(0)

        # create list of rxn id's
        rxnid = []
        for i in range(self.n_initial_points):
            head = self.expl_name
            tail = '_'+str(i)
            string = head+tail
            rxnid.append(string)

        # bind to existing dataframe, reorder, rename and save as .csv
        self.df_rxnid = deepcopy(df_copy)
        self.df_rxnid['Reaction ID'] = rxnid

        rxnid_cols = self.df_rxnid.columns.tolist()
        rxnid_cols = rxnid_cols[-1:] + rxnid_cols[:-1]
        self.df_rxnid = self.df_rxnid[rxnid_cols]

        self.df_rxnid = self.df_rxnid.rename(columns={
            rxnid_cols[0]: rxnid_cols[0],
            rxnid_cols[1]: 'pump1', 
            rxnid_cols[2]: 'pump2',
            rxnid_cols[3]: 'pump3',
            rxnid_cols[4]: 'pump4',
            rxnid_cols[5]: 'pump5',
            rxnid_cols[6]: 'pump6',
            rxnid_cols[7]: 'pump7',
            rxnid_cols[8]: 'pump8',
            rxnid_cols[9]: 'pump9',
            rxnid_cols[10]: 'pump10',
        })

        prime1 = {'Reaction ID':'prime1','pump1':1,'pump2':1,'pump3':1,'pump4':1,'pump5':1,
                'pump6':1,'pump7':1,'pump8':1,'pump9':1,'pump10':1}
        self.df_rxnid = self.df_rxnid.append(prime1, ignore_index=True)

        prime2 = {'Reaction ID':'prime2','pump1':1,'pump2':1,'pump3':1,'pump4':1,'pump5':1,
                'pump6':1,'pump7':1,'pump8':1,'pump9':1,'pump10':1}
        self.df_rxnid = self.df_rxnid.append(prime2, ignore_index=True)

        #self.df_rxnid = self.df_rxnid.reindex([48,49,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
        #                            26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47])
        
        if self.n_initial_points == 48:
            self.df_rxnid = self.df_rxnid.reindex([48,49,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                        26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47])
        elif self.n_initial_points == 24:
            self.df_rxnid = self.df_rxnid.reindex([24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])    
        else:
            raise ValueError('Cannot write legible.csv file with given number of initial points.')   


        self.df_rxnid.to_csv('output_exploration\\CCE_n_legible.csv', index=False)
        self.df_rxnid

        # determine necessary amounts for the run
        for i in self.cols:
            sum_col = sum(df_copy[i])
            print(i, ":", round(sum_col,2), "mL")

        return self.df_rxnid


#####################################################################################


class ISE():
    def __init__(
        self, expl_name:str="ISE_n", prev_res:str=None, package:str='skopt',
        max_volume:float=7.5, cluster_max:float=1.0, template_max:float=0.5,
        n_initial_points:int=48, kappa:float=3, double_weighting:bool=True
    ) -> None:
        
        self.expl_name = expl_name
        self.prev_res = prev_res
        self.package = package   # must be either 'skopt' or 'summit'
        self.max_volume = max_volume
        self.cluster_max = cluster_max
        self.template_max = template_max
        self.n_initial_points = n_initial_points
        self.kappa = kappa
        self.double_weighting = double_weighting


    def summit_domain(self):
        self.domain = Domain()

        # define cluster space
        self.domain += ContinuousVariable(
            name='Cr3O', 
            description='ratio of {Cr3O} intermediate to other components', 
            bounds=[0,self.cluster_max]
        )

        self.domain += ContinuousVariable(
            name='Mn3O', 
            description='ratio of {Mn3O} intermediate to other components', 
            bounds=[0,self.cluster_max]
        )

        self.domain += ContinuousVariable(
            name='Fe3O', 
            description='ratio of {Fe3O} cubane to other components', 
            bounds=[0,self.cluster_max]
        )

        self.domain += ContinuousVariable(
            name='Co3O', 
            description='ratio of {Co3O} cubane to other components', 
            bounds=[0,self.cluster_max]
        )

        # define template space
        self.domain += CategoricalVariable(
            name='T1', 
            description='choice of template from set 1: {Ce, Dy, Yb, None}', 
            levels=['Ce','Dy','Yb','MeOH']
        )

        self.domain += ContinuousVariable(
            name='rT1', 
            description='ratio of template choice 1 to other components', 
            bounds=[0,self.template_max]
        )

        self.domain += CategoricalVariable(
            name='T2', 
            description='choice of template from set 2: {oxalic acid, succinic acid, TMTACN, None}', 
            levels=['OA','SA','TMTACN','MeOH']
        )

        self.domain += ContinuousVariable(
            name='rT2', 
            description='ratio of template choice 2 to other components', 
            bounds=[0,self.template_max]
        )

        # define total volume
        volume_UB = self.max_volume
        volume_LB = self.max_volume/2

        self.domain += ContinuousVariable(
            name='V_TOT', 
            description='total volume of reagents before dilution in MeOH', 
            bounds=[volume_LB,volume_UB]
        )

        # define response variable
        self.domain += ContinuousVariable(
            name='AMF', 
            description='analytical mapping function', 
            bounds=[0,1], 
            is_objective=True, 
            maximise=False
        )

        return self.domain


    def skopt_optimisation(self) -> None:
        self.dimensions = [
            (0.0, self.cluster_max),      # Cr3O
            (0.0, self.cluster_max),      # Mn3O
            (0.0, self.cluster_max),      # Fe3O
            (0.0, self.cluster_max),      # Co3O
            ('Ce','Dy','Yb','MeOH'),      # T1
            (0.0, self.template_max),     # rT1
            ('OA','SA','TMTACN','MeOH'),  # T2
            (0.0, self.template_max),     # rT2
            (3.75, 7.5),                  # V_TOT
        ]

        # define optimisation strategy
        self.optimizer = skopt.Optimizer(
            dimensions=self.dimensions, 
            base_estimator="gp",
            n_initial_points=self.n_initial_points,
            initial_point_generator="LHS",
            acq_func="LCB",
            acq_func_kwargs = {"kappa":self.kappa},
            acq_optimizer="lbfgs",
            random_state=None
        )

        # update algorithm with previous results
        self.previous_results = pd.read_csv(self.prev_res, index_col=0)
        #self.previous_results['V_TOT'] = self.previous_results['V_TOT']*0.75

        x_list = self.previous_results[self.previous_results.columns[:9]].values.tolist()
        y_list = self.previous_results['AMF'].values.tolist()

        # correction as skopt minimizes
        for n in range(len(y_list)):
            y_list[n] = 1 - y_list[n]

        #
        if self.double_weighting == True:
            first_x_iter = x_list[:self.n_initial_points]
            x_list = x_list + first_x_iter

            first_y_iter = y_list[:self.n_initial_points]
            y_list = y_list + first_y_iter

            self.tell = self.optimizer.tell(
                x = x_list,
                y = y_list,
                fit = True
            )

        elif self.double_weighting == False:
            self.tell = self.optimizer.tell(
                x = x_list,
                y = y_list,
                fit = True
            )

        # ask for points
        self.ratio = self.optimizer.ask(self.n_initial_points, strategy='cl_mean')
        self.ratio = pd.DataFrame(self.ratio, columns=['Cr3O','Mn3O','Fe3O','Co3O','T1','rT1','T2','rT2','V_TOT'])

        # print data
        self.ratio.to_csv('output_exploration\\ISE_n_ratio.csv')

        return self.ratio


    def summit_optimisation(self) -> None:
        self.summit_domain()
        # define any datapoints already known
        c1 = []
        c2 = []
        c3 = []
        c4 = []
        T1 = []
        rT1 = []
        T2 = []
        rT2 = []
        V_TOT = []
        AMF = []

        # create Summit DataSet 
        values = {
            'Cr3O': c1, 
            'Mn3O': c2, 
            'Fe3O': c3,
            'Co3O': c4, 
            'T1': T1,
            'rT1': rT1,
            'T2': T2,
            'rT2': rT2,
            'V_TOT': V_TOT,
            'AMF': AMF
        }

        if self.prev_res == None:
            # define optimisation strategy
            strategy = LHS(self.domain, random_state=np.random.RandomState(3))
            
            # ask for points
            self.ratio = strategy.suggest_experiments(self.n_initial_points)

        else:
            # define optimisation strategy
            strategy = SOBO(self.domain, acquisition_type='LCB')

            # update algorithm with previous results
            previous_results_df = pd.read_csv(self.prev_res, index_col=0)
            previous_results = DataSet.from_df(previous_results_df)

            # ask for points
            self.ratio = strategy.suggest_experiments(self.n_initial_points, previous_results)

            # print data
            self.ratio.to_csv('output_exploration\\ISE_n_ratio.csv')

        return self.ratio


    def generate_new_values(self):
        if self.package == 'summit':
            df_to_return = self.summit_optimisation()
        elif self.package == 'skopt':
            df_to_return = self.skopt_optimisation()
        #elif self.package == 'random':
        #    df_to_return = self.random_optimisation()

        return df_to_return


    def calculate_volumes(self):
        # unpack one line at a time into dataframe - calculating actual volumes as we go
        self.cols = ['Cr3O','Mn3O','Fe3O','Co3O','Ce','Dy','Yb','OA','SA','TMTACN','MeOH']
        self.df = pd.DataFrame(columns=self.cols)

        for index, row in self.ratio.iterrows():
            if self.package == 'summit':
                # formatting datatypes
                c1_val = float(row['Cr3O'].DATA)
                c2_val = float(row['Mn3O'].DATA)
                c3_val = float(row['Fe3O'].DATA)
                c4_val = float(row['Co3O'].DATA)
                T1_val = str(row['T1'].DATA)
                rT1_val = float(row['rT1'].DATA)
                T2_val = str(row['T2'].DATA)
                rT2_val = float(row['rT2'].DATA)
                V_TOT_val = float(row['V_TOT'].DATA)

            if self.package == 'skopt':
                # formatting datatypes
                c1_val = float(row['Cr3O'])
                c2_val = float(row['Mn3O'])
                c3_val = float(row['Fe3O'])
                c4_val = float(row['Co3O'])
                T1_val = str(row['T1'])
                rT1_val = float(row['rT1'])
                T2_val = str(row['T2'])
                rT2_val = float(row['rT2'])
                V_TOT_val = float(row['V_TOT'])
            
            # adjust volumes to match V_TOT
            adjusted_volumes=[]
            ratio = [c1_val,c2_val,c3_val,c4_val,rT1_val,rT2_val]
            sum_volume = sum(ratio)
            
            for i in ratio:
                adjusted_individual = i * V_TOT_val/sum_volume 
                adjusted_volumes.append(adjusted_individual)
            
            if round(sum(adjusted_volumes)) != round(V_TOT_val):
                print("'adjusted_volumes' does not match precisely with 'V_TOT_val' for %s. adjusted_volumes=%s, V_TOT_val=%s" % (str(index), round(sum(adjusted_volumes)), round(V_TOT_val)))
        
            # reformat into dictionary for addition to DataFrame
            entry={}
            entry['Cr3O'] = adjusted_volumes[0]
            entry['Mn3O'] = adjusted_volumes[1]
            entry['Fe3O'] = adjusted_volumes[2]
            entry['Co3O'] = adjusted_volumes[3]
            extra_MeOH = self.max_volume - V_TOT_val
            if T1_val == T2_val:
                entry['MeOH'] = extra_MeOH + adjusted_volumes[4] + adjusted_volumes[5]
            elif T1_val == 'MeOH':
                entry['MeOH'] = extra_MeOH + adjusted_volumes[4]
                entry[T2_val]=adjusted_volumes[5]
            elif T2_val == 'MeOH':
                entry[T1_val]=adjusted_volumes[4]
                entry['MeOH'] = extra_MeOH + adjusted_volumes[5]
            else:
                entry[T1_val]=adjusted_volumes[4]
                entry[T2_val]=adjusted_volumes[5]
                entry['MeOH'] = extra_MeOH

            # append to df
            self.df = self.df.append(entry, ignore_index=True)
            self.df = self.df.round(2)

        return self.df


    def write_legible_df(self):
        # write .csv suitable for ClusterBot
        df_copy = deepcopy(self.df)
        df_copy = df_copy.fillna(0)

        # create list of rxn id's
        rxnid = []
        for i in range(self.n_initial_points):
            head = self.expl_name
            tail = '_'+str(i)
            string = head+tail
            rxnid.append(string)

        # bind to existing dataframe, reorder, rename and save as .csv
        self.df_rxnid = deepcopy(df_copy)
        self.df_rxnid['Reaction ID'] = rxnid

        rxnid_cols = self.df_rxnid.columns.tolist()
        rxnid_cols = rxnid_cols[-1:] + rxnid_cols[:-1]
        self.df_rxnid = self.df_rxnid[rxnid_cols]

        self.df_rxnid = self.df_rxnid.rename(columns={
            rxnid_cols[0]: rxnid_cols[0],
            rxnid_cols[1]: 'pump1', 
            rxnid_cols[2]: 'pump2',
            rxnid_cols[3]: 'pump3',
            rxnid_cols[4]: 'pump4',
            rxnid_cols[5]: 'pump5',
            rxnid_cols[6]: 'pump6',
            rxnid_cols[7]: 'pump7',
            rxnid_cols[8]: 'pump8',
            rxnid_cols[9]: 'pump9',
            rxnid_cols[10]: 'pump10',
            rxnid_cols[11]: 'pump11',
        })

        prime1 = {'Reaction ID':'prime1','pump1':1,'pump2':1,'pump3':1,'pump4':1,'pump5':1,
                'pump6':1,'pump7':1,'pump8':1,'pump9':1,'pump10':1,'pump11':1}
        self.df_rxnid = self.df_rxnid.append(prime1, ignore_index=True)

        prime2 = {'Reaction ID':'prime2','pump1':1,'pump2':1,'pump3':1,'pump4':1,'pump5':1,
                'pump6':1,'pump7':1,'pump8':1,'pump9':1,'pump10':1,'pump11':1}
        self.df_rxnid = self.df_rxnid.append(prime2, ignore_index=True)

        if self.n_initial_points == 48:
            self.df_rxnid = self.df_rxnid.reindex([48,49,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                        26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47])
        elif self.n_initial_points == 24:
            self.df_rxnid = self.df_rxnid.reindex([24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])    
        else:
            raise ValueError('Cannot write legible.csv file with given number of initial points.')        

        self.df_rxnid.to_csv('output_exploration\\ISE_n_legible.csv', index=False)
        self.df_rxnid

        # determine necessary amounts for the run
        for i in self.cols:
            sum_col = sum(df_copy[i])
            print(i, ":", round(sum_col,2), "mL")

        return self.df_rxnid