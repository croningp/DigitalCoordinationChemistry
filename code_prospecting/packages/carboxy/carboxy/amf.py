### CODE: DANIEL KOWALSKI

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
from collections import Counter

######################################################################################
class validate_dataset():
######################################################################################
    def __init__(
        self, spec_list, exploration_type:str, ratio_csv_list:list, intensity_cutoff:float=0.05
    ) -> None:

        self.spec_list = spec_list
        self.ratio_csv_list = ratio_csv_list
        self.type = exploration_type   # either 'cce' or 'ise'
        self.intensity_cutoff = intensity_cutoff

        # determine correct input given for 'exploration_type'
        if self.type == 'CCE':
            self.type = 'cce'
        if self.type == 'ISE':
            self.type = 'ise'
        if (self.type != 'cce') and (self.type != 'ise'):
            raise ValueError("Only permissable values for exploration_type are 'cce' or 'ise'.")


    def validate_intensity(self):
        # setup blank lists
        self.max_intensities = []
        self.usable_pn=[]

        # pull max intensity data
        for i in self.spec_list:
            self.max_intensities.append(max(i.data['Intensity']))
        
        # calculate cutoff as a percentage of maximum intensity
        limit = max(self.max_intensities)*self.intensity_cutoff

        # create list indicating which samples are above limit
        for i in self.max_intensities:
            if i >= limit:
                self.usable_pn.append(1)
            else:
                self.usable_pn.append(0)

        return self.usable_pn, limit


    def validate_reagents(self):
        # import dataframes
        df_list = []
        for item in self.ratio_csv_list:
            upload_df = pd.read_csv(item, index_col=0)
            df_list.append(upload_df)

        validation_df = pd.concat(df_list)
        ### NOTE: Need to make this conditional (if statement)
        try:
            validation_df = validation_df.drop('TYPE',axis=0)
            validation_df.reset_index(drop=True, inplace=True)
        except:
            pass

        # cut cluster data out
        if self.type == 'cce':
            sub_df = validation_df[['Co3O(OH)','Co3O','Co4O4']]

        if self.type == 'ise':
            sub_df = validation_df[['Cr3O','Mn3O','Fe3O','Co3O']]
        
        # create array of values
        sub_arr = copy.deepcopy(sub_df.values)

        # iterate through array and replace any non-zero value with a 1
        m=0
        for i in sub_arr:
            n = 0
            for j in i:
                if float(j) > 0.0:
                    sub_arr[m,n] = 1
                n+=1
            m += 1

        # rename
        self.usable_ri = sub_arr

        return self.usable_ri


#########################################################################################################
class Calculate_AMF_MS():
#########################################################################################################    
    def __init__(
        self, ms_spectrum, standards_list:list=[], weighting:int=1,
        threshold:float=0.0, fingerprint_size:int=5, 
        precision_pn:int=2, precision_ri:int=1,
        validation_pn:list=[], validation_ri:list=[],
        verbose:bool=False
    ) -> None:
        
        self.ms_spectrum = ms_spectrum
        self.standards_list = standards_list
        self.weighting = weighting
        self.threshold = threshold
        self.fingerprint_size = fingerprint_size
        self.precision_pn = precision_pn
        self.precision_ri = precision_ri
        self.validation_pn = validation_pn
        self.validation_ri = validation_ri
        self.verbose = verbose
    
    
    def solve(self):
        self.peak_novelty()
        self.compare_standards()
        self.relative_intensities()
    
    
    def report(self, save:bool=False, file_name="AMF_MS_result.txt"):
        stds = []
        for i in self.standards_list:
            stds.append(i.label)
        
        # display
        print(
            "SAMPLE:", self.ms_spectrum.label,
            "\nSTANDARDS:", stds,
            "\n   PEAK NOVELTY:",
            "\n      Unique:", len(self.unique_elements), self.unique_elements,
            "\n   CONTAINING STANDARDS:",
        )
        
        # print to .txt file
        if save == True:
            file_object = open(file_name, "a")

            str1 = (
                "SAMPLE: " + self.ms_spectrum.label +
                "\nSTANDARDS: " + str(stds) +
                "\n   PEAK NOVELTY:" +
                "\n      Unique: " + str(len(self.unique_elements)) + " " + str(self.unique_elements) + 
                "\n   CONTAINING STANDARDS:"
            )
            
            file_object.write(str1)
        
        # duplicate and cut down differences list to determine how many peaks were found
        for i in range(len(self.ri['Stds'])):
            differences_copy = copy.deepcopy(self.ri['Differences'][i])
            differences_copy = [x for x in differences_copy if x != 100]
            
            # display
            print(
                "      % Difference between", self.ri['Stds'][i], 
                "in standard and sample:", round(self.ri['Average Differences'][i],2),
                "based on", len(differences_copy), "matching peaks"
            )
            
            # print to .txt file
            if save == True:
                str2 = (
                    "\n      % Difference between " + str(self.ri['Stds'][i]) +
                    " in standard and sample: " + str(round(self.ri['Average Differences'][i],2)) +
                    " based on " + str(len(differences_copy)) + " matching peaks"
                )
                
                file_object.write(str2)
                
        # display
        print()
        
        #print to .txt file
        if save == True:
            file_object.write("\n \n")
            file_object.close()
                            
    
    def sample_above_threshold(ms_spectrum, threshold):
        # deepcopy DataFrame
        df = copy.deepcopy(ms_spectrum.data)
        
        # flag large DataFrames
        if len(df) > 30000:
            error_prompt = str("DataFrame appears to be quite large. Do you wish to continue? (y/N)" + 
                "\nDATAFRAME LENGTH: " + str(len(df)) + "\n")

            permission = input(error_prompt)
            if permission == 'y':
                pass
            else:
                raise KeyboardInterrupt('Operation aborted. DataFrame too large.')
        
        # remove any values under threshold
        for index, row in df.iterrows():
            if row['Relative Intensity'] < threshold:
                df = df.drop(index, axis='index')

        return df


    def peak_novelty(self):
        # if max intensity less than certain % of highest intensity in the dataset
        # return no unique peaks
        if self.validation_pn == 0.0:
            self.unique_elements = []
            self.repeated_elements = []
        
        else:
            # remove noise
            self.spec_pn = Calculate_AMF_MS.sample_above_threshold(self.ms_spectrum, self.threshold)
            
            # unpack m/z values to array
            spec_pn_mz = self.spec_pn['m/z'].values
            
            # remove noise and unpack m/z values to array for each standard
            std_mz = []
            std_names = []

            for i in self.standards_list:
                std_names.append(i.label)
                std = Calculate_AMF_MS.sample_above_threshold(i, self.threshold)
                for item in list(std['m/z'].values):
                    std_mz.append(item)

            # round data to appropriate precision
            spec_pn_mz = [round(num,self.precision_pn) for num in spec_pn_mz]
            std_mz = [round(num,self.precision_pn) for num in std_mz]
            
            # remove any elements from sample m/z array which appear in stndard arrays
            self.unique_elements = [elem for elem in spec_pn_mz if elem not in std_mz]
            self.repeated_elements = [elem for elem in spec_pn_mz if elem in std_mz]
        
        # display report 
        if self.verbose == True: print(
            "COMPARE PEAK POSITIONS:", self.ms_spectrum.label,
            "\n   AGAINST:", std_names,
            "\n   Total:", len(self.ms_spectrum.data),
            "\n   Over Threshold:", len(spec_pn_mz),
            "\n   Unique:", len(self.unique_elements), self.unique_elements,
            "\n"
        )

        return self.unique_elements, self.repeated_elements


    def compare_standards(self):
        # remove noise
        self.spec_ri = Calculate_AMF_MS.sample_above_threshold(self.ms_spectrum, self.threshold)

        # round m/z values to appropriate precision
        self.spec_ri['m/z'] = self.spec_ri['m/z'].round(self.precision_ri)
        
        if self.verbose == True: print("COMPARE TO STANDARDS:", self.ms_spectrum.label)
        
        # extend validation list by multiplying by the weighting factor
        self.validation_ri_extended = []
        for i in self.validation_ri:
            for _ in range(self.weighting):
                self.validation_ri_extended.append(i)

        if len(self.validation_ri_extended) != len(self.standards_list):
            raise ValueError('Extended validation list does not match length of standards length.')

        # for each standard, sort by relative intensity, round, 
        # and unpack top 'x' values into the fingerprint_list
        fingerprint_list = []

        n = 0
        for obj in self.standards_list:
            std = obj.data
            std = std.sort_values('Relative Intensity', ascending=False)
            std['m/z'] = std['m/z'].round(self.precision_ri)
            fingerprint = std.iloc[0:self.fingerprint_size]

            # if cluster not used in reaction, prep list to be all 100.0's
            #blank = []
            #for _ in range(self.fingerprint_size):
            #    blank.append(100.0)
            #    fingerprint=blank
           
            # add to fingerprint_list
            fingerprint_list.append(fingerprint)
            n+=1
        
        # setup lists to record data
        n=0
        report_stds = []
        report_i_samples = []
        report_i_stds = []
        report_lookup = []

        for i in fingerprint_list:
            # setup temporary lists to permit searching
            i_sample = np.array([])
            hold_list=[]
            
            # for each standard list the values to look up
            lookup = i['m/z'].values
            report_stds.append(self.standards_list[n].label)
            if self.verbose == True: print("   FOR:", self.standards_list[n].label, "\n   Lookup:", lookup)
            n+=1

            i_standard = i['Relative Intensity'].values

            # for each lookup value search for any matching values in the sample
            # add all matching values to the hold_list
            for j in lookup:
                for index, row in self.spec_ri.iterrows():
                    if row['m/z']==j:
                        hold_list.append(row['Relative Intensity'])

                # add a single value to the standard's list (the max)
                if len(hold_list) == 1:
                    i_sample = np.hstack([i_sample, max(hold_list)])
                    hold_list=[]
                elif len(hold_list) > 1:
                    if self.verbose == True: print("      Multiple values found for:", j, hold_list)
                    i_sample = np.hstack([i_sample, max(hold_list)])
                    hold_list=[]
                elif len(hold_list) < 1:
                    i_sample = np.hstack([i_sample, 0.0])
                    hold_list=[]
        

            # unpack temporary lists to permanent
            report_i_stds.append(i_standard.round(2))
            if self.verbose == True: print("   Standard Rel Intensities:", i_standard.round(2))
            report_i_samples.append(i_sample.round(2))
            if self.verbose == True: print("   Sample Rel Intensities:", i_sample.round(2), "\n")
            report_lookup.append(lookup)

        # define dictionary of data
        self.comparison = {
            'Sample': self.ms_spectrum.label,
            'Stds': report_stds,
            'Peak Positions': report_lookup,
            'i_standards': report_i_stds,
            'i_samples': report_i_samples,
            'fingerprint_size': self.fingerprint_size
        }

        return self.comparison


    def relative_intensities(self):
        # unpack relevant data from dictionary generated by 'AMF_MS.compare_standards()'
        standard_list = self.comparison['Stds']
        i_samples = self.comparison['i_samples']
        i_standards = self.comparison['i_standards']

        if self.verbose == True: print("COMPARE RELATIVE INTENSITIES:", self.comparison['Sample'])
        
        # setup lists to record data
        all_samples = []
        all_standards = []
        all_differences = []
        all_av_diffs = []

        # iterate through all comparisons
        for n in range(len(self.standards_list)):
            name = self.comparison['Stds'][n]
            samp = list(i_samples[n])
            std = list(i_standards[n])
            if self.verbose == True: print("   FOR:", name)

            # scale to the same height, based on max peak in the sample spectrum
            max_index = samp.index(max(samp))
            max_samp = samp[max_index]
            max_std = std[max_index]
            differences = []

            # find percentage difference between each corresponding pair of peaks
            if self.validation_ri_extended[n] == 0.0:
                for _ in range(self.fingerprint_size):
                    differences.append(0.0)
                for i in range(len(samp)):
                    std[i] = std[i]/max_std
            else:
                for i in range(len(samp)):
                    # 0 % where 1:1 matching peak
                    if samp[i] == 100.0:
                        differences.append(0.0)
                        std[i] = std[i]/max_std
                    # 100 % where peak not found
                    elif samp[i] == 0.0:
                        differences.append(100.0)
                        std[i] = std[i]/max_std                
                    else:
                        samp[i] = samp[i]/max_samp
                        std[i] = std[i]/max_std

                        differences.append(((std[i]-samp[i])/std[i])*100)
            
            # ensure absolute value used, non-directional
            differences = [abs(i) for i in differences]
           
            av_diff = (sum(differences)/(self.comparison['fingerprint_size']-1))   # '-1' term added as one value 
                                                                                       # will always be 0 % different
            
            # constrain values where (i) no peaks are matched, or (ii) where sample peaks are more intense
            if av_diff > 100:
                av_diff = 100
                if self.verbose == True: print("av_diff value over 100")

            if av_diff < 0:
                av_diff = 0
                if self.verbose == True: print("av_diff value under 0")

            if np.isnan(av_diff):
                av_diff = 0

            # report data
            if self.verbose == True: print(
                "   RELATIVE INTENSITIES:", 
                "\n      Sample:", samp, 
                "\n      Standards:", std,
                "\n      Differences:", differences,
                "\n   AVERAGE % DIFFERENCE BETWEEN SCALED PEAKS:", round(av_diff,2)
            )
            
            # transfer data to permanent lists
            all_samples.append(samp)
            all_standards.append(std)
            all_differences.append(differences)
            all_av_diffs.append(av_diff)
            differences=[]

        # define dictionary of data
        self.ri = {
            'Sample': self.comparison['Sample'],
            'Stds': self.comparison['Stds'],
            'Sample Rel Values': all_samples,
            'Standard Rel Values': all_standards,
            'Differences': all_differences,
            'Average Differences': all_av_diffs
        }

        return self.ri
    


#########################################################################################################
class AMF_MS_CCE:
#########################################################################################################
    def __init__(
        self, ri_list, ue_list, validation_ri, weighting
    ):

        self.ri_list = ri_list
        self.ue_list = ue_list
        self.validation_ri = validation_ri
        self.weighting = weighting


    def populate_lists(self):      
        # permanent lists
        self.list_std1 = []
        self.list_std2 = []
        self.list_std3 = []
        self.list_std4 = []
        self.list_std5 = []
        self.list_std6 = []
        self.diffs = []

        # populate lists
        n = 0
        for i in self.ri_list:
            division = self.weighting * sum(self.validation_ri[n])
            print(n, self.validation_ri[n], division)

            self.list_std1.append(i[0])
            self.list_std2.append(i[1])
            self.list_std3.append(i[2])
            self.list_std4.append(i[3])
            self.list_std5.append(i[4])
            self.list_std6.append(i[5])
            x = (float(i[0])+float(i[1])+float(i[2])+float(i[3])+float(i[4])+float(i[5]))/division
            self.diffs.append(x)

            n+=1
      
    
    def uniqueness(self):
        # setup temporary lists
        all_ue = []
        self.uniques = []

        # combine all values from each set of unique peaks to a single list
        for i in self.ue_list:
            for j in range(len(i)):
                all_ue.append(i[j])

        # sort list
        all_ue = sorted(all_ue)

        # count instances of each value
        self.count = Counter(all_ue)
        print(self.count)

        # bind to DataFrame
        self.df_counts = pd.DataFrame.from_dict(self.count, orient='index', columns=['count']).reset_index()

        # for each sample, lookup counts of each unique peak and add 1/counts to list with peak m/z
        for i in self.ue_list:
            count_list = []
            for j in i:
                for index, row in self.df_counts.iterrows():
                    if j == row['index']:
                        count_list.append([row['index'],(1/row['count'])])

            # setup temporary list
            hold_list = []

            # collect 1/counts terms
            for n in count_list:
                hold_list.append(n[1])

            # sum and collect in list
            unique_val = sum(hold_list)
            self.uniques.append(unique_val)

        # feature scale between 0 and 1
        max_score = max(self.uniques)
        min_score = min(self.uniques)

        self.uniques_arr = np.array(self.uniques)

        self.uniques_prime = (self.uniques_arr - min_score)/(max_score - min_score)
        

    def rel_intensities(self):
        #NEEDS COMMENTS
        self.amf_ms_score = []

        max_diff = max(self.diffs)
        min_diff = min(self.diffs)

        self.diff_arr = np.array(self.diffs)

        if max_diff==min_diff:
            self.diff_prime = self.diff_arr/100
        else:
            self.diff_prime = (self.diff_arr - min_diff)/(max_diff - min_diff)


    def AMF_Score(self):
        #NEEDS COMMENTS
        for i in range(len(self.diff_arr)):   
            diff_score = self.diff_prime[i]
            unique_score = self.uniques_prime[i]

            score = (diff_score + unique_score)/2
            self.amf_ms_score.append(score)


        DICT = {
            "Co3O(OH)": self.list_std1,
            "Co3O(OH)": self.list_std2,
            "Co3O": self.list_std3,
            "Co3O Degradation": self.list_std4,
            "Co4O4": self.list_std5,
            "Co4O4": self.list_std6,
            "Average % Diff": self.diffs,
            "% Diff Scaled": self.diff_prime,
            "Uniqueness": self.uniques,
            "Uniqueness Scaled": self.uniques_prime,
            "Score": self.amf_ms_score
        }
        #pd.set_option('display.max_rows', 10)
        self.result = pd.DataFrame(DICT)
        self.result = self.result.sort_values(by='Score', ascending=False)
        
        return self.result



#########################################################################################################
class AMF_MS_ISE:
#########################################################################################################
    def __init__(
        self, ri_list, ue_list, validation_ri, weighting
    ):

        self.ri_list = ri_list
        self.ue_list = ue_list      
        self.validation_ri = validation_ri
        self.weighting = weighting


    def populate_lists(self):             
        # permanent lists
        self.list_std1 = []
        self.list_std2 = []
        self.list_std3 = []
        self.list_std4 = []
        self.list_std5 = []
        self.list_std6 = []
        self.list_std7 = []
        self.list_std8 = []
        self.diffs = []

        # populate lists
        n = 0
        for i in self.ri_list:
            division = self.weighting * sum(self.validation_ri[n])
            if division == 0:
                division = 1
            print(n, self.validation_ri[n], division)
            
            self.list_std1.append(i[0])
            self.list_std2.append(i[1])
            self.list_std3.append(i[2])
            self.list_std4.append(i[3])
            # NOTE: commented the below out for the test dataset, need to put them BACK for other datasets!!
            #self.list_std5.append(i[4])
            #self.list_std6.append(i[5])
            #self.list_std7.append(i[6])
            #self.list_std8.append(i[7])
            x = (float(i[0])+float(i[1])+float(i[2])+float(i[3]))/division
            #x = (float(i[0])+float(i[1])+float(i[2])+float(i[3])+float(i[4])+float(i[5])+float(i[6])+float(i[7]))/division
            self.diffs.append(x)

            n += 1
      
    
    def uniqueness(self):
        # setup temporary lists
        all_ue = []
        self.uniques = []

        # combine all values from each set of unique peaks to a single list
        for i in self.ue_list:
            for j in range(len(i)):
                all_ue.append(i[j])

        # sort list
        all_ue = sorted(all_ue)

        # count instances of each value
        self.count = Counter(all_ue)
        print(self.count)

        # bind to DataFrame
        self.df_counts = pd.DataFrame.from_dict(self.count, orient='index', columns=['count']).reset_index()

        # for each sample, lookup counts of each unique peak and add 1/counts to list with peak m/z
        for i in self.ue_list:
            count_list = []
            for j in i:
                for index, row in self.df_counts.iterrows():
                    if j == row['index']:
                        count_list.append([row['index'],(1/row['count'])])

            # setup temporary list
            hold_list = []

            # collect 1/counts terms
            for n in count_list:
                hold_list.append(n[1])

            # sum and collect in list
            unique_val = sum(hold_list)
            self.uniques.append(unique_val)

        # feature scale between 0 and 1
        max_score = max(self.uniques)
        min_score = min(self.uniques)

        self.uniques_arr = np.array(self.uniques)

        self.uniques_prime = (self.uniques_arr - min_score)/(max_score - min_score)
        

    def rel_intensities(self):
        #NEEDS COMMENTS
        self.amf_ms_score = []

        max_diff = max(self.diffs)
        min_diff = min(self.diffs)

        self.diff_arr = np.array(self.diffs)

        if max_diff==min_diff:
            self.diff_prime = self.diff_arr/100
        else:
            self.diff_prime = (self.diff_arr - min_diff)/(max_diff - min_diff)


    def AMF_Score(self):
        #NEEDS COMMENTS
        for i in range(len(self.diff_arr)):   
            diff_score = self.diff_prime[i]
            unique_score = self.uniques_prime[i]

            score = (diff_score + unique_score)/2
            self.amf_ms_score.append(score)


        DICT = {
            "Cr3O": self.list_std1,
            "Cr3O": self.list_std2,
            "Mn3O 1": self.list_std3,
            "Mn3O 2": self.list_std4,
            # NOTE: commented the below out for the test dataset, need to put them BACK for other datasets!!
            #"Fe3O 1": self.list_std5,
            #"Fe3O 2": self.list_std6,
            #"Co3O": self.list_std7,
            #"Co3O in MeOH": self.list_std8,
            "Average % Diff": self.diffs,
            "% Diff Scaled": self.diff_prime,
            "Uniqueness": self.uniques,
            "Uniqueness Scaled": self.uniques_prime,
            "Score": self.amf_ms_score
        }

        self.result = pd.DataFrame(DICT)
        self.result = self.result.sort_values(by='Score', ascending=False)
        
        return self.result


#########################################################################################################
class AMF_MS_Combination():
#########################################################################################################  
    def __init__(
        self, result_list, exploration_type:str
    ) -> None:
        
        self.result_list = result_list
        self.type = exploration_type   # either 'cce' or 'ise'

        if self.type == 'CCE':
            self.type = 'cce'
        if self.type == 'ISE':
            self.type = 'ise'
        if (self.type != 'cce') and (self.type != 'ise'):
            raise ValueError("Only permissable values for exploration_type are 'cce' or 'ise'.")

        self.sort_indexes()
        self.diff_from_stds()
        self.uniqueness()
        self.combine_scores()


    def sort_indexes(self):      
        for n in range(len(self.result_list)):
            self.result_list[n] = self.result_list[n].sort_index()

    
    def diff_from_stds(self):
        if self.type == 'cce':
            col = [i for i, x in enumerate(self.result_list[0].columns) if x=='Average % Diff']   # corresponds to 'Average % Diff' column
        if self.type == 'ise':
            col = [i for i, x in enumerate(self.result_list[0].columns) if x=='Average % Diff']   # corresponds to 'Average % Diff' column
        
        col = col[0]

        store_all = []

        for i in self.result_list:
            arr = i.values
            arr = arr.T
            arr = arr[col]

            store_all.append(arr)

        average = sum(store_all)/(len(self.result_list))

        min_av = min(average)
        max_av = max(average)

        self.average_diff_prime = (average - min_av)/(max_av - min_av)

        return self.average_diff_prime


    def uniqueness(self):
        if self.type == 'cce':
            col = [i for i, x in enumerate(self.result_list[0].columns) if x=='Uniqueness']   # corresponds to 'Uniqueness' column
        if self.type == 'ise':
            col = [i for i, x in enumerate(self.result_list[0].columns) if x=='Uniqueness']   # corresponds to 'Uniqueness' column

        col = col[0]

        store_all = []

        for i in self.result_list:
            arr = i.values
            arr = arr.T
            arr = arr[col]

            store_all.append(arr)

        average = sum(store_all)/(len(self.result_list))

        min_av = min(average)
        max_av = max(average)

        self.average_unique_prime = (average - min_av)/(max_av - min_av)

        return self.average_unique_prime


    def combine_scores(self):
        amf_ms_score = []

        for i in range(len(self.average_diff_prime)):
            diff_score = self.average_diff_prime[i]
            unique_score = self.average_unique_prime[i]

            score = (diff_score + unique_score)/2
            amf_ms_score.append(score)

        self.df = pd.DataFrame(amf_ms_score, columns=['Final Score'])
        
        return self.df

    
    def create_csv(self, ratio_csv_list:list, drop_metadata:bool=False):
        # import dataframes
        df_list = []
        for item in ratio_csv_list:
            upload_df = pd.read_csv(item, index_col=0)
            df_list.append(upload_df)

        self.ratio_amf_df = pd.concat(df_list)
        ### NOTE: Need to make this conditional (if statement)
        if drop_metadata == True:
            self.ratio_amf_df = self.ratio_amf_df.drop('TYPE',axis=0)
        self.ratio_amf_df.reset_index(drop=True, inplace=True)

        self.ratio_amf_df.index.name=None

        if self.type == 'cce':
            column_floats = ['Co3O(OH)','Co3O','Co4O4','rT1','rT2','V_TOT']
        if self.type == 'ise':
            column_floats = ['Cr3O','Mn3O','Fe3O','Co3O','rT1','rT2','V_TOT']
        
        for i in column_floats:
            self.ratio_amf_df[i] = self.ratio_amf_df[i].astype(float)
            
        self.ratio_amf_df['AMF'] = self.df.values#[-48:]

        if self.type == 'cce':
            self.ratio_amf_df.to_csv("output_mapping-function\\CCE_n_result_AMF.csv")
        if self.type == 'ise':
            self.ratio_amf_df.to_csv("output_mapping-function\\ISE_n_result_AMF.csv")

        return self.ratio_amf_df