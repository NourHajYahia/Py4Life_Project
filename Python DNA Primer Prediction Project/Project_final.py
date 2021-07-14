# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 21:45:46 2019

@author: nour
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import seaborn as sns

# filling the program with the molar absorption Coefficient of each base
nucleotides_MolarAbsorptionCoeff = {'G': 11700, 'C': 7500, 'T': 9200, 'A': 7500, 'GC': 17785, 'AT': 17183}

# Visualizing the Polynomial Regression results
def viz_polymonial(X,y,y_label):
    plt.plot(X,y, label = y_label, linewidth=3)
    plt.legend(frameon=False,fontsize =10)
    plt.xticks(rotation=60, fontsize=15, color = 'blue')
    plt.yticks(rotation='horizontal',fontsize=15, color = 'blue')
    plt.xlabel('WaveLength [nm]', fontsize = 15)
    plt.ylabel('OD', fontsize = 15)
    sns.despine()
    plt.tight_layout() 
    return

# Polynomial fitting function. Returns the polynomial constants (a0, a1, a2...)
def polymonial_fit(dig,X,y):
    poly_reg = PolynomialFeatures(degree=dig)
    X_poly = poly_reg.fit_transform(X)
    pol_reg = LinearRegression()
    pol_reg.fit(X_poly, y)
    base_coef = pol_reg.coef_
    base_coef[0] = pol_reg.intercept_
    return base_coef

# Comparing the similarity of both reference and sample nucleotides absorption spectra’s.
# Compares one sample array to each one of the spectra arrays of references nucleotide 
# Returns the comparison data, which contains R2 value of the similarity (0-1 values), the prediction and nucleotide concentration.
def Base_ref_to_sample_fit(X_Sample,y_Sample,y_Ref):
    
    # Data dictionary to save the calcuations
    data = {}
    data['Sample'] = np.array([y_Sample.name])
    
    # Ratio between two nucleotide at 260
    Wave_260_location = np.where(X_Sample == 260)[0][0]
    score = {'G': 0, 'C': 0, 'T': 0, 'A': 0, 'A strand': 0, 'C strand': 0, 'T strand': 0, 'G strand': 0, 'GC': 0, 'AT': 0}
    for column_ref in y_Ref.columns:
        y_ref_base = y_Ref[column_ref]
        Ratio = y_Sample.max()/y_ref_base.max()

        # Reference bases polynomial   fitting and plotting
        # dG and dA nucleotide specs splitting for 2 tested groups of polynomial fitting else no need for splitting the data
        if column_ref in ('G','A'):
            X = [X_Sample.iloc[0:Wave_260_location,:],X_Sample.iloc[Wave_260_location:,:]]
            y = [y_ref_base.iloc[0:Wave_260_location],y_ref_base.iloc[Wave_260_location:]]
            y_new = [0,0]
            ref_coef = [0,0]
            ref_polynom = [0,0]
            for j in range(2):
                ref_coef[j] = polymonial_fit(4,X[j].values,y[j].values)
                ref_coef[j] = Ratio * ref_coef[j]
                ref_polynom[j] = np.poly1d(ref_coef[j][::-1])
                y_new[j] = ref_polynom[j](X[j].values)  
            y_new_ref = np.concatenate((y_new[0],y_new[1]))
   
        else:  
            base_coef = polymonial_fit(4,X_Sample.values,y_ref_base.values)
            base_coef = Ratio * base_coef
            base_polynom = np.poly1d(base_coef[::-1])
            y_new_ref = base_polynom(X_Sample.values)
        
        # R^2 (coefficient of determination) regression score function.
        score[column_ref] = [r2_score(y_Sample.values,y_new_ref),y_new_ref]
        
        # If the R2 value is lower than 0 means the prediction is too far from the true, so adjusted to 0.
        if score[column_ref][0] < 0:
            score[column_ref][0]= 0
            
        # Filling the R2 score values to data dictionary
        data[column_ref] = np.array([score[column_ref][0]])
    
    # Finding the nearest prediction to the truth which is the maximum value of r2 score we get from references comparison with the sample   
    Max = max(score.items(), key=lambda k: k[1])
      
    #Filling the nearest prediction name to the data dictionary 
    data['Nuc Prediction'] = np.array([Max[0]])
    
    # Sample nucleotide concentration calcuation
    A_260 = y_Sample.iloc[Wave_260_location]
    L = 1
    concentration = A_260/(nucleotides_MolarAbsorptionCoeff[Max[0][0]]*L)
    
    # Sample nucleotide concentration calculation  saving to the dictionary
    data['Nuc Concentration [M]'] = np.array([concentration])
    return data

# Reading references csv file with pandas.
Ref = pd.read_csv('Referances.csv')
# limiting the pectra data between nucleotide absorption range (240 - 300 nm)
start_Ref = np.where(Ref['Wavelength (nm)'] == 240)[0][0]
end_Ref = np.where(Ref['Wavelength (nm)'] == 300)[0][0]
X_Reference_base = Ref.iloc[start_Ref:end_Ref +1, 0:1]
y_Reference_base = Ref.iloc[start_Ref:end_Ref +1,1:]

# Reading Samples csv file with pandas.
Sample = pd.read_csv('Samples.csv')
# limiting the pectra data between nucleotide absorption range (240 - 300 nm)
start = np.where(Sample['Wavelength (nm)'] == 240)[0][0]
end = np.where(Sample['Wavelength (nm)'] == 300)[0][0]
X_Sample_base = Sample.iloc[start:end +1, 0:1]
y_Sample_base = Sample.iloc[start:end +1,1:]

# for function which passes on each sample send it to the comparison function. and filling the returned data into the data dictionary.
data = Base_ref_to_sample_fit(X_Sample_base,y_Sample_base.iloc[:,0],y_Reference_base)
for column in y_Sample_base.columns[1:]:
    Temp = Base_ref_to_sample_fit(X_Sample_base,y_Sample_base[column],y_Reference_base)
    for keys in Temp:
        if keys in data:
            data[keys]= np.append(data[keys],Temp[keys])

# Converting the data dictionary to pandas DataFrame
Nuc_data = pd.DataFrame.from_dict(data)
Nuc_data_T = Nuc_data.transpose()
Nuc_data_T = Nuc_data_T.reset_index()
# line plot Graph of all the samples combined in one.
sns.set_style("dark", {'axes.axisbelow': True, 'axes.edgecolor': '.2',"axes.facecolor": "0.95"})
fig = plt.figure(1, figsize=(8,5))
for i in range(1,len(Sample.columns[1:])):
    Nuc_prediction = Nuc_data.loc[Nuc_data['Sample'] == Sample.columns[i],'Nuc Prediction']
    viz_polymonial(X_Sample_base.values,y_Sample_base.iloc[:,i],Sample.columns[i] + ' Is: ' + Nuc_prediction[i-1])
    plt.title('Spectroscopy Graph Of All Samples',fontsize=20,color = 'green')
plt.show()


# Seaborn bar plot of each reference nucleatides R2 coeffecient prediction values relative to each sample nucleotides.
nucs_end = np.where(Nuc_data.columns == 'Nuc Prediction')[0][0]
figure_num = int((nucs_end-1)/4)
fig = plt.figure(figure_num, figsize=(12,12))
fig.subplots_adjust(hspace=0.5, wspace=0.4)
sns.set_context("paper", rc={"font.size":20,"axes.titlesize":15,"axes.labelsize":15})
j = 0
for i in range(1,nucs_end):
    if ((i-1)/3)%1 == 0 and (i-1) != 0:
        j = i-1
    ax1 = plt.subplot(figure_num+1, 4, i, title=Nuc_data.iloc[:,i].name + ' Predictions')
    sns.barplot(x=Nuc_data.iloc[:,i].values ,y=Nuc_data.iloc[:,0].values, data=Nuc_data, palette='Spectral')
    plt.suptitle('R2 Referance Nucleotides Prediction Values Of Each Sample', fontsize=20, color = 'green')
    plt.xticks(color = 'blue')
    plt.yticks(color = 'blue')
    plt.xlabel('R2 Score', fontsize = 12)
    plt.ylabel('Samples', fontsize = 12)
plt.show()

#nucleotide pairing possibility prediction, Calculating the volume fraction to get equivalent concentration of nucleotides for best primer creation and annealing and plot prediction for the annealing product.
#AT pairing possibility.
if 'A strand' in Nuc_data['Nuc Prediction'].values and 'T strand' in Nuc_data['Nuc Prediction'].values:
    A_strand = Nuc_data.loc[Nuc_data['Nuc Prediction'] == 'A strand']
    T_strand = Nuc_data.loc[Nuc_data['Nuc Prediction'] == 'T strand']
    volume_ratio = A_strand['Nuc Concentration [M]'].iloc[0]/T_strand['Nuc Concentration [M]'].iloc[0]
    A_in_mL_volume = 1 / (1 + volume_ratio )
    print('\033[0;32mAT Pairing \033[0;30m')
    print('The mixture solution volume ratio between', A_strand['Sample'].iloc[0],'and',T_strand['Sample'].iloc[0],'to reach equivalent concentration is:','\033[0;33m%.2f' %volume_ratio,'\033[0;30m')
    print('For polyApolyT primer incubation of 1 mL mixture of','\033[0;33m%.2f' %A_in_mL_volume,'\033[0;30mmL of', A_strand['Sample'].iloc[0],'and','\033[0;33m%.2f' %(1 - A_in_mL_volume),'\033[0;30mmL of',T_strand['Sample'].iloc[0],'in 90 degrees temperature for 10 min followed by slow cooling to 30 degrees temperature and finally moving the sample to the fridge.')
    AT_concentration = A_in_mL_volume * A_strand['Nuc Concentration [M]'].iloc[0]
    OD_260_AT = AT_concentration * 1 * nucleotides_MolarAbsorptionCoeff['AT']
    Wave_260_location = np.where(X_Sample_base == 260)[0][0]
    ratio = OD_260_AT/y_Reference_base.iloc[Wave_260_location]['AT']
    base_coef_AT = polymonial_fit(4,X_Reference_base,y_Reference_base['AT'])
    base_coef_AT = ratio * base_coef_AT
    base_polynom = np.poly1d(base_coef_AT[::-1])
    y_AT_predicted = base_polynom(X_Sample_base.values)
    fig = plt.figure(1, figsize=(8,5))
    viz_polymonial(X_Reference_base,y_AT_predicted,'AT')
    plt.title('Pairs Prediction',fontsize=20,color = 'green')

#GC pairing possibility.     
if 'G strand' in Nuc_data['Nuc Prediction'].values and 'C strand' in Nuc_data['Nuc Prediction'].values:
    G_strand = Nuc_data.loc[Nuc_data['Nuc Prediction'] == 'G strand']
    C_strand = Nuc_data.loc[Nuc_data['Nuc Prediction'] == 'C strand']
    volume_ratio = G_strand['Nuc Concentration [M]'].iloc[0]/C_strand['Nuc Concentration [M]'].iloc[0]
    G_in_mL_volume = 1 / (1 + volume_ratio )
    print('\n\033[0;32mGC Pairing \033[0;30m')
    print('The volume ratio between', G_strand['Sample'].iloc[0],'and',C_strand['Sample'].iloc[0],'is:','\033[0;33m%.2f' %volume_ratio,'\033[0;30m')
    print('For polyGpolyC primer incubation of 1 mL mixture of','\033[0;33m%.2f' %G_in_mL_volume,'\033[0;30mmL of', G_strand['Sample'].iloc[0],'and','\033[0;33m%.2f' %(1 - G_in_mL_volume),'\033[0;30m1mL of',C_strand['Sample'].iloc[0],'in 90 degrees temperature for 10 min followed by slow cooling to 30 degrees temperature and finally moving the sample to the fridge.')
    GC_concentration = G_in_mL_volume * G_strand['Nuc Concentration [M]'].iloc[0]
    OD_260_GC = GC_concentration * 1 * nucleotides_MolarAbsorptionCoeff['GC']
    ratio = OD_260_GC/y_Reference_base.iloc[Wave_260_location]['GC']
    base_coef_GC = polymonial_fit(4,X_Reference_base,y_Reference_base['GC'])
    base_coef_GC = ratio * base_coef_GC
    base_polynom = np.poly1d(base_coef_GC[::-1])
    y_GC_predicted = base_polynom(X_Sample_base.values)
    fig = plt.figure(1, figsize=(10,5))
    viz_polymonial(X_Reference_base,y_GC_predicted,'GC')
    plt.show()
     
else:
    print('There is no AT ot GC pairs in your Data.')                   
          
#writing the dataframe to csv file names Result.csv       
Nuc_data.to_csv('Results.csv', index=False)         
          