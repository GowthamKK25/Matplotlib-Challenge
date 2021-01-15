#!/usr/bin/env python
# coding: utf-8

# In[171]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np 


# Study data files
mouse_metadata_path = "../data/Mouse_metadata.csv"
study_results_path = "../data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path) 

# Combine the data into a single dataset and cleaning data
mouse = pd.merge(mouse_metadata,study_results,on=["Mouse ID"])
mouse = mouse.drop_duplicates(subset = ["Mouse ID","Timepoint"], keep = False) 
mouse


# In[3]:


# Count mouse 
mouse["Mouse ID"].nunique()


# In[4]:


#Summary Stats
Drugs = mouse.groupby("Drug Regimen")  
x = Drugs["Tumor Volume (mm3)"].mean() 
y = Drugs["Tumor Volume (mm3)"].median() 
z = Drugs["Tumor Volume (mm3)"].var() 
a = Drugs["Tumor Volume (mm3)"].std()  
b = Drugs["Tumor Volume (mm3)"].sem()  

sumstats = pd.DataFrame({"Mean Tumor Volume": x,"Median Tumor Volume": y,"Tumor Volume Variance":z,"Tumor Volume Std. Dev.":a,"Tumor Volume Std. Err.":b})


# In[5]:


sumstats


# In[6]:


# Using Pandas to plot drug reg.
drugplotpoints = Drugs["Mouse ID"].count() 
drugplotpoints.plot(kind="bar",figsize=(10,5)) 
plt.title("Drug Regimen")
plt.xlabel("Drug Type")
plt.ylabel("Mouses per Drug") 


# In[7]:


#Using mat plot to plot drug reg.
x_axis = mouse["Drug Regimen"].unique()
x_axis 
y_axis = Drugs["Mouse ID"].count()  
y_axis 

 
plt.title("Drug Regimen")
plt.xlabel("Drug Type")
plt.ylabel("Mouses per Drug") 
plt.bar(x_axis, y_axis,alpha=1, align='center',width = .5)
plt.xticks(rotation="vertical")
plt.figure(figsize=(20, 5))


# In[8]:


#Using Pandas for pie chart 

gendermice = mouse.groupby(["Mouse ID","Sex"]) 
mouse_gender_df = pd.DataFrame(gendermice.size())
mouse_gender = pd.DataFrame(mouse_gender_df.groupby(["Sex"]).count())
mouse_gender.columns = ["Total Count"] 
mouse_gender.plot.pie(y= 'Total Count',figsize=(5,5),shadow = True, autopct="%1.1f%%")


# In[9]:


#Using Mat plot lib for pie chart 

values = [49.8, 50.2]
labels = ["Female","Male"]
plt.pie(values,autopct="%1.1f%%", shadow=True,labels=labels) 

plt.axis("equal")


# In[386]:


drugs = mouse.groupby('Mouse ID').max()['Timepoint']

drugs = pd.DataFrame(drugs)
drugs = pd.merge(drugs, mouse, on=("Mouse ID","Timepoint"),how="left")
x = drugs 
x


# In[379]:


# Filering my data toward the specfied drugs
Capomulin = x.loc[x["Drug Regimen"] == "Capomulin"] 
Ramicane = x.loc[x["Drug Regimen"] == "Ramicane"] 
Infubinol = x.loc[x["Drug Regimen"] == "Infubinol"] 
Ceftamin = x.loc[x["Drug Regimen"] == "Ceftamin"] 


# In[380]:


# Quartile analysis for each drug
#Capomulin
quartiles = Capomulin["Tumor Volume (mm3)"].quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Values below {lower_bound} could be outliers for Capomulin.")
print(f"Values above {upper_bound} could be outliers for Capomulin.") 


# In[381]:


#Ramicane
quartiles = Ramicane["Tumor Volume (mm3)"].quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Values below {lower_bound} could be outliers for Ramicane.")
print(f"Values above {upper_bound} could be outliers for Ramicane.")  


# In[382]:


#Infubinol
quartiles = Infubinol["Tumor Volume (mm3)"].quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Values below {lower_bound} could be outliers for Infubinol.")
print(f"Values above {upper_bound} could be outliers for Infubinol.")  


# In[383]:


#Ceftamin 
quartiles = Ceftamin["Tumor Volume (mm3)"].quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Values below {lower_bound} could be outliers for Ceftamin.")
print(f"Values above {upper_bound} could be outliers for Ceftamin.")  


# In[384]:


data = [Ceftamin["Tumor Volume (mm3)"], Infubinol["Tumor Volume (mm3)"], Capomulin["Tumor Volume (mm3)"], Capomulin["Tumor Volume (mm3)"]]

fig1, ax1 = plt.subplots()
ax1.set_title('Tumors')
ax1.set_ylabel('Final Tumor Volume (mm3)')
ax1.set_xlabel('Drug Regimen')

ax1.boxplot(data, labels=["Ceftamin","Infubinol","Capomulin","Capomulin"])

plt.ylim(20,80)
plt.show()


# In[385]:


#Select a mouse that was treated with Capomulin and generate a line plot of tumor volume vs. time point for that mouse.


# In[318]:


#List of Mouse Id's to select from
Capomulin


# In[338]:


#Selecting mouse ID for capomulin drug
Mouse_Id = input(f'Insert one of following mouse ids ') 


# In[342]:


#Getting data for x and y axes
Timepoint = mouse.loc[Mouse_Id,["Timepoint"]] 

Tumorvol = mouse.loc[Mouse_Id,["Tumor Volume (mm3)"]] 


# In[345]:


#Plot Chart
plt.plot(Tumorvol, Timepoint, color="green") 


plt.title(Mouse_Id + "'s tumor volume vs. time point graph")
plt.xlabel("tumor volume")

plt.ylabel("time point")


plt.show()


# In[347]:


#Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin treatment regimen.


# In[354]:


#Get data for both axes
Capoweight = Capomulin.loc[:,"Weight (g)"] 
Capovol =  Capomulin.loc[:,"Tumor Volume (mm3)"] 


# In[355]:


#Plot chart
plt.scatter(Capoweight, Capovol, color="green") 


plt.title("mouse weight versus average tumor volume")
plt.xlabel("tumor volume")

plt.ylabel("time point")


plt.show()


# In[393]:


#Plot linear Regression
from scipy.stats import linregress

x_values = Capoweight
y_values = Capovol
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(x_values,y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(6,10),fontsize=15,color="green")
plt.title("mouse weight versus average tumor volume")
plt.xlabel("tumor volume")
plt.ylabel("time point")
plt.show()


# In[ ]:




