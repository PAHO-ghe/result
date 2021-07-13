# COVID-19 impact on HIA

## Background 
COVID19 impact is measured in excess death in PAHO region using K-nearest neighbor method by machine learning. Impact of COVID19 is categorized by 13 UHC indicators.

## Method 

In the study, focus of impact measure is seleted as below among the UHC indicators with previosuly reported COVID19 relation

- Acute lymphoid leukemia treatment
- Asthma treatment
- Art coverage (Other immunosuppresive conditions such as HIV)
- Breast cancer treatment
- COPD treatment
- Colon/Rectum cancer treatment
- Cervical cancer treatment
- Diabetes treatment
- IHD treatmentt
- LRI treatment
- Stroke treatment
- TB treatment
- Utherine cancer treatment
- Full description of index : (https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30750-9/fulltext)

Each of UHC indicators can be interprete in the following example: 
A score of 80 in a given country-year reflects a relatively low LRI case-fatality and a score of 10 in a given country-year reflects a relatively high LRI case-fatality when compared to LRI case-fatality ratios observed in all other country-years.
Coverage of the LRI indicator is calculated using a mortality-incidence ratio (which is a proxy for case-fatality). We calculate mortality-incidence ratios for all country-years and then scale those values to range from 0 to 100, where 100 represents the best (or lowest) mortality-incidence ratio observed across all country-years and 0 is the worst (or highest) mortality-incidence ratio observed across all country-years. Most of the indicators included in the index utilize either a mortality-incidence or mortality-prevalence ratio for calculating coverage scores.

The ART indicator uses crude coverage of ART among people living with HIV. So a score of 80 in a given country-year means that we estimate that 80% of people living with HIV are on ART in that country-year. There are four indicators in the index that use crude coverage for the coverage component (both of the vaccine coverage indicators, ART coverage, and met demand for family planning).

The LRI coverage score can be compared across countries, but not across indicators within a given country, since the LRI coverage scores are scaled across countries. So we can compare Mexico’s LRI coverage score to Peru’s LRI coverage score, but we cannot compare LRI and ART treatment within Mexico.
Further details to explain the indicators can be found from Part 2, Section 2 of the reference paper (https://www.thelancet.com/cms/10.1016/S0140-6736(20)30750-9/attachment/7b01f380-587a-49d3-aae2-b6ea93898367/mmc1.pdf)

### ICD 10 codes on UHC indicators
Coming up. 

### Set environement 
```python
import os
import warnings
warnings.filterwarnings('ignore')
import sys
import csv 
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
%matplotlib inline  
from matplotlib import style
style.use("ggplot")

import seaborn as sns; sns.set()
import datetime
from datetime import datetime

from sklearn import metrics
from sklearn.cluster import KMeans

import pylab as pl

# Load data
path=(https://github.com/PAHO-ghe/result/edit/gh-pages/index.md)
filepathis=path
os.chdir(path)

df_index=pd.read_csv('./PAHO_final.csv')

df_index.isnull().sum()
df_index.columns

# EDA 
df_index['indicator_marker'] = np.where(df_index['indicator'] == 'ExcessDeathRatio', 'target',
                               np.where(df_index['indicator'].str.contains('0'), 'GHE', 'UHC'))
df_index = df_index[(df_index['indicator_marker'] != 'GHE')] # drop all GHE indicators   
df_index = df_index[(df_index['year'] == 2020) | (df_index['year'] == 2019)] # trim dataset to only include 2019 and 2020 for analysis 
df_UHCExc=df_index[(df_index['indicator'] == 'ExcessDeathRatio') | (df_index['indicator'] == 'UHC effective coverage index')]
df_index = df_index[(df_index['indicator'] != 'UHC effective coverage index')&(df_index['indicator'] != 'Covid19DeathPercent')]
df_index=df_index[['year', 'country', 'indicator', 'value']] #Select key columns 


df=df_index.pivot(index='country', columns='indicator',values='value') # reshape df_index to wide format to bind it with df_weight
df=df.dropna()

# Call GHE catalog 
df_catalog=pd.read_csv((https://github.com/PAHO-ghe/result/edit/gh-pages/index.md))
df_catalog=df_catalog.rename(columns={"ghecause": "indicator","causename":"causename" })
df_catalog['indicator'] = df_catalog['indicator'].apply(str)

# Build a merged dataset with PAHO_final & GHE catalog to patch up the GHE cause names to the GHE cause code in the original dataset
df_merged = pd.merge(df_index, df_catalog,how="left", on=["indicator"])
df_merged['finalindicator']=df_merged['causename'].combine_first(df_merged['indicator'])
#Select key columns 
df_merged=df_merged[['year', 'country', 'finalindicator', 'value']]

#Preserve old dfs
df_index_old=df_index
df_index=df_merged.rename(columns={"finalindicator": "indicator"})
# reshape df_index to wide format
df=df_index.pivot(index='country', columns='indicator',values='value')
df=df.dropna()
```


### Target: Set output as Excess death with COVID19 
```python
targetoutput='ExcessDeathRatio'
df_original=df #Use this to display final outcome

# Distance based model. Need feature scaling. 
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MaxAbsScaler
from sklearn.preprocessing import RobustScaler

#scaler = MinMaxScaler() 
#scaler = StandardScaler()
#scaler = MaxAbsScaler()
scaler = RobustScaler()

df_norm = pd.DataFrame(scaler.fit_transform(df_original.values), columns=df_original.columns, index=df_original.index)
df_original['ExcessDeathRatio'] = pd.qcut(df_original[targetoutput], 10,labels = False)
df=df_original

```
### Data : Distribution and correlation 
```python
df_UHCExc=df_UHCExc.pivot(index='country', columns='indicator',values='value') #Correlation between UHC coverage index and Excess death ratio
column_1 = df_UHCExc["ExcessDeathRatio"]
column_2 = df_UHCExc["UHC effective coverage index"]
correlation = column_1.corr(column_2)

plt.figure(figsize=(16, 10))

corr=df.corr()
mask = np.triu(np.ones_like(corr, dtype=np.bool))
cut_off = 0.7  # only show cells with abs(correlation) at least this value
extreme_1 = 0.80  # show with a star
extreme_2 = 0.90  # show with a second star
extreme_3 = 0.95  # show with a third star
mask |= np.abs(corr) < cut_off
corr = corr[~mask]  # fill in NaN in the non-desired cells

remove_empty_rows_and_cols = True
if remove_empty_rows_and_cols:
    wanted_cols = np.flatnonzero(np.count_nonzero(~mask, axis=1))
    wanted_rows = np.flatnonzero(np.count_nonzero(~mask, axis=0))
    corr = corr.iloc[wanted_cols, wanted_rows]

annot = [[f"{val:.2f}"
          + ('' if abs(val) < extreme_1 else '*')  # add one star if abs(val) >= extreme_1
          + ('' if abs(val) < extreme_2 else '*')  # add an extra star if abs(val) >= extreme_2
          + ('' if abs(val) < extreme_3 else '*')  # add yet an extra star if abs(val) >= extreme_3
          for val in row] for row in corr.to_numpy()]

heatmap = sns.heatmap(corr, vmin=0, vmax=1, annot=annot, fmt='', cmap='coolwarm',linewidths=1.0)
heatmap.set_title('Triangle Correlation Heatmap', fontdict={'fontsize': 18}, pad=16)
plt.show()


```
![graph1](https://user-images.githubusercontent.com/81782228/125119701-53cbf300-e0a6-11eb-9e87-fc5af95ad40d.png)

There are several relations between indicators and GHE causes with high correlation to each others however, none of them shows strong correlation (>0.6) with COVID19 induced death percentage per population.

```python
corr_all=df.corr()
df_highcorr=pd.DataFrame(corr_all[targetoutput])
corr_out=df_highcorr.round(2).head(len(df_highcorr))
corr_out.sort_values(by=targetoutput)
```
![Screen Shot 2021-07-09 at 11 12 43 AM](https://user-images.githubusercontent.com/81782228/125119962-a60d1400-e0a6-11eb-96d2-1771da282aad.png)


### Model 
KNN model computes relation to all indicators in the model and return unsupervised classification result of those indicator's as in format of clusters. In the analysis, the ML model computes distance of each possible combinations of UHC indicators and identify hidden patterns of class. Further details of model validation and result follows as below:
#### Cluster size 
Get SSE to check optimal cluster size
![sse](https://user-images.githubusercontent.com/81782228/125120152-eec4cd00-e0a6-11eb-90bc-908c2e4c51ca.png)

Based on SSE check among clusters sized from 2 to 15, the optimal size of cluster is either 2 or 3 with its delta ibiggest among all tested intervals. The pattern is also observed in the above graph with kinked point at where k=2 and k=3.


#### Number of countries in each clusters 

![graph4](https://user-images.githubusercontent.com/81782228/125120282-174cc700-e0a7-11eb-85db-a4ae125a607f.png)
![graph4-1](https://user-images.githubusercontent.com/81782228/125120315-26cc1000-e0a7-11eb-8825-145392bc64f2.png)
![graph4-2](https://user-images.githubusercontent.com/81782228/125120327-292e6a00-e0a7-11eb-8075-c7cd983961e7.png)

As in process to further validation of the cluster sanity, balance between clusters have checked. Balance between the cluster means how many number of samples allocated in each clusters. When the ML model built a model with 2 clusters, its balabce was xx% and xx% where the majority of case (xx%) were allocated in single cluster. Next ML model was built on k=3 scenario, where all countries are distributed into three clusters with balance of xx%, xx%, and 20%. However, final decision on what # of clusters to use for the analysis should be decided by how the target index is distributed per cluster. Therefore, the analysis will be done for both k=3 and k=2 scenario.

In case with 4 cluster is check for further accuracy, however, as shown in the graph, all clusters similarly distributed between xx% to xx%.


#### Define KNN model and test scenarios of K=2 and K=3

##### UHC coverage index mean value per cluster

| Cluster label      | Mean COVID19 impact as in excess death ratio (K=2) |Mean COVID19 impact as in excess death ratio (K=3) |
| ----------- | ----------- | ----------- |
| 0      | 3.7      | 4.2|
| 1   | 6.6        |7.1 |
| 2   |-        |3.2|


##### UHC indicator mean per cluster 

| Cluster label      | Cluster 0 | Cluster 1 |Cluster 2 |
| ----------- | ----------- | ----------- |----------- |	
|ART	|63.2	|78.4|	65.8|
|Acute lymphoid leukaemia	|8.8	|49.3	|20.7|
|Asthma	|48.0	|78.6|	68.0|
|Breast cancer	|55.5|	87.4	|73.0|
|CKD	|14.1	|50.9|	25.8|
|COPD	|52.7|	78.4|	66.2|
|Cervical cancer	|43.9	|76.6|	60.0|
|Colon/rectum cancer	|40.2|	82.3|	59.7|
|Diabetes	|35.8	|75.0	|40.4|
|IHD	|46.4|	81.6|	68.5|
|LRI	|82.7	|96.1|	92.8|
|Stroke	|41.0	|83.3|	60.4|
|TB	|65.8	|85.9|	81.9|
|Uterine cancer	|48.7	|86.3|	69.4|




| Cluster label      | Cluster 0 | Cluster 1 |
| ----------- | ----------- | ----------- |
|ART							|64.8	|77.0|
|Acute lymphoid leukaemia	|15.3	|48.0|
|Asthma						|60.4	|76.9|
|Breast cancer				|66.4	|85.2|
|CKD						|	20.3|	50.4|
|COPD						|60.8	|77.4|
|Cervical cancer			|	53.2|	76.1|
|Colon/rectum cancer		|	52.2|	79.6|
|Diabetes					|37.9	|72.9}
|IHD						|	58.9|	82.5|
|LRI						|	88.6|	96.2|
|Stroke						|52.8	|80.8|
|TB							|75.2	|86.4|
|Uterine cancer				|61.1	|84.9|


![download](https://user-images.githubusercontent.com/81782228/125319211-2aa69f00-e300-11eb-86c9-665136f34b94.png)

|             | High score cluster| Low score cluster|
| ----------- | -----------       | ----------- |
|9							| 	Chile, Panama, Peru, USA|   |
|8							| 	|   Argentina, Bolivia (Plurinational state of), Ecuador|
|7							| 	Canada, Colombia|   Brazil|
|6							| 	|   Paraguary, Belize, Honduras|
|5							| 	Cuba, |   Nicaragua, Domunical Republic|
|4							| 	|   Guyana, Haiti, Venezuela (Bolivarian Republic of)|
|3							| 	|  Guatemala, Jamaica, Mexico |
|2							| 	Costa Rica, Uruguay|   Saint Vincent and the Grenadines|
|1							| 	|   Bahamas, Trinidad and Tobago, Antigua and Barbuda|
|0							| 	|   Grenada, Barbados, Saint Lucia, Suriname|

![download (1)](https://user-images.githubusercontent.com/81782228/125319461-65a8d280-e300-11eb-816e-032b8e6ecef9.png)

|               | High score cluster| Mid score cluster                  |Low score cluster|
| -----------   | -----------       | -----------                        |----------- |
|9							| Chile, Panama, Peru, USA	|                            |              |
|8							| 	                |  Bolivia (Plurinational State Of)  |  Argentina, Ecuador|
|7							| 	Canada, Colombia|                                     |  Brazil|
|6							| 	                | Belize, Honduras                     |  Paraguay|
|5							| Cuba	            | Dominican Republic                   |  Nicaragua|
|4							| 	                |  Guyana, Haiti                      |  Venezuela (Bolivarian Republic of)|
|3							| 	                |   Guatemala                         |  Jamaica, Mexico|
|2							| Costa Rica	      | Saint Vincent and the Grenadines     |  Uruguay|
|1							| 	                |                                      |  Antigua and Barbuda, Bahanas, Trinidad and Tobago|
|0							| 	                |   Suriname                            |  Barbados, Grenada, Saint Lucia|


Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here.
## Result
### Distribution of COVID19 impact score between clusters
![download (3)](https://user-images.githubusercontent.com/81782228/125319974-e071ed80-e300-11eb-929c-e9dbbf2e15ec.png)
Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here.

### Pattern of UHC scores impact on high/low COVID19 impact 
![download (4)](https://user-images.githubusercontent.com/81782228/125320678-90475b00-e301-11eb-84a1-2c442e0bca43.png)
Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here.


![download (5)](https://user-images.githubusercontent.com/81782228/125322782-be2d9f00-e303-11eb-8dc8-ed8b35398e0f.png)
Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here. Graph description goes here.

## Conclusion 
Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. Conclusion message goes here. 












###### housekeeping matters
[editor on GitHub](https://github.com/PAHO-ghe/result/edit/gh-pages/index.md) layout and styles from the Jekyll theme  [repository settings](https://github.com/PAHO-ghe/result/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

