# COVID-19 impact on HIA

## Background 
Impact of COVID19 is measured in the Americas. Impact is expressed in the size of excess death by COVDI19 in 2020 in PAHO region. For the analysis K-nearest neighbor based unsupervised classification method was used by machine learning. The impact of COVID19 is evaluated by fourteen UHC indicator categories.


## Method 

In the study, focus of impact measure is categorized as below fourteen areas, which is among the UHC indicators developed by IHME.
- Acute lymphoid leukemia treatment
- Asthma treatment
- Art coverage (Other immunosuppressive conditions such as HIV)
- Breast cancer treatment
- COPD treatment
- Colon/Rectum cancer treatment
- Cervical cancer treatment
- Diabetes treatment
- IHD treatment
- LRI treatment
- Stroke treatment
- TB treatment
- Uterine cancer treatment


Each of UHC indicators can be interpreted in the following example: 
A score of 80 in a given country-year reflects a relatively low LRI case-fatality and a score of 10 in a given country-year reflects a relatively high LRI case-fatality when compared to LRI case-fatality ratios observed in all other country-years.
Coverage of the LRI indicator is calculated using a mortality-incidence ratio (which is a proxy for case-fatality). We calculate mortality-incidence ratios for all country-years and then scale those values to range from 0 to 100, where 100 represents the best (or lowest) mortality-incidence ratio observed across all country-years and 0 is the worst (or highest) mortality-incidence ratio observed across all country-years. Most of the indicators included in the index utilize either a mortality-incidence or mortality-prevalence ratio for calculating coverage scores.

The ART indicator uses crude coverage of ART among people living with HIV. So a score of 80 in a given country-year means that we estimate that 80% of people living with HIV are on ART in that country-year. There are four indicators in the index that use crude coverage for the coverage component (both of the vaccine coverage indicators, ART coverage, and met demand for family planning).

The LRI coverage score can be compared across countries, but not across indicators within a given country, since the LRI coverage scores are scaled across countries. So we can compare Mexico’s LRI coverage score to Peru’s LRI coverage score, but we cannot compare LRI and ART treatment within Mexico.
Further details to explain the indicators can be found from Part 2, Section 2 of the reference paper [paper link](https://www.thelancet.com/cms/10.1016/S0140-6736(20)30750-9/attachment/7b01f380-587a-49d3-aae2-b6ea93898367/mmc1.pdf)

### ICD 10 codes on UHC indicators
 The ICD codes is pulled from IHME’s Global Health Data Exchange [IHME data exchange](http://ghdx.healthdata.org/record/ihme-data/gbd-2019-cause-icd-code-mappings) and refer to ICD codes used for the different causes of death and disability for the GBD 2019 study. Below table lists the ICD codes associated with all of the different causes of death and disability included in the UHC effective coverage index and therefore the ICD codes associated with each indicator.
 For each indicator, there are two components, 1) the coverage and 2) the weight. Taking the ART coverage indicator as an example, the coverage is crude coverage of ART among people living with HIV and the weight is the total potential health gain of the ART coverage intervention in each country-year. One of the inputs to the weight is the disability adjusted life-years due to HIV/AIDS. GBD 2019 estimates are the inputs to the UHC effective coverage index. Full description of index : [paper link](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30750-9/fulltext)


|	UHC ECI indicator					|	ICD10	|
| ------------------------------------- | ----------- |
|	Acute lymphoid leukemia treatment	|	C91.0-C91.02, C91.2-C91.32, C91.6-C91.62	|
|	Asthma treatment					|	J45-J46.0, Z82.5	|
|	ART coverage						|	B20-B23.8, B24-B24.0, B97.81, C46-C46.52, C46.7-C46.9, F02.4, O98.7-O98.73, Z11.4, Z20.6, Z21, Z83.0	|
|	Breast cancer treatment				|	C50-C50.629, C50.8-C50.929, Z12.3-Z12.39, Z80.3, Z85.3, Z86.000	|
|	COPD treatment						|	J41-J42.4, J43-J44.9	|
|	Colon and rectum cancer treatment	|	C18-C19.0, C20, C21-C21.8, Z12.1-Z12.13, Z85.03-Z85.048, Z86.010	|
|	Cervical cancer teratmnnt			|	C53-C53.9, Z12.4, Z85.41	|
|	Diabetes treatment					|	E08-E08.11, E08.3-E08.9, E10-E10.11, E10.3-E11.1, E11.3-E12.1, E12.3-E13.11, E13.3-E14.1, E14.3-E14.9, R73-R73.9, Z13.1, Z83.3	|
|	IHD treatment						|	I20-I21.6, I21.9-I25.9, Z82.4-Z82.49	|
|	LRI treatment						|	A48.1, A70, B96.0-B96.1, B97.21, B97.4-B97.6, J09-J18.2, J18.8-J18.9, J19.6-J22.9, J85.1, J91.0, P23-P23.9, U04-U04.9, Z25.1	|
|	Stroke treatment					|	G45-G46.8, I60-I62, I62.9-I64, I64.1, I65-I69.998, Z82.3	|
|	TB treatment						|	A10-A14, A15    -A18.89, A19-A19.9, B90-B90.9, K67.3, K93.0, M49.0, N74.0-N74.1, P37.0, U84.3, Z03.0, Z11.1, Z20.1, Z23.2	|
|	Uterine cancer treatment			|	C54-C54.3, C54.8-C54.9, Z85.42, Z86.001	|


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


### Set target: Excess death by COVID-19 in year 2020
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
### Data correlation
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

There are several indicators (CKD, Cervical cancer, Uterine cancer, Colon/rectum cancer, Breast cancer, with high correlation (>0.85*) to each other however, none of them shows significant strong correlation with COVID19 induced excess death as in the below table. 

```python
corr_all=df.corr()
df_highcorr=pd.DataFrame(corr_all[targetoutput])
corr_out=df_highcorr.round(2).head(len(df_highcorr))
corr_out.sort_values(by=targetoutput)
```

| Indicator  | Correlation coefficient to excess death rate |
| ----------- | ----------- | 
| TB     | -0.01     | 
| ART   | 0.02       |
| LRI   |0.10       |
| COPD   |0.13       |
|    Breast cancer |     0.18   |
|   CKD  |       0.24 |
|   Cervical cancer  |       0.30 |
| Acute lymphoid leukaemia    |       0.33 |
|  Colon/rectum cancer   |       0.33 |
|   Uterine cancer  |       0.33 |
|  IHD   |        0.34|
|  Asthma   |       0.36 |
|  Stroke   |        0.36|
|   Diabetes  |        0.37|

### Data distribution
UHC score distribution across countries: Below graph shows fourteen UHC indicators value across the 32 PAHO region countries. It shows that most countries have LRI, TB, and breast cancer as the first three high scored UHC indicators. Among these three countries, the rank of countries changes but stays relatively persistent as in the order that is shown in the graph. High score in UHC indicator means high case-fatality in that specific country at year 2020. However, this insight itself cannot clarify what area of UHC indicators contributed to explain high/low impact by COVID19 excess death. Therefore, k-nearest neighbor model classifies the dataset for this analysis. 
![download (6)](https://user-images.githubusercontent.com/81782228/125655094-ae42000d-d6c0-443c-b165-6d0c119d4f39.png)


### Model 
K-nearest neighbor model computes relation of all indicators in the dataset and returns its classes based on similarity. The result of the analysis is in format of clusters which is driven thru unsupervised classification effort by machine, therefore, free of human biases or perception to draw the result. In this analysis, the K-nearest neighbor model computes the distance of each data point and its possible combinations of UHC indicators and identifies hidden patterns of its distribution. Further details of model validation and result follows as below:

#### Model validation 1. Cluster size 
K-nearest neighbor models with different cluster sizes were tested. For this validation, a total of fourteen models were tested with two to fifteen centroids set as a seed of the cluster. Based on the sum of square error (SSE) value among clusters sized from 2 to 15, the optimal size of cluster can be decided either 2 or 3 as shown in the below graph. The graph shows a sudden kinked line when the centroid is set as 2 or 3. For a more detailed view, its difference of SSE was computed and the biggest difference among all tested intervals were also shown between when the centroid is set as 2 or 3.

![sse](https://user-images.githubusercontent.com/81782228/125120152-eec4cd00-e0a6-11eb-90bc-908c2e4c51ca.png)

| Centroids | Difference between SSE (%) |
| ----- | ----- | 
| 2  | -41    | 
| 3  | -18    |
| 4  |-19     |
| 5  |  -8    |
| 6  | -10    |
| 7  |  -12   |
| 8  |   -13  |
| 9  |   -12  |
| 10 |    -9  |
| 11 |    -4  |
| 12 |   -15  |
| 13 |    -5  |
| 14 |  -15   |
| 15 |   -    |


#### Model validation 2. Number of countries in each cluster
As in process to further validation of the cluster sanity, balance between clusters have been checked. Balance between the clusters means how many samples are allocated in each cluster. Ideal case for the cluster is each class containing a similar size of samples, such as 50% if k=2 case and 33% on k=3 case. In our dataset, when the model built a model with 2 clusters, its balance was 72% (23 countries out of a total of 32 countries) and 28% (9 countries/32 countries, where the majority of cases were allocated in a single cluster. Next model was built on k=3 scenario, where all countries are distributed into three clusters with a balance of 28%, 25%, and 47%. However, the final decision on cluster size for the analysis should be decided by how the target index is distributed per cluster on top of the number of countries per cluster. Therefore, the analysis will be done for both k=3 and k=2 scenarios. A model with four clusters is checked for further accuracy, however, as shown in the graph, one of the clusters contains less than 5% of data, therefore, it will not be considered further. 

![graph4](https://user-images.githubusercontent.com/81782228/125120282-174cc700-e0a7-11eb-85db-a4ae125a607f.png)
![graph4-1](https://user-images.githubusercontent.com/81782228/125120315-26cc1000-e0a7-11eb-8825-145392bc64f2.png)
![graph4-2](https://user-images.githubusercontent.com/81782228/125120327-292e6a00-e0a7-11eb-8075-c7cd983961e7.png)



#### Model validation 3. Distribution of indicators between two models (K=2 vs K=3)

##### Mean COVID19 impact mean value per cluster
Since the model computes similarity of data as in distance-based method, data should be standardized. Therefore, excess death rate data is standardized into 10 buckets of countries, sorted, and ranked for its distance-based analysis.
When the model is built by two centroids, it returns two clusters with the first cluster including countries with an excess death rate score 3.7 and 6.6 in the second cluster. Which shows, two clusters that are driven based on UHC similarity can be classed into a group with high (>=6.6) mean score of excess death compared to a group with low (<=3.7) mean score of excess death with clear division between the clusters. However, the model with three clusters, returns three clusters with each cluster having mean of the class as 4.2, 7.1 and 3.2. However, the separation between the first and third cluster with its mean as 4.2 and 3.2 can be not as clear as the previous model with two clusters and might have an overlap between clusters. Therefore, in addition to it, each indicator's mean value should be checked per cluster. 

| Cluster label | Mean COVID19 impact as in excess death ratio (K=2) |Mean COVID19 impact as in excess death ratio (K=3) |
| ----------- | ----------- | ----------- |
| 0      | 3.7   (low impact)   | 4.2 (mid impact)|
| 1   | 6.6     (high impact)   |7.1 (high impact)|
| 2   |-        |3.2 (low impact)|


#### Model validation 4. UHC indicators mean per cluster 

|     | Cluster 0 (Low impact group) | Cluster 1 (High impact) |
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


|     | Cluster 0 (Mid impact)| Cluster 1 (High impact) |Cluster 2 (Low impact) |
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



#### Model validation 5. Distribution of countries per clusters
When model is built with two clusters, it returns two classes with high/low impacted by excess death score. First cluster has its mean score on 6.6 and includes 9 countries. Out of that, four countries (44%) with highly impacted by COVID19 as it can be seen on excess death score above eight. Two countries are between 8-6, 1 country in 5, two countries in range of 2-4. A class with lower average excess death score includes countries with a rather similar proportion (~13%) across 10 groups, however the bottom score bucket with range 0-1 includes four countries out of 23 countries (17%). Further detail of countries in each class can be found in the below table.


![download](https://user-images.githubusercontent.com/81782228/125497587-bd77dcdd-50ec-4783-8faa-6a113502ef60.png)

|   Decile rank          | High score cluster| Low score cluster|
| ----------- | -----------       | ----------- |
|9							| 	Chile, Panama, Peru, USA|   |
|8							| 	|   Argentina, Bolivia (Plurinational state of), Ecuador|
|7							| 	Canada, Colombia|   Brazil|
|6							| 	|   Paraguary, Belize, Honduras|
|5							| 	Cuba, |   Nicaragua, Domunican Republic|
|4							| 	|   Guyana, Haiti, Venezuela (Bolivarian Republic of)|
|3							| 	|  Guatemala, Jamaica, Mexico |
|2							| 	Costa Rica, Uruguay|   Saint Vincent and the Grenadines|
|1							| 	|   Bahamas, Trinidad and Tobago, Antigua and Barbuda|
|0							| 	|   Grenada, Barbados, Saint Lucia, Suriname|


When a model is built with three clusters, it returns three classes with high/mid/low impacted by excess death score. It can be compared with the above model and it is shown that most of the countries from the two cluster model's highly impacted group are repeated on the three cluster model's highly impacted group. Mid/low impacted group is shown as a split of the low impacted group from the two cluster models. Further detail of countries in each class can be found in the below table.

![download (2)](https://user-images.githubusercontent.com/81782228/125497901-780b24d4-bb06-407a-8137-75329ec81649.png)

|   Decile rank            | High score cluster| Mid score cluster                  |Low score cluster|
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




## Result
### Pattern of UHC scores impact on high/low COVID19 impact 
The model classified nine countries within PAHO regions as a group of countries that is highly impacted by COVID19, which is measured by excess death. The order of high case-fatality(average) across this group's countries are as below : 
- LRI
- TB
- Breast cancer
- Uterine cancer
- IHD
- Stroke
- Colon/rectum cancer
- COPD
- ART
- Asthma
- Cervical cancer
- Diabetes
- CKD
- Acute lymphoid leukemia 

In terms of the country wise comparison, Canada and USA were top ranked in its highly impacted countries. These two countries also have high case-fatality (expressed in dark navy color cube) by different types of cancers and stroke, unlike other countries where most of countries in this group shows relatively low case-fatality (expressed in light blue color cube) by other UHC indicator area than the LRI and TB. However, it requires attention to interpret this outcome, that health surveillance and data quality can be different between countries and diagnosis sensitivity can be different between countries. 

![highimpactedcountries_heatmap](https://user-images.githubusercontent.com/81782228/125704178-593dddd6-12c8-4d10-a12a-bf2953be5745.png)

The model classified twenty-three countries within PAHO regions as a group of countries that is less impacted by COVID19, which is measured by excess death. The order of high case-fatality (average) across this group's countries are as below : 
- LRI
- TB
- Breast cancer
- ART 
- Uterine cancer 
- COPD 
- Asthma 
- IHD
- Cervical cancer
- Stroke
- Colon/rectum cancer
- Diabetes
- CKD
- Acute lymphoid leukemia 

As the graph shows, the order of top three high case-fatality UHC indicator areas and bottom three high case-fatality UHC indicator areas are the same, however with different intensity. This shows that those seven UHC indicator areas have a similar pattern across all PAHO region countries but with less of intensity among the countries with less of COVDI19 excess death. Also, this group of countries does not show persistently high case-fatality among certain UHC indicator areas. 

![lessimpactedcountries_heatmap](https://user-images.githubusercontent.com/81782228/125704202-74580836-2ca8-4417-9868-ef026878c563.png)

Below box graph shows distribution of each UHC indicator across countries within each cluster. This box graph is ordered in ascending order of mean value on UHC indicator average within the highly impacted countries. It allows comparison between two classes on its case-fatality. Both clusters have the same order of case-fatality among the first three and last three indicators. However, you can compare eight UHC indicator areas that are acting differently between the clusters. Among less impacted countries, uterine cancer and asthma case-fatality was higher than highly impacted countries and ART, COPD, and cervical cancer showed higher case-fatality. 
![download (5)](https://user-images.githubusercontent.com/81782228/125322782-be2d9f00-e303-11eb-8dc8-ed8b35398e0f.png)


## Conclusion 
COVID-19 impacted countries in the PAHO region differently. Two groups of countries have been identified that are highly impacted or less impacted, which is measured in excess death by COVDID-19. Nice countries with high COVID-19 impact include two countries in North America. And those countries highly impacted by COVID-19 showed higher case-fatality in all UHC indicator areas. Particularly, countries with higher impact by COVID-19 showed higher burden of uterine cancer, IHD, and stroke in the year before the pandemic compared to less impacted countries, excluding LRI, TB, and breast cancer burden, which was similar pattern as the less impacted countries.

This finding shows COVID19 impact was higher when countries already carried the burden of case-fatality of LRI, TB, and breast cancer across all PAHO region countries. It is supported by the result that shows, countries with lower case-fatality of these three areas of UHC indicator had lower excess death score due to COVID19. Therefore, lowering LRI, TB, and breast cancer case-fatality can lead to less impact from COVID-19 like pandemic.

Lastly, the noteworthy point of this machine learning study shows that there are countries with strong similarity of health system performance that explains higher impact by COVID-19.









###### housekeeping matters
[editor on GitHub](https://github.com/PAHO-ghe/result/edit/gh-pages/index.md) 
[repository settings](https://github.com/PAHO-ghe/result/settings/pages)  
`_config.yml` 

