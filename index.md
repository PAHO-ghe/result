# COVID-19 impact on HIA

## Background 

## Method 
### Set environement 
```markdown
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
```


### Target: Set output as Excess death with COVID19 
```markdown

```
### Data : Distribution and correlation 
```markdown

```

### Model 
#### Cluster size 
```markdown

```

#### Number of countries in each clusters 
```markdown

```

#### Define KNN model and test scenarios of K=2 and K=3
```markdown

```

##### UHC coverage index mean value per cluster
```markdown

```

##### UHC indicators per cluster 
```markdown

```


## Result
### Distribution of COVID19 impact score between clusters
```markdown

```

### Distribution patter on all GHE and UHC scores impact on high/low COVID19 impact 
```markdown

```

### Pattern of UHC scores in 2019 amonghighly impacted countries by COVID 19 at 2020 
```markdown

```

### Indicator pattern in countries with less impacted by covid19 
```markdown

```

## Conclusion 






You can use the [editor on GitHub](https://github.com/PAHO-ghe/result/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/PAHO-ghe/result/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
