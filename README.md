# TwinGeneticFactor

- The source code named "calc_genetic_nogender_mix.py" calculates GEIs for each methylation site.

## What we calculated:

We calculate std (standard deviations) of two distributions; twin and non-twin pairs' methylation levels.
 
GEIs is [std. of non-twin] / [std. of twin]

We adopt Levene's test to assess the equality of variance of twin and non-twin pairs.


## Requirement

Run on python >3

We used to calculate statistics by 

'''
from scipy.stats import levene
import statistics
'''

## Data for this code (Data Availability)

### Data used
Monozygonal twin pair: 	twin_mono_intersection.pickle
Methylation file:	methyl_T_intersection.pickle

### Result file
df_analysis_res.pickle

### How to get the data

Data are not publicly available.

The methylation levels of each site in each subject were used under license in the current study and are, thus, not publicly available. However, the data are available from the "Twin Research Center, Osaka University, Japan," with permission from the center.

The details of the result are at http://www.res.kutc.kansai-u.ac.jp/~takenaka/twin/
(Not open public before the article is accepted. Ask Yoichi Takenaka: takenaka@kansai-u.ac.jp for ID and Password)
