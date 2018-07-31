# Consistent-thresholds-RTCC
We developed a new methodological (trait-based) approach to identify distinct drivers of ecological community assembly. The method can first separate the combined effects of biotic and abiotic factors on changes of community composition along an environmental gradient, and, second, find consistent thresholds in the value of an environmental variable along the gradient indicating substantial structural changes in the leading drivers. Although we have developed this approach for microbial communities, it is completely general and can be applied to macro-organismal communities with the same easiness. To carry out this approach, a dataset of (phenotypic and/or genotypic) community traits and metadata from the variables are necessary. This method can help understand how the environment affects community assembly, which enables quantitative predictions for the ecological impact of environmental changes. This research can be found under the title "Consistent environmental thresholds in community assembly".

## The scripts 
Script 1: gen-Metacom.py;
Script 2: boxplot.py;
Script 3: hindex-repetitions.py;
Script 4: Chow-test.R;
RTCC C code: permut.c;
OSX shared library: permlib.so

## Additional data tables 
- Additional data table S1 (traits_plus_range.csv):Example data of traits and environmental ranges of the pool of species in the metacommunity.
- Additional data table S2 (table_presence_absence.csv): Example presence-absence data frame.
- Additional data table S3 (metadata.csv): Example metadata: a wide salinity gradient from 0.1 to 40% of dissolved salts.
- Additional data table S4 (metacomm_sum.csv): Dataset of relative abundances (between 0 and 1) for each species, normalized to the total abundance observed in the metacommunity.
- Additional data table S5 (clustering_index_vs_env_gradient.csv): Dataset for the curve of the clustering index values for a given trait along the environmental gradient of a variable resulting from the STEP 3 in the protocol (see supplementary information of the related article).
