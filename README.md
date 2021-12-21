# Oil supply news
Matlab code for SVAR using the oil supply surprise series constructed using high-frequency oil price variation around OPEC announcements as an external instrument to identify a structural oil supply news shock, compute the dynamic causal effects as well as the historical decomposition

**Reference**: Känzig (2021), "The macroeconomic effects of oil supply news: Evidence from OPEC announcements", https://www.aeaweb.org/articles?id=10.1257/aer.20190964&&from=f (article), https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3185839 (working paper+appendix)

Tested in: Matlab 2019b on Windows 10 (64-bit)

# Contents

**[mainAnalysisOilSupplyNews.m](mainAnalysisOilSupplyNews.m)**: Main shell to reproduce all results. For external instruments
estimator, set `estType = 'proxy'`; for heteroskedasticity-based set `estType = 'hetero'`

**[\subfiles](subfiles):** Subscripts called in `mainAnalysisOilSupplyNews.m`
- [transformAndPlotData.m](subfiles/transformAndPlotData.m): script to transform (and plot) raw data
- [loadProxy.m](subfiles/loadProxy.m): script to read in external instrument
- [runProxyVAR.m](subfiles/runProxyVAR.m): script to estimate proxy VAR; bands are computed using bootstrapping techniques
- [runRigobonVAR.m](subfiles/runRigobonVAR.m): script to estimate heteroskedasticity-based VAR; bands are computed using bootstrapping techniques

**[\auxfiles](auxfiles):** Matlab functions to estimate VARX using OLS, as well as other subroutines

**[\data](data):** Data for analysis
- [OilDataM.mat](data/OilDataM.mat): raw data used in VAR
- [OilSurprisesMLog.mat](data/OilSurprisesMLog.mat): oil supply surprise series constructed using high-frequency approach
- [OilSurprisesMLogControl.mat](data/OilSurprisesMLogControl.mat): control series 

**[\figures](figures):** Stores results from analysis


**[paper](Känzig&#32;2020&#32;-&#32;The&#32;macroeconomic&#32;effects&#32;of&#32;oil&#32;supply&#32;news.pdf):** Pdf containing paper and online appendix
