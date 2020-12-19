# Oil supply news
Matlab code for SVAR using the oil supply surprise series constructed using high-frequency oil price variation around OPEC announcements as an external instrument to identify a structural oil supply shock, compute the dynamic causal effects as well as the historical decomposition.

**Reference**: Känzig (2020), "The macroeconomic effects of oil supply news: Evidence from OPEC announcements", https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3185839 (paper+appendix).

Tested in: Matlab 2019b on Windows 10 (64-bit)

# Contents

**[mainAnalysisOilSupplyNews.m](mainAnalysisOilSupplyNews.m)**: Main shell to reproduce all results. For external instruments
estimator, set `estType = 'proxy'`; for heteroskedasticity-based set `estType = 'hetero'`.

**[subfiles](subfiles):** Subroutines called in `mainAnalysisOilSupplyNews.m`
- [transformAndPlotData.m](subfiles/transformAndPlotData.m): function to transform (and plot) raw data
- [loadProxy.m](subfiles/loadProxy.m): function to read in external instrument
- [runProxyVAR.m](functions/runProxyVAR.m): function to estimate proxy VAR; bands are computed using bootstrapping techniques
- [runRigobonVAR.m](functions/runRigobonVAR.m): function to estimate heteroskedasticity-based VAR; bands are computed using bootstrapping techniques

**[auxfiles](auxfiles):** Matlab routines to estimate VARX using OLS, as well as other subroutines.

**[data](data):** Data for analysis
- [OilDataM.mat](data/OilDataM.mat): raw data used in VAR
- [OilSurprisesMLog.mat](data/OilSurprisesMLog.mat): oil supply surprise series constructed using high-frequency approach
- [OilSurprisesMLogControl.mat](data/OilSurprisesMLogControl.mat): control series 

**[figures](figures):** Stores results from analysis
