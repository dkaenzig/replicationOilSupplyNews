# Oil supply news
Matlab code for SVAR using the oil supply surprise series constructed using high-frequency oil price variation around OPEC announcements as an external instrument to identify a structural oil supply shock, compute the dynamic causal effects as well as the historical decomposition.

**Reference**: Känzig (2020), "The macroeconomic effects of oil supply news: Evidence from OPEC announcements", https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3185839 (paper+appendix).

Tested in: Matlab 2019b on Windows 10 (64-bit)

# Contents

**[mainAnalysisOilSupplyNews.m](mainAnalysisOilSupplyNews.m)**: Main shell to reproduce all results. For external instruments
estimator, set estType = 'proxy'; for heteroskedasticity-based set estType    = 'hetero';

**[subfiles](subfiles):** Subroutines called in mainAnalysisOilSupplyNews.m

**[auxfiles](auxfiles):** Matlab routines to estimate VARX using OLS, as well as other subroutines.

**[data](data):** Data for analysis

**[figures](figures):** Results