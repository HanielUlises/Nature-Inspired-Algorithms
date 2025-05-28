import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

data = pd.read_csv('ACO_Results.csv')

# Fitting ANOVA model
model = ols('Fit ~ C(Num_Ants) + C(Alpha) + C(Beta) + C(Num_Ants):C(Alpha) + C(Num_Ants):C(Beta) + C(Alpha):C(Beta)', data=data).fit()

# Analysis and results
anova_results = sm.stats.anova_lm(model, typ=2)
print(anova_results)
