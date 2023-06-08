import pandas as pd
import plotly.express as px

df = pd.read_excel("/mnt/c/Users/0/Documents/EAFIT/Tesis/results.xlsx",
                   sheet_name='Results')

px.bar(df, 'Underlying model', 'Outliers', color='Model', barmode='group')
px.bar(df, 'Underlying model', 'No Outliers', color='Model', barmode='group')
df_c = pd.concat(
    [
        df[['Underlying model', 'Outliers']],
        df[['Underlying model', 'No Outliers']].rename(columns={'No Outliers': 'Outliers'})
    ]
).rename(columns={'Outliers': 'Percentage of Preference', 'Underlying model': 'Underlying_Model'})
models = pd.concat((df['Model'] , df['Model'] + '_o'))
df_c['Model'] = models

px.bar(df_c, 'Underlying_Model', 'Percentage of Preference', color='Model', barmode='group')
px.bar(df_c.query('Underlying_Model.str.startswith("A")'), 'Underlying_Model', 'Percentage of Preference', color='Model', barmode='group')
px.bar(df_c.query('Underlying_Model.str.startswith("M") and ~Underlying_Model.str.endswith("M")'), 'Underlying_Model', 'Percentage of Preference', color='Model', barmode='group')
px.bar(df_c.query('Underlying_Model.str.startswith("M") and Underlying_Model.str.endswith("M")'), 'Underlying_Model', 'Percentage of Preference', color='Model', barmode='group')

df['Text'] = df['Outliers'].apply(lambda x: f"{x:,.2f}") + ' (' + df['No Outliers'].apply(lambda x: f"{x:,.2f}") + ')'
display(
    df[['Underlying model', 'Text', 'Model']]
    .pivot(index='Underlying model', columns='Model', values='Text')
    )
display(
    df[['Underlying model', 'Outliers', 'No Outliers', 'Model']]
    .pivot(index='Underlying model', columns='Model', values=['Outliers', 'No Outliers'])
    )
1