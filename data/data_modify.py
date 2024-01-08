import pandas as pd
import numpy as np
import os


data = pd.read_csv('chepe_data.csv')

print(data.head())



data.to_csv('chepe_data.csv', index=False)