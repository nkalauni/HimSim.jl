

import pandas as pd
import numpy as np
import os


data = pd.read_csv('other_params.csv')


data.to_csv('chepe_data_modify.csv', index=False)