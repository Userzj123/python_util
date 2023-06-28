from turb.lesgo_utils import lesgo_data
import numpy as np
import turb.lesgo_utils
import matplotlib.pyplot as plt

result_dir = '/home/zyou6474/tasks/channel_flow'
dims = (128, 128, 64)
domain = (2*np.pi, np.pi, 1)

ldata = lesgo_data(domain, dims, result_dir, ntheta=3)


ldata.read_data(1)
ldata.data['theta'].shape