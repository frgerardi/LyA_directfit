import numpy as np
import matplotlib.pyplot as plt


def callback(sampler):
    sample = sampler.collection    #[sampler.last_point_callback:]
    idx = np.linspace(1,sample['H0'].shape[0],sample['H0'].shape[0])
    plt.set_cmap('coolwarm')
    plt.scatter(sample['H0'],sample['bias_LYA'],c=idx)
    plt.savefig('sampling_callback.png')