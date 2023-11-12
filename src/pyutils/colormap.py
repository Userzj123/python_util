from matplotlib import cm
def plot_examples(colormaps):
    """
    Helper function to plot data with associated colormap.
    """
    np.random.seed(19680801)
    data = np.random.randn(30, 30)
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            layout='constrained', squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()
    
    
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

start_ind = 80
length = 80

twilight_shifted = mpl.colormaps['twilight_shifted'].resampled(256)
newcolors = twilight_shifted(np.linspace(00.5, 1, 256))

white = np.array([256/256, 256/256, 256/256, 1])
end = newcolors[start_ind, :]


color2white = np.zeros((length, 4))
for i in range(3):
    color2white[:, i] = np.linspace(white[i], end[i], length)
color2white[:, 3] = 1

newcolors = np.concatenate((color2white, newcolors[start_ind:, :]), axis=0)


positive_twilight = ListedColormap(newcolors)

plot_examples([twilight_shifted, positive_twilight])