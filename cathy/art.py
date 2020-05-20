
import numpy
import scipy.ndimage
from matplotlib import pyplot

import catheter_utils.projections


class ProjectionsAnimator:

    def __init__(self, projections_meta, projections_raw):
        self.meta = projections_meta
        self.projections = [catheter_utils.projections.reconstruct(raw) for raw in projections_raw]
        self.signal_max = max(numpy.amax(p) for p in self.projections)
        self.signal_min = min(numpy.amin(p) for p in self.projections)
        max_fov = projections_meta["fov"].max()
        self.domain_max = max_fov/2.0
        self.domain_min = -self.domain_max

    def _pick(self, i, axis):
        fov = self.meta["fov"][i]
        fs = self.projections[i][axis]
        xs = numpy.linspace(-fov / 2, fov / 2, len(fs))
        return xs, fs

    def plot(self, i, axis, ax=None, plot_args=None, plot_kwargs=None):
        if not plot_args:
            plot_args = []
        if not plot_kwargs:
            plot_kwargs = {}

        xs, fs = self._pick(i, axis)
        if ax:
            return ax.plot(xs, fs, *plot_args, **plot_kwargs)
        else:
            return pyplot.plot(xs, fs, *plot_args, **plot_kwargs)

    def animate(self, i, axis, plot=None):
        xs, fs = self._pick(i, axis)
        if plot:
            plot.set_data(xs, fs)

    @property
    def xlim(self):
        return self.domain_min, self.domain_max

    @property
    def ylim(self):
        return self.signal_min, self.signal_max

