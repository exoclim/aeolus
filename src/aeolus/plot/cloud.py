"""Utilities to plot cloud-related output of UM."""
import warnings

import iris

import matplotlib as mpl

import numpy as np

from .cm_custom import cloudtypes_cmap


__all__ = ("CloudPlotter", "cloudtypes_cmap")


class CloudPlotter:
    """Factory to create a composite plot of low, medium, high cloud area fractions."""

    _stash_items = {"lo": "m01s09i203", "me": "m01s09i204", "hi": "m01s09i205"}
    _scales = {"lo": 100, "me": 200, "hi": 400}
    cloud_labels = [
        "No\nCloud",
        "L",
        "M",
        "L + M",
        "H",
        "H + L",
        "H + M",
        "All",
    ]

    def __init__(self, cubelist):
        """Initialise CloudPlotter from `iris.cube.CubeList` containing cloud fractions."""
        self.factor = 10.0  # scaling factor
        self.cubes = {}
        for key, stash in self._stash_items.items():
            try:
                self.cubes[key] = cubelist.extract_cube(iris.AttributeConstraint(STASH=stash))
            except iris.exceptions.ConstraintMismatchError:
                warnings.warn(f"Warning!\n{key} ({stash}) is not found in\n\n{cubelist}")

    def scale_data(self, threshold=0.1):
        """
        Scale cloud levels and add them together to make a composite cloud field.

        Values below `threshold` are set to zero.
        """
        for i, (key, cube) in enumerate(self.cubes.items()):
            cp_cube = cube.copy(
                data=np.where(
                    cube.data >= threshold, cube.data * self.factor + self._scales[key], 0.0
                )
            )
            if i == 0:
                self.aggr_cld = cp_cube
            else:
                self.aggr_cld += cp_cube
        self.aggr_cld.rename("aggregated_cloud_fraction")
        self.aggr_cld.attributes["source"] = " + ".join(
            [
                (f"({self._stash_items[k]}*{self.factor}" f"+{self._scales[k]})")
                for k in self.cubes.keys()
            ]
        )

    def _make_norm(self, cmap, nsteps):
        """
        Create `BoundaryNorm` for cloud levels.

        Uses `nsteps` gradations of colour and `cmap` colormap.
        Also returns `midpoints` used for colorbar labelling.
        """
        base = np.linspace(0, self.factor, nsteps) + 1
        assert nsteps > 0, "number of intervals should be positive"
        lvls = [base]  # no clouds
        lvls.append(base + self._scales["lo"])  # low clouds (101-110)
        lvls.append(base + self._scales["me"])  # medium clouds (201-210)
        lvls.append(base * 2 + self._scales["lo"] + self._scales["me"])  # low + medium (302-320)
        lvls.append(base + self._scales["hi"])  # high clouds (401-410)
        lvls.append(base * 2 + self._scales["lo"] + self._scales["hi"])  # high + low (502-520)
        lvls.append(base * 2 + self._scales["hi"] + self._scales["me"])  # high + medium (602-620)
        lvls.append(
            base * 3 + self._scales["lo"] + self._scales["me"] + self._scales["hi"]
        )  # high + medium + low (703-730)
        lvls.append([1000])
        midpoints = []
        if nsteps == 1:
            for lvl1, lvl2 in zip(lvls[:-1], lvls[1:]):
                midpoints.append(0.5 * (lvl1[-1] + lvl2[0]))
        else:
            for lvl in lvls[:-1]:
                if len(lvl) % 2 == 0:
                    midpoints.append(lvl[len(lvl) // 2])
                else:
                    midpoints.append(0.5 * (lvl[len(lvl) // 2] + lvl[(len(lvl) // 2) + 1]))
        midpoints = np.asarray(midpoints)
        norm = mpl.colors.BoundaryNorm(np.concatenate(lvls), cmap.N)
        return norm, midpoints

    def get_all(self, nsteps=10, cmap=None):
        """Get cloud types cube and kwargs for external plotting."""
        if cmap is None:
            cmap = cloudtypes_cmap  # colormap from a separate file
        norm, cb_ticks = self._make_norm(cmap=cmap, nsteps=nsteps)
        plt_kw = {"norm": norm, "cmap": cmap}
        cb_labels = self.cloud_labels
        return self.aggr_cld, plt_kw, cb_ticks, cb_labels

    def pcolormesh(self, ax, xy=(), nsteps=10, cmap=None, cb_kwargs={}, **pc_kwargs):
        """
        Display composite cloud plot using `matplotlib.pyplot.pcolormesh`.

        Parameters
        ----------
        ax: matplotlib axes object
            Axes where to create the plot in
        xy: tuple, optional
            Sequence of X and Y coordinate arrays (passed to pcolormesh)
        nsteps: integer, optional
            Number of colour steps for each of cloud category
        cmap: matplotlib.colors.ListedColormap, optional
            Colormap instance. If not given, custom cloud colormap is used.
        cb_kwags: dict, optional
            Dictionary of kwargs passed to colorbar
        **pc_kwargs: optional
            additional pcolormesh arguments
        Returns
        -------
        h: `matplotlib.collections.QuadMesh`
            pcolormesh instance
        cb: `matplotlib.colorbar.Colorbar`
            colorbar
        """
        fig = ax.figure
        if cmap is None:
            cmap = cloudtypes_cmap  # colormap from a separate file
        norm, midpoints = self._make_norm(cmap=cmap, nsteps=nsteps)
        h = ax.pcolormesh(*xy, self.aggr_cld.data, cmap=cmap, norm=norm, **pc_kwargs)
        if "cax" in cb_kwargs:
            ax_kw = {}
        else:
            ax_kw = {"ax": ax}
        cb = fig.colorbar(h, **ax_kw, **cb_kwargs)
        cb.set_ticks(midpoints)
        cb.set_ticklabels(self.cloud_labels)
        return h, cb
