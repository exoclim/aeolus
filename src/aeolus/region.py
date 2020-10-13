"""Subsetting variables over geographical regions."""
import iris
from iris.analysis.cartography import wrap_lons

from .exceptions import BoundaryError
from .model import um
from .plot.text import fmt_lonlat


__all__ = ("Region",)


class BoundsRect:
    """Bounding longitudes and latitudes of a given lon-lat rectangle."""

    _side_names = ("west", "east", "south", "north")

    def __init__(self, west_bound, east_bound, south_bound, north_bound, model=um):
        """Initialise BoundsRect."""
        self.west = west_bound
        self.east = east_bound
        self.south = south_bound
        self.north = north_bound
        self.model = model
        self.west_coord = model.x
        self.east_coord = model.x
        self.south_coord = model.y
        self.north_coord = model.y
        if self.south > self.north:
            raise BoundaryError(
                f"South boundary value ({self.south}) should be less than north ({self.north})"
            )
        self.sides = [(name, getattr(self, f"{name}_coord")) for name in self._side_names]

    def __repr__(self):  # noqa
        return (
            f"BoundsRect(west={self.west}, east={self.east}, south={self.south}, "
            f"north={self.north})"
        )


class Region:
    """
    Rectangular geographical region.

    Attributes
    ----------
    name : str
        The region's name
    description : str
        A description of the region
    constraint: iris.Constraint
        A constraint object associated with the region
    """

    def __init__(
        self, west_bound, east_bound, south_bound, north_bound, name="", description="", model=um
    ):
        """
        Instantiate a `Region` object.

        Parameters
        ----------
        west_bound, east_bound, south_bound, north_bound : scalar, optional
            The western, eastern, southern, and northern boundaries, respectively, of the
            region.
        name: str, optional
            The region's name.
        description : str, optional
            A description of the region.
        model: aeolus.model.Model, optional
            Model class with a relevant coordinate names.
        """
        self.name = name
        self.description = description
        self.model = model

        self.bounds = BoundsRect(west_bound, east_bound, south_bound, north_bound, model)
        self._sides = self.bounds.sides

        self.lon_size = abs(self.bounds.east - self.bounds.west)
        self.lat_size = self.bounds.north - self.bounds.south

    def __repr__(self):  # noqa
        txt = (
            f"Geographical region '{self.name}' (west={self.bounds.west}, "
            f"east={self.bounds.east}, south={self.bounds.south}, north={self.bounds.north})"
        )
        if self.description:
            txt += "\n\n"
            txt += self.description
        return txt

    def __getitem__(self, index):  # noqa
        return {
            "value": getattr(self.bounds, self._sides[index][0]),
            "name": self._sides[index][0],
            "coord": self._sides[index][1],
        }

    def _perpendicular_side_limits(self, side):
        """Get minimum and maximum values of the region boundary perpendicular to the given one."""
        if side in ["west", "east"]:
            coord_name = self.model.y
            _min, _max = self.bounds.south, self.bounds.north
        elif side in ["south", "north"]:
            coord_name = self.model.x
            _min, _max = self.bounds.west, self.bounds.east
        else:
            raise BoundaryError(f"Boundary name '{side}' is not valid")
        return coord_name, (_min, _max)

    def to_str(self, sep="_"):  # noqa
        return sep.join([fmt_lonlat(i["value"], i["coord"]) for i in self])

    @classmethod
    def from_cube(
        cls, cube, name=None, margin=None, margin_units="points", shift_lons=False, model=um
    ):
        """
        Create a Region from limits of longitude and latitude of the cube.

        Parameters
        ----------
        cube: iris.cube.Cube
            Source cube.
        name: str, optional
            Name for the region. If not given, created automatically from `cube`'s name.
        margin: scalar, optional
            Use `margin` number of points or degrees to create a region smaller than the cube.
        margin_units: str, optional
            Units of margin. Can be "points" or "degrees".
        shift_lons: bool, optional
            Shift longitudes to -180...180.
        model: aeolus.model.Model, optional
            Model class with a relevant coordinate names.

        Returns
        -------
        aeolus.region.Region
        """
        if name is None:
            name = f"extent_of_{cube.name()}"
        lons = cube.coord(model.x).points
        if shift_lons:
            lons = sorted(wrap_lons(lons, -180, 360))
        lats = cube.coord(model.y).points
        idx0, idx1 = 0, -1
        if margin is not None:
            if margin_units == "points":
                idx0 += margin
                idx1 -= margin
                lon0 = lons[idx0]
                lon1 = lons[idx1]
                lat0 = lats[idx0]
                lat1 = lats[idx1]
            else:
                lon0 = lons[idx0] + margin
                lon1 = lons[idx1] - margin
                lat0 = lats[idx0] + margin
                lat1 = lats[idx1] - margin
        else:
            lon0 = lons[idx0]
            lon1 = lons[idx1]
            lat0 = lats[idx0]
            lat1 = lats[idx1]

        return cls(lon0, lon1, lat0, lat1, name=name)

    @property
    def constraint(self):
        """Constraint to select data within the region."""
        cnstr = iris.Constraint(
            latitude=lambda x: self.bounds.south <= x.point <= self.bounds.north
        )
        if self.bounds.west < self.bounds.east:
            # Western boundary is to the west
            cnstr &= iris.Constraint(
                longitude=lambda x: self.bounds.west <= x.point <= self.bounds.east
            )
        else:
            # Region wrapping around dateline (180deg)
            cnstr &= iris.Constraint(
                longitude=lambda x: (self.bounds.west <= x.point) or (x.point <= self.bounds.east)
            )
        return cnstr

    def add_to_ax(self, ax, **kwargs):
        """Add a Rectangle patch to matplotlib axes `ax` with given keyword arguments `kwargs`."""
        from matplotlib.patches import Rectangle  # noqa

        xy = (self.bounds.west, self.bounds.south)
        width = self.lon_size
        height = self.lat_size
        ax.add_patch(Rectangle(xy, width, height, **kwargs))
