"""Subsetting variables over geographical regions."""
from dataclasses import dataclass, field

import iris

from .exceptions import BoundaryError
from .util import fmt_lonlat


__all__ = ("Region",)


@dataclass
class BoundsRect:
    """Bounding longitudes and latitudes of a given lon-lat rectangle."""

    west: float = field(metadata={"coord": "longitude"})
    east: float = field(metadata={"coord": "longitude"})
    south: float = field(metadata={"coord": "latitude"})
    north: float = field(metadata={"coord": "latitude"})

    def __post_init__(self):  # noqa
        # if self.west > self.east:
        #     raise BoundaryError("West boundary value should be less than east")
        if self.south > self.north:
            raise BoundaryError("South boundary value should be less than north")

    def __repr__(self):  # noqa
        return (
            f"BoundsRect(west={self.west}, east={self.east}, south={self.south}, "
            f"north={self.north})"
        )


class Region(object):
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

    def __init__(self, west_bound, east_bound, south_bound, north_bound, name="", description=""):
        """
        Instantiate a `Region` object.

        Parameters
        ----------
        name: str
            The region's name.
        description : str, optional
            A description of the region.
        west_bound, east_bound, south_bound, north_bound : scalar, optional
            The western, eastern, southern, and northern boundaries, respectively, of the
            region.
        """
        self.name = name
        self.description = description

        self.bounds = BoundsRect(west_bound, east_bound, south_bound, north_bound)
        self._sides = [
            (key, f.metadata["coord"]) for key, f in self.bounds.__dataclass_fields__.items()
        ]

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
            coord_name = "latitude"
            _min, _max = self.bounds.south, self.bounds.north
        elif side in ["south", "north"]:
            coord_name = "longitude"
            _min, _max = self.bounds.west, self.bounds.east
        else:
            raise BoundaryError(f"Boundary name '{side}' is not valid")
        return coord_name, (_min, _max)

    def to_str(self, sep="_"):  # noqa
        return sep.join([fmt_lonlat(i["value"], i["coord"]) for i in self])

    @classmethod
    def from_cube(cls, cube, margin=None, margin_units="points"):
        """
        Create a Region from limits of longitude and latitude of the cube.

        Parameters
        ----------
        cube: iris.cube.Cube
            Source cube
        margin: scalar, optional
            Use `margin` number of points or degrees to create a region smaller than the cube
        margin_units: str, optional
            Units of margin. Can be "points" or "degrees".

        Returns
        -------
        aeolus.region.Region
        """
        name = f"extent_of_{cube.name()}"
        lons = cube.coord("longitude").points
        lats = cube.coord("latitude").points
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
        cnstr = iris.Constraint(latitude=lambda x: self.bounds.south <= x <= self.bounds.north)
        if self.bounds.west < self.bounds.east:
            # Western boundary is to the west
            cnstr &= iris.Constraint(longitude=lambda x: self.bounds.west <= x <= self.bounds.east)
        else:
            # Region wrapping around dateline (180deg)
            cnstr &= iris.Constraint(
                longitude=lambda x: (self.bounds.west <= x) or (x <= self.bounds.east)
            )
        return cnstr

    def add_to_ax(self, ax, **kwargs):
        """Add a Rectangle patch to matplotlib axes `ax` with given keyword arguments `kwargs`."""
        from matplotlib.patches import Rectangle  # noqa

        xy = (self.bounds.west, self.bounds.south)
        width = self.lon_size
        height = self.lat_size
        ax.add_patch(Rectangle(xy, width, height, **kwargs))
