"""Input and output functionality."""
import iris


__all__ = ("load_multidir",)


def load_multidir(path_mask, labels, label_name="run"):
    """Load cubelists from multiple directories and merge."""
    joint_cl = iris.cube.CubeList()
    for label in labels:
        cl = iris.load(str(path_mask).format(label))
        for cube in cl:
            cube.attributes["um_version"] = ""  # FIXME
            cube.add_aux_coord(iris.coords.AuxCoord([label], long_name=label_name))
            joint_cl.append(cube)
    return joint_cl.merge()
