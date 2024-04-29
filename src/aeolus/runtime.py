"""Runtime configuration."""
import threading

__all__ = ("RUNTIME",)


class RuntimeOpts(threading.local):
    """Run-time configuration controller."""

    def __init__(self, figsave_stamp=False, figsave_reduce_size=False):
        """Initialise run-time configuration."""
        self.__dict__["figsave_stamp"] = figsave_stamp
        self.__dict__["figsave_reduce_size"] = figsave_reduce_size

    def __repr__(self):
        """Repr of run-time configuration."""
        keyval_pairs = [f"{key}={self.__dict__[key]}" for key in self.__dict__]
        msg = f"RuntimeOpts({', '.join(keyval_pairs)})"
        return msg

    def __setattr__(self, name, value):
        """Set attributes."""
        if name not in self.__dict__:
            msg = f"'RuntimeOpts' object has no attribute {name!r}"
            raise AttributeError(msg)
        self.__dict__[name] = value


RUNTIME = RuntimeOpts()
