# -*- coding: utf-8 -*-
"""Exceptions specific to aeolus package."""


class AeolusWarning(UserWarning):
    """Base class for warnings in aeolus package."""

    pass


class DeprecatedWarning(AeolusWarning):
    """Warning for a deprecated feature."""

    pass


class AeolusError(Exception):
    """Base class for errors in aeolus package."""

    pass


class NotYetImplementedError(AeolusError):
    """
    Raised by missing functionality.

    Different meaning to NotImplementedError, which is for abstract methods.
    """

    pass


class ArgumentError(AeolusError):
    """Raised when argument type is not recognized."""

    pass


class LoadError(AeolusError):
    """Raised when input files or directories are not found."""

    pass


class BoundaryError(AeolusError):
    """Raised when there is an error with geographical regions."""

    pass
