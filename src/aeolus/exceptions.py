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
    """Raised when argument type is not recognised or is not allowed."""

    pass


class LoadError(AeolusError):
    """Raised when input files or directories are not found."""

    pass


class BoundaryError(AeolusError):
    """Raised when there is an error with geographical regions."""

    pass


class NotFoundError(AeolusError):
    """Raised when metadata is not found."""

    pass


class UnitFormatError(AeolusError):
    """Raised when cube units cannot be formatted."""

    pass


class MissingCubeError(AeolusError):
    """Raised when cubes required for a calculation are missing."""

    pass


class BadCoordinateError(AeolusError):
    """Raised when coordinates are inconsistent."""

    pass
