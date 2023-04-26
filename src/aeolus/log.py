# -*- coding: utf-8 -*-
"""Logging utilities."""
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, Union

import loguru
from loguru import logger


__all__ = ("create_logger",)


def create_logger(
    script_path: Path, subdir: Optional[Union[str, Path]] = "logs"
) -> loguru._logger.Logger:
    """Create a logger using loguru."""
    logpath = script_path.parent / subdir
    logpath.mkdir(exist_ok=True)
    logger.configure(
        handlers=[
            {"sink": sys.stdout, "level": "INFO"},
            {
                "sink": logpath / f"log_{script_path.stem}_{datetime.now():%Y-%m-%d_%H%M%S}.log",
                "level": "DEBUG",
            },
        ]
    )
    return logger
