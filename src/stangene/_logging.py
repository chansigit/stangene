"""Structured logging for stangene."""

import logging


def get_logger(name: str) -> logging.Logger:
    """Get a logger under the stangene namespace."""
    logger = logging.getLogger(f"stangene.{name}")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(name)s] %(levelname)s: %(message)s")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger
