import logging
from rolypoly.utils.logging.loggit import get_logger, setup_logging


def caller_func():
    logger = get_logger()
    logger.info("Log from caller_func")


if __name__ == "__main__":
    # Configure logging for this test using the caller's logger name
    setup_logging(log_file=None, log_level="info", logger_name="__main__")
    caller_func()
