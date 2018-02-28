import logging

FORMAT = '%(asctime)s %(filename)-15s %(levelname)-10s: %(message)s'
logging.basicConfig(format=FORMAT)

def get_logger(verbosity=logging.INFO):
    logger = logging.getLogger(__name__)
    logger.setLevel(verbosity)
    return logger

