import logging

'''
Convenient logging.

HANDL package exposes a `get_logger()` to allow users to convieniently
log and print with pre-formated logging statements.
'''

# Logging format
FORMAT = '%(asctime)s %(filename)-15s %(levelname)-10s: %(message)s'
logging.basicConfig(format=FORMAT)

def get_logger(verbosity=logging.INFO):
    '''
    Returns logger object
    '''
    logger = logging.getLogger(__name__)
    logger.setLevel(verbosity)
    return logger
