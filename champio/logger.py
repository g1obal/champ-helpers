"""
Logger initializer

Author: Gokhan Oztarhan
Created date: 06/12/2021
Last modified: 04/12/2022
"""

import sys
import logging


def set_logger(
    verbose_file=1, verbose_console=1, 
    filename='logfile', filemode='w', terminator=''
):
    log_formatter = logging.Formatter('%(message)s')
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    if verbose_file:
        file_handler = logging.FileHandler(filename, mode=filemode)
        file_handler.terminator = terminator
        file_handler.setFormatter(log_formatter)
        file_handler.setLevel(logging.INFO)
        logger.addHandler(file_handler)
    
    if verbose_console:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.terminator = terminator
        stream_handler.setFormatter(log_formatter)
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)


def unset_logger():
    logger = logging.getLogger()
    
    for handler in logger.handlers:
        handler.close()
        logger.removeHandler(handler)

    logging.shutdown()

