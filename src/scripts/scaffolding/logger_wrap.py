#!/usr/bin/env python3
import logging
import time

class UpdatedFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%', start_time=None):
        super().__init__(fmt, datefmt, style)
        self.start_time = start_time or time.time()

    def format(self, record):
        elapsed_seconds = time.time() - self.start_time
        elapsed_time = self.formatTime(record, datefmt="%H:%M:%S")
        elapsed_seconds = record.relativeCreated / 1000
        hours, remainder = divmod(elapsed_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        record.elapsed_seconds = "%02d:%02d:%02d" % (hours, minutes, seconds)
        record.elapsed_time = elapsed_time
        return super().format(record)
    
class UpdatedAdapter(logging.LoggerAdapter):
    def __init__(self, logger, classname):
        super().__init__(logger, {'classname': classname})
    
    def process(self, msg, kwargs):
        return f'{self.extra["classname"]} - {msg}', kwargs


def initLogger(log_file):
    logger = logging.getLogger('HiCPipeline')
    
    # Clear any existing handlers to prevent duplication and unexpected behavior
    logger.handlers.clear()
    
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = UpdatedFormatter('%(elapsed_seconds)s %(levelname)s:  %(message)s', datefmt='%H:%M:%S')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    # Prevent propagation to root logger (which might have its own handlers)
    logger.propagate = False
    
    return logger