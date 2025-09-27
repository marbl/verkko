#!/usr/bin/env python3
import sys
import logging
import time
import os
import subprocess


def setup_logging(log_file, file_log_level = logging.DEBUG):
    log_level = logging.INFO
    # Configure logging with runtime from program start
    start_time = time.time()
    
    class RuntimeFormatter(logging.Formatter):
        def __init__(self, fmt=None, datefmt=None, start_time=None):
            super().__init__(fmt, datefmt)
            self.start_time = start_time
            
        def format(self, record):
            record.runtime = time.time() - self.start_time
            return super().format(record)
        
    datefmt = '%H:%M:%S'

    def format_runtime(seconds):
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        seconds = int(seconds % 60)
        return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    
    class RuntimeFormatter(logging.Formatter):
        def __init__(self, fmt=None, datefmt=None, start_time=None):
            super().__init__(fmt, datefmt)
            self.start_time = start_time
            
        def format(self, record):
            runtime_seconds = time.time() - self.start_time
            record.runtime = format_runtime(runtime_seconds)
            return super().format(record)

    log_format = '%(runtime)s - %(levelname)s - [%(filename)s:%(funcName)s:%(lineno)d] - %(message)s'
    formatter = RuntimeFormatter(log_format, datefmt)
        
    # Always log to both file and console
    #log_file = os.path.join(args.outdir, f"{args.basename}.log")

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    # File handler
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(RuntimeFormatter(log_format, datefmt=datefmt, start_time=start_time))
    file_handler.setLevel(file_log_level)
    root_logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(RuntimeFormatter(log_format, datefmt=datefmt, start_time=start_time))
    console_handler.setLevel(logging.INFO)
    root_logger.addHandler(console_handler)
    
    # Log the GitHub commit hash if available
    try:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=sys.path[0]).strip().decode('utf-8')
        logging.info(f"GitHub Commit Hash: {commit_hash}")
    except Exception as e:
        logging.warning(f"Failed to retrieve GitHub commit hash: {e}")

    # Log the command-line arguments
    logging.info(f"Command-line arguments: {' '.join(sys.argv)}")
    logging.info(f"Logging to file: {log_file}")


def log_assert(condition, message, logger=None):
    """Assert a condition and log an error message if it fails."""
    if not condition:
        error_msg = f"Assertion failed: {message}"
        if logger:
            logger.error(error_msg)
        else:
            logging.error(error_msg)
        exit(1)
        #raise AssertionError(error_msg)
