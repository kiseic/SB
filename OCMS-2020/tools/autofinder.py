#!/usr/bin/env python
# -*- Coding: utf-8 -*-

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.


import argparse
import copy
import datetime
import logging
import logging.handlers
import os
import re
import subprocess
import sys
import textwrap
import threading
import time

from concurrent.futures import ThreadPoolExecutor, Future
from queue import Queue, Empty
from string import Template
from typing import List, Dict, TypeVar, Tuple

#----------------------------------#
# Definitions of global parameters #
#----------------------------------#

# NOTE We internally assume Python built-in types for these parameters.
#      Thus, they are initialized according to this assumption.
executable = ''
alpha_x = alpha_y = 1.1
beta_x = beta_y = 0.1
epsilon_x = epsilon_y = 10 ** -5
keep_all = False
unify = True
reuse = False
report_header = ''
report_filename = 'data.resonances'
sorting = 0

# String template for the job output file name
outfile_t = Template('dat.det.cx=$cx.cy=$cy.dx=$dx.dy=$dy')

# String template for the autofinder's output file
header_t = Template(
    '## =======================================\n' +
    '## Resonance data output by autofinder.py\n' +
    '## =======================================\n' +
    '## Command-line for autofinder.py:\n' +
    '## $cmd\n' +
    '## NOTE: missing parameters are set to default\n' +
    '## =======================================\n' +
    '## autofinder.py header ends here\n' +
    '## =======================================\n' +
    '# Resonance data\n' +
    '# polarization=$polarization\n' +
    '# nin=$nin, nout=$nout\n' +
    '# nbe=$nbe\n' +
    '$optional_a' +
    '$optional_b' +
    '# ----------------------------------------\n' +
    '# [Re k] [Im k] [det]\n' +
    '# ----------------------------------------\n'
)


def configure_logger(file_log_level: str = 'INFO',
                     log_directory: str = '/tmp/ocms-logs',
                     log_name: str = 'ocms') -> None:
    """
    Configure the root logger.

    This logger outputs logs into a stream (stderr) and into a file.
    The log levels are set differently: The file logger's level is set to
    INFO by default and is tunable, while the stream logger only outputs
    at the WARNING level.

    By default, log files are stored in <log_directory> and are of the form:
    YYYYMMDD-<log_name>.log, with YYYYY: year, MM: month and DD:day.
    The size of the log files is limited to 10 MB; if more data have to be
    written, a new log file is created with the default name while the name
    of the previous log file is suffixed by an increasing digit,e.g.:
    At time t, the log file 'logfile.log' exceeds 10MB: 'logfile.log' is
    renamed 'logfile.log.0' and a new 'logfile.log' is created, the latter
    always contains the last information.

    Parameters:
    -----------
    file_log_level: str
        Log level used by the file logger

    log_directory: str
        Full-path to the log directory

    log_name: str
        String to append to the log file name

    Ref.:
    -----
    Based on 'duallog' module (https://github.com/acschaefer/duallog)
    by Alexander Schaefer.
    """
    # Create the root logger.
    logger = logging.getLogger()
    logger.setLevel('DEBUG')

    # Create the formatters
    stream_formatter = logging.Formatter(
        '[%(asctime)s] %(funcName)s %(levelname)s -> %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_formatter = logging.Formatter(
        '[%(asctime)s] %(threadName)s %(module)s.%(funcName)s %(levelname)s ' +
        '-> %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Set up logging to the console.
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel('WARNING')
    stream_handler.setFormatter(stream_formatter)
    logger.addHandler(stream_handler)

    # Validate the given directory.
    log_directory = os.path.normpath(log_directory)

    # if output directory is an existing file
    if os.path.isfile(log_directory):
        logger.critical('Output directory is an existing file')
        raise FileExistsError
    # Create a folder for the log files.
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)

    # Log file name
    log_file_format = '{year:04d}{month:02d}{day:02d}-{name}.log'

    # Construct the log file name.
    t = datetime.datetime.now()
    log_file = log_file_format.format(
        year=t.year, month=t.month, day=t.day, name=log_name)
    log_file = os.path.join(log_directory, log_file)

    # Set up logging to the logfile.
    file_handler = logging.handlers.RotatingFileHandler(
        filename=log_file,
        maxBytes=10 * 1024 * 1024,
        backupCount=100
    )
    file_handler.setLevel(file_log_level)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)


def process_inputs(inputs: List[str]) -> Dict:
    """
    Process the input parameters.

    This function transforms a list of string inputs, each input being of the
    form 'key=value' into a dictionary of the form {'key': 'value'}.

    During the process, several checks are performed:
        * Avoid duplications (the first value is kept and used).
        * Avoid symmetry mismatch, e.g. only 'a' is provided instead of 'b' for
          SYM=1.
        * Avoid unexpected values on 'a' and 'b'.

    Parameter:
    ----------
    inputs: list(str)
        List of string of the form 'key=value'

    Returns:
    --------
        Dictionary of the form {'key': 'value'}

    Throw:
    ------
        ValueError
    """
    min_params_set = ('nin', 'nout', 'cx', 'cy', 'dwx', 'dwy', 'dx', 'dy')
    allowed_params = set(min_params_set)
    allowed_params.update('a', 'b')

    input_params = {}
    sym = 0
    ab = 0
    for elt in inputs:
        m = re.search(r'(?P<key>\w+)=(?P<value>.*)', elt)
        if m:
            key = m.group('key')
            value = m.group('value')
            # Process only meaningful parameters and ignore unknown ones
            if key in allowed_params:
                # Avoid duplications
                if not key in input_params.keys():
                    input_params[key] = value
                    # NOTE 'a' and 'b' are special because they are optional
                    #      but can be mandatory depending on SYM value.
                    #      To insure that, we use bitwise operations:
                    #        a  = 0b10 = 2
                    #         b = 0b01 = 1
                    #        ab = 0b11 = 3
                    if key == 'a':
                        sym += 1
                        ab = ab | 2
                    if key == 'b':
                        sym += 1
                        ab = ab | 1
                else:
                    logging.warning(
                        f'Multiple \'{key}\' defined: skip {key}={value} ' +
                        f'and keep {key}={input_params[key]}')
            else:
                logging.warning(
                    f'Skip unknown parameter \'{key}={value}\'')
    # Check consistency on the parity vs. symmetry
    # NOTE SYM=4 has the same validation conditions as SYM=2 (both 'a' and 'b'
    #      defined).
    if not ab & sym and ab | sym > 0:
        raise ValueError(
            'Symmetry mismatch: expect SYM=0: no \'a\', no \'b\'; ' +
            'SYM=1: only \'b\'; SYM=2 or SYM=4: \'a\' and \'b\'')

    # Set 'a' and 'b' to default if not provided/needed
    if sym <= 1:
        input_params['a'] = ''
    if sym == 0:
        input_params['b'] = ''

    # Check that the provided 'a' and 'b' are -1 or 1
    to_check = []
    if ab & 2:
        to_check.append('a')
    if ab & 1:
        to_check.append('b')
    for key in to_check:
        if not JobParams.to_python_int(input_params[key]) in (-1, 1):
            raise ValueError(
                f'Provided \'{key}\' value is not -1 or 1: ' +
                f'\'{key}={input_params[key]}\'')

    # Stop if not all the mandatory parameters are not provided
    for param in min_params_set:
        if not param in input_params.keys():
            raise ValueError(f'Parameter \'{param}\' is not provided')

    return input_params


def generate_report_params(cmd: str, inputs: Dict, basename: str) \
        -> Tuple[str, str]:
    """
    Generate the output file parameters.

    This method generates the report header using the input parameters and the
    executable name, and creates the report file name by appending parity
    characters if needed ('a' first, then 'b').

    Possible returned filename values:
      * no 'a', no 'b': filename = basename
      * only 'b': filename = basename.e (b = 1) or
                             basename.o (b = -1)
      * 'a' and 'b': filename = basename.eo (a =  1, b = -1) or
                                basename.oe (a = -1, b =  1) or
                                basename.ee (a =  1, b =  1) or
                                basename.oo (a = -1, b = -1)


    Parameters:
    -----------
    cmd: str
        Raw typed command-line
    inputs: dict
        Processed inputs

    Returns:
    --------
    header, filename: tuple(str, str)
        Formatted header and processed file name
    """
    global executable, header_t

    header = ''

    # Extract information from the executable name
    polarization = ''
    nbe = ''
    matched = re.match(r'.*det.*\.(?P<polarization>\w+)\.NBE\=(?P<nbe>\d+)',
                       executable)
    if matched:
        polarization = matched.group('polarization')
        nbe = matched.group('nbe')
    else:
        logging.warning('Cannot extract polarization and NBE information ' +
                        'from the executable name')

    optionals = {}
    optionals['a'] = ''
    optionals['b'] = ''
    suffix = ''
    for key in ('a', 'b'):
        # NOTE: 'a' and 'b' are always present in the 'inputs' array
        if inputs[key]:
            optionals[key] = f'# {key}={inputs[key]}\n'
            suffix += 'e' if JobParams.to_python_int(inputs[key]) == 1 else 'o'
    if suffix:
        suffix = '.' + suffix

    # Wrap the original command to ease the readability while being usable as
    # is (e.g. using copy/paste and removing the leading '#' symbols)
    wrapped_cmd = textwrap.wrap(cmd, width=42)
    glue = ' \\\n##   '
    glued_cmd = glue.join(wrapped_cmd)

    header_params = {'cmd': glued_cmd, 'polarization': polarization, 'nbe': nbe,
                     'optional_a': optionals['a'],
                     'optional_b': optionals['b']}
    header_params.update(inputs)

    header = header_t.safe_substitute(header_params)

    filename = basename + suffix

    return header, filename


IF = TypeVar('IF', int, float)


class JobParams(object):
    """
    JobParams class.

    This class stores the job parameters as strings compliant with Fortran
    formats, i.e. floats are stored in strings of the form '1.0d0'. Such a
    string is a Fortran type-compliant string, or simpy a 'Fortran' string.

    This class then provides methods for performing basic transformations such
    as casting a 'Fortran' string to a Python float, or vice versa.
    """
    precision: int = 9

    def __init__(self, nin: str, nout: str, a: str, b: str, cx: str, cy: str,
                 dwx: str, dwy: str, dx: str, dy: str) -> None:
        self.nin = str(nin)
        self.nout = str(nout)
        self.a = str(a)
        self.b = str(b)
        self.cx = str(cx)
        self.cy = str(cy)
        self.dwx = str(dwx)
        self.dwy = str(dwy)
        self.dx = str(dx)
        self.dy = str(dy)

    def __str__(self):
        return f'{self.__dict__}'

    def dump(self) -> Dict:
        """Dump the internal parameters into a dictionary."""
        return self.__dict__

    @staticmethod
    def to_fortran_string(variable: IF) -> str:
        """
        Cast a Python variable to a 'Fortran' string.

        Parameter:
        ----------
        variable: int or float
            Variable to cast

        Returns:
        --------
            'Fortran' string

        Examples:
        ---------
            * JobParams.to_fortran_string(4.2) -> '4.2d0'
            * JobParams.to_fortran_string(42) -> '42'

        Throws:
        -------
            ValueError  
        """
        precision = JobParams.precision
        if isinstance(variable, int):
            return str(variable)
        elif isinstance(variable, float):
            str_variable = f'{round(variable, precision)}'
            decimal, exponent = JobParams.filter_exponent_chars(str_variable)
            return f'{decimal}d{exponent}'
        else:
            raise ValueError(f'Unsupported type: {variable}: {type(variable)}')

    @staticmethod
    def to_fortran_string_from_vars(*variables: List[IF]) -> List[str]:
        """
        Cast a Python variables to a 'Fortran'strings.

        Parameter:
        ----------
        variables: list(int) or list(float)
            List of variables to cast

        Returns:
        --------
            List of 'Fortran' strings

        Example:
        --------
            * JobParams.to_fortran_string_from_vars(4.2, 42) -> ['4.2d0', '42']
        """
        strings = []
        for variable in variables:
            strings.append(JobParams.to_fortran_string(variable))
        return [s for s in strings]

    @staticmethod
    def to_python_int(variable: str) -> int:
        """
        Cast a 'Fortran'-compliant string to a Python integer.

        Parameter:
        ----------
        variable: str
            String to cast

        Returns:
        --------
            Python integer

        Example:
        --------
            * JobParams.to_python_int('42') -> 42
        """
        return int(variable)

    @staticmethod
    def to_python_float(variable: str) -> float:
        """
        Cast a 'Fortran'-compliant string to a Python float.

        Parameter:
        ----------
        variable: str
            String to cast

        Returns:
        --------
            Python float

        Example:
        --------
            * JobParams.to_python_int('4.2d0') -> 4.2
        """
        precision = JobParams.precision
        decimal, exponent = JobParams.filter_exponent_chars(variable)
        variable = float(decimal) * 10.0 ** int(exponent)
        return round(variable, precision)

    @staticmethod
    def filter_exponent_chars(variable: str) -> Tuple[float, float]:
        """
        Filter the exponent characters from a string representing a float.

        Parameter:
        ----------
        variable: str
            String representing a Fortran float

        Returns:
        --------
        decimal, exponent: tuple(float, float)
            Tuple consisting of the decimal part and exponent part of the float

        Example:
        --------
            * JobParams.filter_exponent_chars('4.2d-2') -> (4.2, -2) 

        Throws:
        -------
            ValueError   
        """
        pattern = r'(?P<decimal>[+-]?\d+[.,]?\d*)(?P<char>[ed])?(?P<exponent>[+-]?\d*)?'
        matched = re.fullmatch(pattern, variable, flags=re.IGNORECASE)
        if matched:
            mismatch_cne = matched.group('char') and \
                not matched.group('exponent')
            mismatch_nce = not matched.group('char') and \
                matched.group('exponent')
            mismatch = mismatch_cne or mismatch_nce
            if not mismatch:
                decimal = matched.group('decimal')
                exponent = '0' if not matched.group('exponent') \
                    else matched.group('exponent')
                return decimal, exponent

        raise ValueError(f'Unsupported pattern: {variable}')


class Job(object):
    """
    Job class.

    This class is in charge of all the complex computations:
        * Run the executable witht its internal parameters
        * Find the minimum/minima
    """
    ordered_params = ('nin', 'nout', 'a', 'b', 'cx', 'cy', 'dwx', 'dwy', 'dx',
                      'dy')
    # Global lock to access/modify the other class parameters
    rlock = threading.RLock()
    # Record of the processed files
    processed_output_files = set()
    # File locks for thread-safe R/W operations
    file_locks = {}

    def __init__(self, job_params: JobParams) -> None:
        self.job_params = job_params

    def run(self) -> Tuple[List[JobParams], List[str]]:
        """
        Run the job to do the hard-work.

        Its workflow is:
            1. Run the executable
            2. Find the minimum/minima (i.e. resonance(s))
            3. Compute the new job parameters using step 2 information
            4. Report the new job parameters and the minimum/minima
        """
        global executable, epsilon_x, epsilon_y, keep_all, reuse, outfile_t

        logging.info(f'Run job with params: {self.job_params}')

        string_input = ''
        for param_name in Job.ordered_params:
            string_input += f'{getattr(self.job_params, param_name)} '
        string_input = ' '.join(string_input.split())
        logging.debug(f'string_input={string_input}')

        # Get the binary input string (needed by subprocess.run)
        binary_input = string_input.encode('ascii')
        logging.debug(f'binary_input={binary_input}')

        outfile = outfile_t.substitute(cx=self.job_params.cx,
                                       cy=self.job_params.cy,
                                       dx=self.job_params.dx,
                                       dy=self.job_params.dy)

        processed_output_files = set()
        file_lock = None
        with Job.rlock:
            # NOTE processed_output_files class attribute is needed because
            #      we cannot use copy.deepcopy(Job.file_locks.keys()) which
            #      provides the same information. We can however make a deep
            #      copy of Job.file_locks but then we also copy the locks
            #      which might be troublesome
            processed_output_files = copy.deepcopy(Job.processed_output_files)
            logging.debug(f'processed_output_files={processed_output_files}')

            if not outfile in processed_output_files:
                # Add it to the processed files to avoid concurrency
                Job.processed_output_files.add(outfile)
                Job.file_locks[outfile] = threading.RLock()
            file_lock = Job.file_locks[outfile]

        logging.debug(f'acquire file lock: {outfile}')
        file_lock.acquire()
        if not outfile in processed_output_files:
            # Create the output file (or read it)
            if not os.path.exists(outfile) or not reuse:
                logging.info(f'Generate {outfile}')
                with open(outfile, 'w') as of:
                    # Compute the determinant value distribution
                    subprocess.run([executable],
                                   input=binary_input,
                                   stdout=of)
            else:
                logging.info(f'File {outfile} exists: read it')
        else:
            logging.debug(f'release file lock: {outfile}')
            file_lock.release()
            # Cancel the job because it has already been processed
            return [], []

        logging.info(f'Find the minima from {outfile}')
        cmd = ['minfinder.py', outfile, '--nodisplay']
        # Find the minima
        complete_process = subprocess.run(cmd, stdout=subprocess.PIPE)
        binary_minima = complete_process.stdout
        logging.debug(f'binary_minima={binary_minima}')

        # Determine the new job parameters
        new_job_params = []
        minima = []
        to_remove = not copy.deepcopy(keep_all)
        # NOTE The regular expression naturally filters the case where no
        #      no minima are found.
        for matching in re.finditer(rb'(.*?)\s(.*?)\s(.*?)\s', binary_minima):
            logging.debug(f'matching={matching.groups()}')
            cx, cy, scaled_det = [bin_value.decode('ascii')
                                  for bin_value in matching.groups()]
            minima.append((cx, cy, scaled_det))
            # Check the stop conditions
            if JobParams.to_python_float(self.job_params.dx) > epsilon_x and \
                    JobParams.to_python_float(self.job_params.dy) > epsilon_y:
                job_params = copy.deepcopy(self.job_params)
                job_params.cx, job_params.cy = \
                    JobParams.to_fortran_string_from_vars(
                        JobParams.to_python_float(cx),
                        JobParams.to_python_float(cy)
                    )
                job_params = Job._update_job_params(job_params)
                new_job_params.append(job_params)
                if to_remove:
                    logging.debug(f'remove {outfile}')
                    os.remove(outfile)
                    to_remove = False
        logging.debug(f'release file lock: {outfile}')
        file_lock.release()

        logging.debug(
            f'new_job_params: {[p.dump() for p in new_job_params]}')
        logging.info(f'Minima: {minima}')

        return new_job_params, minima

    @staticmethod
    def _apply_linear_transformation(coef: float, value: str) -> float:
        """
        Return f(x) = a * x.

        Parameters:
        -----------
        coef: float
            Multiplicative coefficient

        value: str
            Input value

        Returns:
        --------
            coef * value
        """
        return coef * JobParams.to_python_float(value)

    @staticmethod
    def _update_job_params(params: JobParams) -> JobParams:
        """
        Update the input job parameters to obtain the new job parameters.

        This method updates dwx, dwy, dx and dy according to the rule:
            x_i+1 = K_x * y_i
        where i represents the iteration, x = {dwx, dwy, dx, dy}, y is {dx, dy}
        and K_x is the multiplicative coefficient associated to x (alpha for
        dwx and dwy, beta for dx, dy).

        Parameters:
        -----------
        params: JobParams
            Job parameters to update.

        Returns:
        --------
        params: JobParams
            Updated job parameters
        """
        global alpha_x, alpha_y, beta_x, beta_y, epsilon_x, epsilon_y, unify

        dwx = Job._apply_linear_transformation(alpha_x, params.dx)
        dwy = Job._apply_linear_transformation(alpha_y, params.dy)

        dx = JobParams.to_python_float(params.dx)
        if dx > epsilon_x:
            dx = Job._apply_linear_transformation(beta_x, params.dx)
        dy = JobParams.to_python_float(params.dy)
        if dy > epsilon_y:
            dy = Job._apply_linear_transformation(beta_y, params.dy)

        if unify:
            # NOTE Assume dx = dy after the first iteration
            dx = min(dx, dy)
            dy = dx

        params.dwx, params.dwy, params.dx, params.dy = \
            JobParams.to_fortran_string_from_vars(dwx, dwy, dx, dy)

        return params


class AutoFinder(object):
    """
    AutoFinder class.

    This is the main class that triggers the jobs, monitors them and processes
    the results that are asynchroneously received.

    The main thread monitors the overall progress. Once any arbitrary number of
    jobs is pushed into the job queue, it triggers the asynchroneous job
    handling. This task is devoted to an auxiliary thread that controls the
    workers pool. The job results are asynchroneously processed and displayed.

    The overall procedure ends when all the enqueued jobs have been processed.
    """

    def __init__(self, workers: int = 5) -> None:
        threading.current_thread().name = self.__class__.__name__
        self.rlock = threading.RLock()
        self.running = False
        self.job_queue = Queue()
        self.enqueued = 0
        self.last_enqueued = 0
        self.done = 0
        self.writer = ThreadPoolExecutor(max_workers=1)
        self.executor = ThreadPoolExecutor(max_workers=workers,
                                           thread_name_prefix='Job')
        self.assign_trigger = False
        self.job_hdlr_thread = threading.Thread(name='JobHandler',
                                                target=self._handle_jobs)
        self.minima = set()

    def start(self) -> None:
        """Start the auto-finding procedure"""
        logging.info(f'Start auto-finding procedure')
        self.running = True

        self.job_hdlr_thread.daemon = True
        logging.debug(f'start daemon: {self.job_hdlr_thread}')
        self.job_hdlr_thread.start()

        self.monitor()

    def stop(self) -> None:
        """Start the auto-finding procedure"""
        logging.info(f'Stop auto-finding procedure')
        self.running = False

        logging.debug(f'stop daemon {self.job_hdlr_thread}')
        self.job_hdlr_thread.join()

        logging.debug(f'stop job queue')
        self.job_queue.join()

        logging.debug(f'shutdown the pool executor')
        self.executor.shutdown()

        logging.debug(f'shutdown the writer')
        self.writer.shutdown()

    def _update_enqueued(self, val: int = 1) -> None:
        """
        Update the number of enqueued jobs.

        This method is thread-safe.

        Parameter:
        ----------
        val: int
            Value to add to the number of enqueued jobs
        """
        with self.rlock:
            logging.debug(f'enqueued={self.enqueued}, val={val}')
            self.enqueued += val

    def _get_enqueued(self) -> int:
        """Return the number of enqueued jobs (thread-safe)."""
        with self.rlock:
            logging.debug(f'enqueued={self.enqueued}')
            return self.enqueued

    def _update_done(self, val: int = 1) -> None:
        """
        Update the number of processed jobs.

        This method is thread-safe.

        Parameter:
        ----------
        val: int
            Value to add to the number of processed jobs
        """
        with self.rlock:
            logging.debug(f'done={self.done}, val={val}')
            self.done += val

    def _get_done(self):
        """Return the number of processed jobs (thread-safe)."""
        with self.rlock:
            logging.debug(f'done={self.done}')
            return self.done

    def _assign(self, assign: bool) -> None:
        """
        Set the asynchroneous "assign" trigger.

        This method is thread-safe.

        Parameter:
        ----------
        assign: bool
            Trigger value
        """
        with self.rlock:
            logging.debug(f'assign={assign}')
            self.assign_trigger = assign

    def enqueue(self, job_params: JobParams) -> None:
        """
        Enqueue a job.

        Parameter:
        ----------
        job_params: JobParams
            Job parameters
        """
        logging.debug(f'job_params={job_params}')
        self._update_enqueued()
        self.job_queue.put(job_params)

    def monitor(self) -> None:
        """
        Monitor the overall progress.

        This method triggers the job assignment task whenever one or more jobs
        are added. Otherwise it periodically checks the overall progress.
        The whole process ends when all enqueud jobs have been processed. 
        """
        while True:
            enqueued = self._get_enqueued()
            done = self._get_done()
            if enqueued == done:
                break
            # Assign the new jobs
            if self.last_enqueued < enqueued:
                self.last_enqueued = enqueued
                self._assign(True)
            time.sleep(0.1)
        self.stop()

    def _handle_jobs(self) -> None:
        """
        Handle the jobs.

        This method is triggered by the main thread to dispatch the jobs
        between the workers.
        """
        while self.running:
            if self.assign_trigger:
                self._assign_jobs()
                self._assign(False)
            else:
                time.sleep(0.1)

    def _assign_jobs(self) -> None:
        """
        Assign jobs to the workers.

        This method transfers all the available jobs in the queue to the
        workers. The jobs start asynchroneously. It stops when the queue is
        empty.
        """
        while True:
            try:
                job_params = self.job_queue.get_nowait()
            except Empty:
                logging.debug('no jobs')
                break
            job = Job(job_params)
            with self.rlock:
                logging.debug(f'job_params={job_params}')
                future = self.executor.submit(job.run)
                future.add_done_callback(self._process_future)

    def _process_future(self, future: Future) -> None:
        """
        Process the job result.

        This method is triggered when a job result (a future) is available. 

        Parameter:
        ----------
        future: Future
            Job result
        """
        job_params_list, minima = future.result()
        logging.debug(
            f'job_params_list={[p.dump() for p in job_params_list]}')
        logging.debug(f'minima={minima}')

        # Assign new jobs (if any) or display the minima (if any)
        if job_params_list:
            for job_params in job_params_list:
                self.job_queue.put(job_params)
            # Update the number of enqueued jobs only once
            self._update_enqueued(len(job_params_list))
        elif minima:
            self._report_minima(minima)
        else:
            logging.debug('canceled job')

        self.job_queue.task_done()
        self._update_done()

    def _write_report_file(self, minima: List[Tuple[float, float, float]]) \
            -> None:
        """
        Write the output report file.

        The minima are sorted before being written.

        Parameters:
        -----------
        minima: list(tuple(float, float, float))
            Unsorted list of minima
        """
        global report_header, report_filename, sorting

        sorted_minima = sorted(minima, key=lambda minimum: minimum[sorting])

        with open(report_filename, 'w') as of:
            of.write(report_header)
            for cx, cy, scaled_det in sorted_minima:
                of.write(f'{cx} {cy} {scaled_det}\n')

    def _report_minima(self, minima: List[Tuple[float, float, float]]) -> None:
        """
        Report the found minima.

        This methods stores the found minima and triggers the writing into the
        output file.

        Parameter:
        ----------
        minima: list(tuple(float, float, float))
            List of found resonances, each of them consisting of the following
            characteristics:
              (real part, imaginary part, scaled determinant value) 
        """
        for cx, cy, scaled_det in minima:
            minimum = (cx, cy, scaled_det)
            with self.rlock:
                if not minimum in self.minima:
                    self.minima.add(minimum)
                current_minima = self.minima
        # NOTE No specific cautious is needed because there a single writer
        #      and the associated executor has its own internal thread-safe
        #      FIFO queue
        self.writer.submit(self._write_report_file, minima=current_minima)


def main(argv):
    global executable, alpha_x, alpha_y, beta_x, beta_y, epsilon_x, epsilon_y
    global keep_all, unify, reuse, report_header, report_filename, sorting

    parser = argparse.ArgumentParser(
        description=str('Tool for automatically finding minima of a ' +
                        'determinant value distribution in the complex ' +
                        'wave number space.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('executable', type=str,
                        help='Full-path name to the executable to run')
    parser.add_argument('inputs', type=str, nargs='*',
                        help='Input parameters of the form \
                            \'nin=3.3d0 nout=1.0d0 a=1 b=-1 ...\'')
    parser.add_argument('--workers', '-w', type=int, default=5,
                        help='Number of concurrent workers')
    parser.add_argument('--alphax', '-ax', type=str, default='1.1d0',
                        help='Multiplicative coefficient for determining new \
                            dwx')
    parser.add_argument('--alphay', '-ay', type=str, default='1.1d0',
                        help='Multiplicative coefficient for determining new \
                            dwy')
    parser.add_argument('--betax', '-bx', type=str, default='0.1d0',
                        help='Multiplicative coefficient for determining new \
                            dx')
    parser.add_argument('--betay', '-by', type=str, default='0.1d0',
                        help='Multiplicative coefficient for determining new \
                            dy')
    parser.add_argument('--epsilonx', '-ex', type=str, default='1.0d-5',
                        help='Stop criteria on dx')
    parser.add_argument('--epsilony', '-ey', type=str, default='1.0d-5',
                        help='Stop criteria on dy')
    parser.add_argument('--output', '-o', type=str, default='data.resonances',
                        help='Base output file name that might be suffixed \
                            with parity characters e/o (suffix format: \'.b\' \
                            for SYM=1, \'.ab\' for SYM=2)')
    parser.add_argument('--sorting', '-s', type=int, choices=(0, 1, 2), default=0,
                        help='Field index to used for sorting the resonances \
                        (ascending order): 0=Re(k), 1=Im(k), 2=Det')
    parser.add_argument('--ununify', action='store_true',
                        help='Allow dx and dy to be set independently')
    parser.add_argument('--reuse', action='store_true',
                        help='Allow to reuse (i.e. read) existing files')
    parser.add_argument('--keepall', action='store_true',
                        help='Keep all intermediate files')
    parser.add_argument('--loglevel', type=str, default='INFO',
                        help='Set the log level of the log file')
    parser.add_argument('--logdir', type=str, default='/tmp/ocms-logs',
                        help='Set the log directory')
    parser.add_argument('--logname', type=str, default='ocms',
                        help='Set the string to append to the log file name \
                        (pattern: YYYYMMDD-<logname>.log)')
    parser.add_argument('--cleanlogs', action='store_true',
                        help='Erase all log files in the log directory BEFORE \
                            running the process. Make sure to make a copy of \
                            the log files you want to keep before using this \
                            command')
    args = parser.parse_args()

    # Record the command-line to run autofinder.py
    raw_cmd = sys.argv
    raw_cmd[0] = os.path.basename(raw_cmd[0])
    cmd = ' '.join(raw_cmd)

    # Configure the logger
    if args.cleanlogs:
        cmd = ['rm', '-rf', args.logdir]
        subprocess.run(cmd)
    configure_logger(file_log_level=args.loglevel.upper(),
                     log_directory=args.logdir, log_name=args.logname)

    # Assign the mandatory parameters
    executable = args.executable
    if not os.path.isfile(executable):
        logging.critical(f'\'{executable}\' does not exist')
        raise SystemExit(1)

    try:
        inputs = process_inputs(args.inputs)
    except ValueError as e:
        logging.critical(e)
        raise SystemExit(1)

    # Prepare the output file parameters
    report_header, report_filename = generate_report_params(cmd, inputs,
                                                            args.output)

    # Assign the optional parameters
    workers = args.workers
    alpha_x = JobParams.to_python_float(args.alphax)
    alpha_y = JobParams.to_python_float(args.alphay)
    beta_x = JobParams.to_python_float(args.betax)
    beta_y = JobParams.to_python_float(args.betay)
    epsilon_x = JobParams.to_python_float(args.epsilonx)
    epsilon_y = JobParams.to_python_float(args.epsilony)
    keep_all = args.keepall
    unify = not args.ununify
    reuse = args.reuse
    sorting = args.sorting

    # Start the auto-finding procedure
    master = AutoFinder(workers)
    job_params = JobParams(**inputs)
    master.enqueue(job_params)
    try:
        master.start()
    except KeyboardInterrupt:
        sys.stderr.write("\nKeyboard interrupted.\n")
        master.stop()
        raise SystemExit(1)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
