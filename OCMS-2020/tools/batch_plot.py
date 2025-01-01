#!/usr/bin/env python
# -*- Coding: utf-8 -*-

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.


import os
import re
import sys
import shlex
import pathlib
import operator
import argparse
import subprocess
import multiprocessing
from typing import List, Dict, TypeVar, Tuple

APPLICATION_NAME    = "Batch Plot"
APPLICATION_VERSION = "1.1.0"


# ================================================================================
# - Definitions of file/directorory/parameter names, commands, etc. -
# ================================================================================

BATCH_PROCESSING_DIRECTORY_NAME_PREFIX = "mode{}k="  # parity code (.ee., etc) will be embedded in {}. 
BATCH_PROCESSING_DIRECTORY_NAME_PARITY = { "1": "e", "-1": "o" }

PLOT2D_RELATIVE_PATH_FROM_BATCH_PLOT    = "plot2d.py"
NEARFIELD_RELATIVE_PATH_FROM_BATCH_PLOT = "nearfield.plt"
FARFIELD_RELATIVE_PATH_FROM_BATCH_PLOT  = "farfield.plt"

CONFIG_PARAM_NAME_SYM                    = "SYM"
CONFIG_PARAM_NAME_PARITY_A               = "a"
CONFIG_PARAM_NAME_PARITY_B               = "b"
CONFIG_PARAM_NAME_HUSIMI_EXE_FILE        = "husimi_exe"
CONFIG_PARAM_NAME_WFUNC_EXE_FILE         = "wfunc_exe"
CONFIG_PARAM_NAME_BOUNDARY_DATA_FILE     = "boundary_data"
CONFIG_PARAM_NAME_HUSIMI_PLOT_OPTION = "husimi_plot_options"
CONFIG_PARAM_NAME_WFUNC_PLOT_OPTION  = "wfunc_plot_options"
CONFIG_PARAM_NAME_WFUNC_LOGPLOT_OPTION   = "wfunc_logplot_options"

SHELL_SCRIPT_SHEBANG              = "#!/bin/bash -f"
SHELL_SCRIPT_NAME_HUSIMI          = "husimi.sh"
SHELL_SCRIPT_NAME_WFUNC           = "wfunc.sh"
SHELL_SCRIPT_PRIORITY_HUSIMI      = 2  # assign smaller number to give higher priority
SHELL_SCRIPT_PRIORITY_WFUNC       = 1  # (higher one will be computed before lower one)
SHELL_SCRIPT_PARAM_NAME_EXECFILE  = "execfile"
SHELL_SCRIPT_PARAM_NAME_REAL_K    = "kx"
SHELL_SCRIPT_PARAM_NAME_IMAG_K    = "ky"

# List of names of file-specifying-parameters to be resolved from relative-paths to absolute-paths.
CONFIG_PARAM_NAMES_OF_FILES = [ 
    CONFIG_PARAM_NAME_HUSIMI_EXE_FILE, 
    CONFIG_PARAM_NAME_WFUNC_EXE_FILE, 
    CONFIG_PARAM_NAME_BOUNDARY_DATA_FILE
]

# List of names of parameters defined in the config file which should not be written in shell scripts.
SHELL_SCRIPT_PARAM_NAMES_EXCLUDED = [
    CONFIG_PARAM_NAME_HUSIMI_EXE_FILE,
    CONFIG_PARAM_NAME_WFUNC_EXE_FILE,
    CONFIG_PARAM_NAME_WFUNC_PLOT_OPTION,
    CONFIG_PARAM_NAME_WFUNC_LOGPLOT_OPTION,
    CONFIG_PARAM_NAME_HUSIMI_PLOT_OPTION
]

def get_exec_command_in_wfunc_shell_script(symmetry_code: str) -> str:
    """
    Returns the content of the command to invoke executable program in wfunc.sh
    """
    if symmetry_code == "0":
        return "echo $nin $nout $kx $ky $xmin $xmax $ymin $ymax $ixmax $iymax | $execfile"
    elif symmetry_code == "1":
        return "echo $nin $nout $b $kx $ky $xmin $xmax $ymin $ymax $ixmax $iymax | $execfile"
    elif symmetry_code == "2":
        return "echo $nin $nout $a $b $kx $ky $xmin $xmax $ymin $ymax $ixmax $iymax | $execfile"
    elif symmetry_code == "4":
        return "echo $nin $nout $a $b $kx $ky $xmin $xmax $ymin $ymax $ixmax $iymax | $execfile"
    else:
        raise ValueError("Unexpected symmetry code: " + symmetry_code)

def get_exec_command_in_husimi_shell_script(symmetry_code: str) -> str:
    """
    Returns the content of the command to invoke executable program in husimi.sh
    """
    if symmetry_code == "0":
        return "echo $nin $nout $kx $ky | $execfile"
    elif symmetry_code == "1":
        return "echo $nin $nout $b $kx $ky | $execfile"
    elif symmetry_code == "2":
        return "echo $nin $nout $a $b $kx $ky | $execfile"
    elif symmetry_code == "4":
        return "echo $nin $nout $a $b $kx $ky | $execfile"
    else:
        raise ValueError("Unexpected symmetry code: " + symmetry_code)

def get_plot_command_in_wfunc_shell_script(plot2d_file_path: str, plot_option: str, logplot_option: str, \
        neadfield_file_path: str, farfield_file_path: str) -> str:
    """
    Returns the content of the command to call plot2d.py, nearfield.plt, and farfield.plt in wfunc.sh
    """
    plot2d_file_path = shlex.quote(plot2d_file_path)
    nearfield_file_path = shlex.quote(neadfield_file_path)
    farfield_file_path = shlex.quote(farfield_file_path)
    return plot2d_file_path + " dat.wfunc --vscale log --savefig " + logplot_option + "\n" \
         + "mv dat.wfunc.png dat.wfunc_logscale.png\n" \
         + plot2d_file_path + " dat.wfunc --savefig " + plot_option + "\n" \
         + nearfield_file_path + "\n" \
         + farfield_file_path

def get_plot_command_in_husimi_shell_script(plot2d_file_path: str, plot_option: str) -> str:
    """
    Returns the content of the command to calling plot2d.py in husimi.sh
    """
    plot2d_file_path = shlex.quote(plot2d_file_path)
    return plot2d_file_path + " dat.husimi --savefig " + plot_option

def convert_parity_value_to_parity_name(parity_value: str) -> str:
    """
    Converts a parity value ("1" or "-1") to a parity name "e" or "o".
    """
    if parity_value in BATCH_PROCESSING_DIRECTORY_NAME_PARITY:
        return BATCH_PROCESSING_DIRECTORY_NAME_PARITY[parity_value]
    else:
        raise ValueError("Unexpected parity value for a: " + parity_value)

def get_batch_processing_directory_name(real_k: str, imag_k: str, params: dict) -> str:
    """
    Returns the name of the batch-processing directory for the specified resonance k.
    """

    # prepare the parity code (".ee.", ".eo.", ".o.", and so on) to be embedded in the directory name
    parity_code = "."
    if CONFIG_PARAM_NAME_PARITY_A in params:
        parity_code += convert_parity_value_to_parity_name( params[CONFIG_PARAM_NAME_PARITY_A] )
    if CONFIG_PARAM_NAME_PARITY_B in params:
        parity_code += convert_parity_value_to_parity_name( params[CONFIG_PARAM_NAME_PARITY_B] )
    if parity_code != ".":
        parity_code += "."

    # if the value of imag_k(type: str) has sign:
    if imag_k[0] == "-" or imag_k[0] == "+":
        return BATCH_PROCESSING_DIRECTORY_NAME_PREFIX.format(parity_code) + real_k + imag_k + "i"

    # if the value of imag_k(type: str) has no sign:
    else:
        # it is necessary to insert "+" between real_k and imag_k to separate them in the directory name.
        return BATCH_PROCESSING_DIRECTORY_NAME_PREFIX.format(parity_code) + real_k + "+" + imag_k + "i"

def split_multi_param_def_line(multi_param_def_line: str) -> list:
    """
    Splits a multi-parameter-defining line into the list of single-parameter-defining-lines,  and returns it.
    """
    splitted_line_list = []
    last_split_point_index = 0
    is_in_quated_sector = False
    # loop for characters in the line
    for char_index in range(len(multi_param_def_line)):
        # split line by "," except it in quoted sectors, and stores splitted contents into the list
        if multi_param_def_line[char_index] == "\"":
            is_in_quated_sector = not is_in_quated_sector
        if multi_param_def_line[char_index] == "'":
            raise ValueError("Use double-quotation in config lines, not single-quotation: " + multi_param_def_line)
        if not is_in_quated_sector and multi_param_def_line[char_index] == ",":
            splitted_line_list.append( multi_param_def_line[last_split_point_index:char_index] )
            last_split_point_index = char_index + 1
    splitted_line_list.append( multi_param_def_line[last_split_point_index:-1] )
    return splitted_line_list


# ================================================================================
# - Main function / top-level flow of this script -
# ================================================================================

def main(argv: list) -> None:

    # parse arguments
    arg_parser = argparse.ArgumentParser(
        description="Tool to perform parallel-batch-processing of wfunc/husimi.sh for each resonance in the list."
    )
    arg_parser.add_argument(
        "--config", action="store", type=str, required=True, 
        help="specify the batch-configuration file"
    ) 
    arg_parser.add_argument(
        "--data", action="store", type=str, required=True, 
        help="specify the resonance-list data file"
    ) 
    arg_parser.add_argument(
        "--parallel", action="store", type=int, default=-1, 
        help="specify the number of processes to be executed in parallel [default: use all cpu cores]"
    ) 
    arg_parser.add_argument(
        "-v", "--version", action="version", version=APPLICATION_NAME + " " + APPLICATION_VERSION
    ) 
    args = arg_parser.parse_args()

    # get values of arguments (paths of input files, number of parallel processes, etc.)
    config_file_path = args.config
    data_file_path = args.data
    num_parallels = args.parallel
    max_available_cores = multiprocessing.cpu_count()
    if num_parallels == -1:
        num_parallels = max_available_cores
    
    print("\n" + APPLICATION_NAME + ": " + str(num_parallels) + "-cores are assigned for batch-processing" \
            + " (max: " + str(max_available_cores) + "-cores)")

    # load config/data files
    config = BatchConfiguration()
    config.loadFile(config_file_path)
    data = ResonanceData()
    data.loadFile(data_file_path)

    print("\n" + APPLICATION_NAME + ": generate batch resources...")

    # generate directories and shell-scripts for batch-processing
    batch_resource_generator = BatchResourceGenerator(config, data)
    batch_resource_generator.generate_directories()
    batch_resource_generator.generate_shell_scripts()

    print("\n" + APPLICATION_NAME + ": start batch processing...")

    # start batch-processing
    batch_process_controller = BatchProcessController(config, data, num_parallels)
    batch_process_controller.schedule()
    batch_process_controller.run()


# ================================================================================
# - Component classes -
# ================================================================================

class BatchConfiguration(object):
    """
    A class storing values loaded from a batch-configuration file (e.g. "config_ee.txt").
    """
    def __init__(self):
        self.parameter_dict = None
    
    def loadFile(self, config_file_path: str) -> None:
        """
        Loads configuration values from a batch-configuration file (e.g. "config_ee.txt").
        """
        self.parameter_dict = dict()

        # open the file and read all lines
        with open(config_file_path, "r") as file:
            for line in file:

                # comment lines
                if line[0] == "#":
                    pass

                # parameter lines (syntax: "name = value" or "name1 = value1, name2 = value2, ...")
                elif "=" in line:
                    for param_def in split_multi_param_def_line(line):  # multiple parameters can be defined in a line
                        name = param_def.split('=',1)[0].strip()  # extract left-hand-side of the first "="
                        value = param_def.split('=',1)[1].strip() # extract right-hand-side of the first "="
                        value = value.strip("\"").strip("'")      # unquote the value
                        self.parameter_dict[name] = value

        # resolve absolute paths of file-path-parameters, e.g.: wfunc_exe, husimi_exe, boundary_data
        self.resolve_absolute_paths_of_file_params()

    def resolve_absolute_paths_of_file_params(self) -> None:
        """
        Replaces the value of parameters specifying paths of executable programs with absolute paths.
        """
        for parameter_name in CONFIG_PARAM_NAMES_OF_FILES:
            if parameter_name in self.parameter_dict:
                self.parameter_dict[parameter_name] \
                    = os.path.join(os.getcwd(),  self.parameter_dict[parameter_name])

    def get_parameter(self, parameter_name: str) -> str:
        """
        Returns the value of the specified parameter, defined in the loaded file.
        """
        if self.parameter_dict is None:
            raise LookupError("Parameters of the Config object are not loaded yet.")

        if parameter_name in self.parameter_dict :
            return self.parameter_dict[parameter_name]
        else:
            raise LookupError("The parameter \"" + parameter_name + "\" is not defined in the loaded config file.")

    def get_parameter_dict(self) -> dict:
        """
        Returns the dict storing parameters defined in the loaded file.
        """
        if self.parameter_dict is None:
            raise LookupError("Parameters of the Config object are not loaded yet.")
        else:
            return self.parameter_dict

    def dump(self) -> None:
        """
        Dumps contents for debugging.
        """
        print("\n- Dump of Config object -\n")
        for key in self.parameter_dict.keys():
            print("parameter_dict[" + key + "] = " + self.parameter_dict[key])


class ResonanceData(object):
    """
    A class storing data loaded from a resonance-list data file (e.g. "data.resonances.ee").
    """
    def __init__(self):
        self.parameter_dict = None
        self.re_k_list = None
        self.im_k_list = None
        self.det_list = None

    def loadFile(self, config_file_path: str) -> None:
        """
        Loads data and parameters from a resonance-list data file (e.g. "data.resonances.ee").
        """
        self.parameter_dict = dict()
        self.re_k_list = []
        self.im_k_list = []
        self.det_list = []
        
        # open the file and read all lines
        with open(config_file_path, "r") as file:
            for line in file:

                # comment lines
                if (line[0] == "#"):
                    pass

                # blank lines
                elif not line.strip():
                    pass
                
                # data lines
                else:
                    columns = line.split()
                    if len(columns) == 3:
                        self.re_k_list.append(columns[0])
                        self.im_k_list.append(columns[1])
                        self.det_list.append(columns[2])
                    else:
                        raise ValueError(
                            "The number of columns in \"" + config_file_path 
                            + "\" should be 3 but detected " + str(len(columns)) + ". line: " + line
                        )

    def get_parameter(self, parameter_name: str) -> str:
        """
        Returns the value of the specified parameter, defined in a header line in the loaded file.
        """
        if self.parameter_dict is None:
            raise LookupError("Parameters of the Data object are not loaded yet.")

        if parameter_name in self.parameter_dict :
            return self.parameter_dict[parameter_name]
        else:
            raise LookupError("The parameter \"" + parameter_name + "\" is not defined in the loaded data file.")

    def get_parameter_dict(self) -> dict:
        """
        Returns the dict storing parameters defined in the loaded file.
        """
        if self.parameter_dict is None:
            raise LookupError("Parameters of the Config object are not loaded yet.")
        else:
            return self.parameter_dict

    def get_real_k_list(self) -> list:
        """
        Returns the list of values of re(k), defined as the 1st column in the loaded file.
        """
        if self.re_k_list is not None:
            return self.re_k_list
        else:
            raise LookupError("Data of re(k) could not be loaded expectedly from the data file.")

    def get_imag_k_list(self) -> list:
        """
        Returns the list of values of im(k), defined as the 2nd column in the loaded file.
        """
        if self.im_k_list is not None:
            return self.im_k_list
        else:
            raise LookupError("Data of im(k) could not be loaded expectedly from the data file.")

    def get_determinant_list(self) -> list:
        """
        Returns the list of determinat values, defined as the 3rd column in the loaded file.
        """
        if self.det_list is not None:
            return self.det_list
        else:
            raise LookupError("Data of determinant could not be loaded expectedly from the data file.")

    def dump(self) -> None:
        """
        Dumps contents for debugging.
        """
        print("\n- Dump of Data object -\n")
        for key in self.parameter_dict.keys():
            print("parameter_dict[" + key + "] = " + str(self.parameter_dict[key]))
        print("\ndata (re_k, im_t, det):")
        for i in range(len(self.re_k_list)):
            print(str(self.re_k_list[i]) + "\t" + str(self.im_k_list[i]) + "\t" + str(self.det_list[i]))


# Migrated from "JobParams" class in "autofinder.py", and modified
class FortranValueConverter(object):
    """
    A class for converting the format of floating-point numbers to it available in Fortran.
    """

    @staticmethod
    def to_fortran_float_string(variable: str) -> str:
        """
        Cast a Python variable to a 'Fortran' string.

        Parameter:
        ----------
        variable: str (should be interpretable as a floating-point number)
            variable to cast

        Returns:
        --------
            'Fortran' string

        Examples:
        ---------
            * FortranValueConverter.to_fortran_string(4.2) -> '4.2d0'
            * FortranValueConverter.to_fortran_string(42) -> '42'

        Throws:
        -------
            ValueError  
        """
        decimal, exponent = FortranValueConverter.filter_exponent_chars(variable)
        return f'{decimal}d{exponent}'

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
            * FortranValueConverter.filter_exponent_chars('4.2d-2') -> (4.2, -2) 

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


class BatchResourceGenerator(object):
    """
    A class for generating resources for batch-processing, e.g. directories and shell-scripts.
    """
    def __init__(self, batch_config: BatchConfiguration, resonance_data: ResonanceData):
        self.batch_config = batch_config
        self.resonance_data = resonance_data
        self.real_k_list = self.resonance_data.get_real_k_list()
        self.imag_k_list = self.resonance_data.get_imag_k_list()
    
    def generate_directories(self) -> None:
        """
        Creates batch-processing directories for all resonances.
        """
        # loop for each pair of elements stored in real_k_list and imag_k_list with the same index
        for real_k, imag_k in zip(self.real_k_list, self.imag_k_list):

            # create a directory for a resonance k (a pair of real_k and imag_k)
            dir_path = get_batch_processing_directory_name(real_k, imag_k, self.batch_config.get_parameter_dict())
            os.makedirs(dir_path, exist_ok=True)

    def generate_shell_scripts(self) -> None:
        """
        Generates wfunc.sh and husimi.sh in all batch-processing directories.
        """
        symmetry_code = self.batch_config.get_parameter(CONFIG_PARAM_NAME_SYM)
        resource_index = 0 # use for printing progress

        # get the absolute file path of plot scripts
        plot2d_file_path = os.path.join(os.path.dirname(__file__), PLOT2D_RELATIVE_PATH_FROM_BATCH_PLOT)
        nearfield_file_path = os.path.join(os.path.dirname(__file__), NEARFIELD_RELATIVE_PATH_FROM_BATCH_PLOT)
        farfield_file_path = os.path.join(os.path.dirname(__file__), FARFIELD_RELATIVE_PATH_FROM_BATCH_PLOT)

        # loop for each pair of elements stored in real_k_list and imag_k_list with the same index
        for real_k, imag_k in zip(self.real_k_list, self.imag_k_list):

            # create a path of a shell script to be generated
            dir_path = get_batch_processing_directory_name(real_k, imag_k, self.batch_config.get_parameter_dict())
            wfunc_script_path = os.path.join(dir_path, SHELL_SCRIPT_NAME_WFUNC)
            husimi_script_path = os.path.join(dir_path, SHELL_SCRIPT_NAME_HUSIMI)

            # generate wfunc.sh
            with open(wfunc_script_path, mode="w") as wfunc_script_file:
                exec_file_path = self.batch_config.get_parameter(CONFIG_PARAM_NAME_WFUNC_EXE_FILE)
                plot_option = self.batch_config.get_parameter(CONFIG_PARAM_NAME_WFUNC_PLOT_OPTION)
                logplot_option = self.batch_config.get_parameter(CONFIG_PARAM_NAME_WFUNC_LOGPLOT_OPTION)
                exec_command = get_exec_command_in_wfunc_shell_script(symmetry_code)
                plot_command = get_plot_command_in_wfunc_shell_script(
                    plot2d_file_path, plot_option, logplot_option, nearfield_file_path, farfield_file_path
                )
                self.write_contents_in_shell_script(
                    wfunc_script_file, real_k, imag_k, exec_file_path, exec_command, plot_command
                )
            
            # generate husimi.sh
            with open(husimi_script_path, mode="w") as husimi_script_file:
                exec_file_path = self.batch_config.get_parameter(CONFIG_PARAM_NAME_HUSIMI_EXE_FILE)
                plot_option = self.batch_config.get_parameter(CONFIG_PARAM_NAME_HUSIMI_PLOT_OPTION)
                exec_command = get_exec_command_in_husimi_shell_script(symmetry_code)
                plot_command = get_plot_command_in_husimi_shell_script(plot2d_file_path, plot_option)
                self.write_contents_in_shell_script(
                    husimi_script_file, real_k, imag_k, exec_file_path, exec_command, plot_command
                )

            # add execute-parmissions to generated shell scripts
            subprocess.call("chmod +x " + shlex.quote(wfunc_script_path), shell=True, cwd=os.getcwd())
            subprocess.call("chmod +x " + shlex.quote(husimi_script_path), shell=True, cwd=os.getcwd())

            # print progress
            print("[" + str(resource_index+1) + "] " + wfunc_script_path)
            print("[" + str(resource_index+2) + "] " + husimi_script_path)
            resource_index += 2

    def write_contents_in_shell_script(self, shell_script_file: str, real_k: str, imag_k: str, \
            exec_file_path: str, exec_command: str, plot_command: str) -> None:
        """
        Writes contents of wfunc.sh and husimi.sh.
        """
        # write the shebang (e.g. "#!/bin/bash -f" )
        shell_script_file.write(SHELL_SCRIPT_SHEBANG + "\n")

        # set current directory
        shell_script_file.write("cd $(dirname $0)\n")

        # set to exit immediately when any error occurred
        shell_script_file.write("set -e\n")

        # write the definition line of the executable file
        shell_script_file.write(SHELL_SCRIPT_PARAM_NAME_EXECFILE + "=" + shlex.quote(exec_file_path) + "\n")

        # write all parameters defined in the batch-configuration file
        parameter_dict = self.batch_config.get_parameter_dict()
        for parameter_name in parameter_dict:
            if parameter_name not in SHELL_SCRIPT_PARAM_NAMES_EXCLUDED:
                shell_script_file.write(parameter_name + "=" + shlex.quote(parameter_dict[parameter_name]) + "\n")
        
        # convert the format of resonance parameters to it of declaration of floating-point numbers in fortran77
        try:
            real_k_converted = FortranValueConverter.to_fortran_float_string(real_k)
            imag_k_converted = FortranValueConverter.to_fortran_float_string(imag_k)
        except ValueError:
            raise ValueError(f"Uninterpretable floating-point format are used for resonance data: ({real_k}, {imag_k})")

        # write resonance parameters
        shell_script_file.write(SHELL_SCRIPT_PARAM_NAME_REAL_K + "=" + shlex.quote(real_k_converted) + "\n")
        shell_script_file.write(SHELL_SCRIPT_PARAM_NAME_IMAG_K + "=" + shlex.quote(imag_k_converted) + "\n")

        # write execution commands
        shell_script_file.write(exec_command + "\n")
        shell_script_file.write(plot_command + "\n")


class BatchProcessController(object):
    """
    A class to controlling processes of parallel-batch-processing.
    """
    def __init__(self, batch_config: BatchConfiguration, resonance_data: ResonanceData, num_parallel_processes: int):
        self.batch_config = batch_config
        self.resonance_data = resonance_data
        self.process_command_list = None  #type: list of dict, and each element has values with keys: "priority" and "command".
        self.num_total_processes = -1
        self.num_available_cores = num_parallel_processes

    def schedule(self) -> None:
        """
        Prepares for batch-processing.
        """
        self.process_command_list = []

        # loop for each pair of elements stored in real_k_list and imag_k_list with the same index
        for real_k, imag_k in zip(self.resonance_data.get_real_k_list(), self.resonance_data.get_imag_k_list()):
            directory_name = get_batch_processing_directory_name(real_k, imag_k, self.batch_config.get_parameter_dict())
            wfunc_path = os.path.join("./", directory_name, SHELL_SCRIPT_NAME_WFUNC)
            husimi_path = os.path.join("./", directory_name, SHELL_SCRIPT_NAME_HUSIMI)

            # register commands to execute as parallel processes
            self.process_command_list.append( { "priority":SHELL_SCRIPT_PRIORITY_WFUNC,  "command": shlex.quote(wfunc_path)  } )
            self.process_command_list.append( { "priority":SHELL_SCRIPT_PRIORITY_HUSIMI, "command": shlex.quote(husimi_path) } )

        self.num_total_processes = len(self.process_command_list)

        # sort elements in process_command_list by values of priorities
        self.process_command_list.sort(key=operator.itemgetter("priority"))


    def run(self) -> None:
        """
        Runs all processes in parallel.
        """
        if self.process_command_list is None or len(self.process_command_list) == 0:
            raise LookupError(
                "No commands to be processed are registered. Note that schedule() should be called before calling run()."
            )

        running_process_list = [ None ] * self.num_available_cores  # index: [core_index]
        process_index = 0
        core_index = 0

        # loop for all processes(commands)
        while True:

            # dispatch processes to all cores, and execute them
            for core_index in range(self.num_available_cores):
                if process_index < self.num_total_processes:

                    # execute a process on a core
                    command = self.process_command_list[process_index]["command"]
                    running_process_list[core_index] = subprocess.Popen(command, shell=True)

                    # print progress information
                    print(
                        "[{}/{}]".format(process_index+1,self.num_total_processes) \
                        + "(core-{})".format(core_index+1) \
                        + ": " + command
                    )
                else:
                    running_process_list[core_index] = None
                process_index += 1
            
            # synchronize processes of all cores
            for core_index in range(self.num_available_cores):
                if running_process_list[core_index] is not None:
                    running_process_list[core_index].wait()

            # error detections of processes: cancel all remained processes if any errors occurred
            for core_index in range(self.num_available_cores):
                if running_process_list[core_index] is not None and running_process_list[core_index].returncode != 0:
                    print("\n\n" + APPLICATION_NAME + ": Some errors occurred in batch-processes. Exit...\n")
                    process_index = self.num_total_processes
                    break

            # if there is no more processes to be executed, break from the while-loop.
            if self.num_total_processes <= process_index:
                break


# ================================================================================
# - Entry point -
# ================================================================================

try:
    if __name__ == "__main__":
        main(sys.argv)
        sys.exit(0)

except KeyboardInterrupt:
    sys.stderr.write("\n" + APPLICATION_NAME + ": Exit...\n")
    sys.exit(1)
