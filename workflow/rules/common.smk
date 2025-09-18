#!/bin/env python

import os
import warnings

################### Helper Functions ###########################################

def get_sample_files(wildcards):
    """
    Example helper function to get sample files.
    Modify according to your needs.
    """
    return config["samples"][wildcards.sample]

def get_all_samples():
    """
    Get all sample names from config.
    """
    return list(config["samples"].keys())

def is_odd(file):
    """
    Example helper function to check if a sample length is odd.
    Modify according to your needs.
    """

    with open(file, "r") as f:
        length = len(f.read())
    
    if length % 2 == 1:
        return "odd"
    else:
        return "even"

# Add your common functions here
