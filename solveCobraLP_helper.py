#############################################################################################################
# Python minumal implementation of helper function to the solveCobraLP function, from Matlab COBRA Toolbox
# Daniel Ribeiro 2021
#############################################################################################################

import numpy as np
from libsmop import *

@function
def columnVector(vec, *args, **kwargs):
    # Daniel Ribeiro 2021: Python implementation of the function columnVector from Matlab COBRA Toolbox
    # Converts a vector to a column vector
    #
    # USAGE:
    #
    #   vecT = columnVector(vec)
    #
    # INPUT:
    #   vec:     a vector
    #
    # OUTPUT:
    #   vecT:    a column vector
    #
    # .. Authors:
    #     - Original file: Markus Herrgard
    #     - Minor changes: Laurent Heirendt January 2017
    varargin = columnVector.varargin
    nargin = columnVector.nargin

    n, m = vec.shape

    if (m != 1 and n < m) or n == 1:
        vecT = vec.T
    else:
        vecT = vec

    return vecT

@function
def updateStructData(origStruct, updateStruct, *args, **kwargs):
    # Daniel Ribeiro 2021: Python implementation of the function updateStructData from Matlab COBRA Toolbox
    # Update the struct in origStruct with the data from updateStruct
    #
    # USAGE:
    #    updatedStruct = updateStruct(origStruct,updateStruct)
    #
    # INPUTS:
    #    origStruct:        The original Struct
    #    updateStruct:      The struct to update the information in origStruct.
    #
    # OUTPUT:
    #    updatedStruct:     The struct with the information from origStruct
    #                       updated by the info from updateStruct
    varargin = updateStructData.varargin
    nargin = updateStructData.nargin
    
    updateFieldNames = updateStruct.__dict__.keys(); # get all field names from the original struct that contains the additional information

    # set up the updated struct
    updatedStruct = origStruct

    for field in updateFieldNames:
        # for each field, if it is not a struct, replace the information. If it
        # is a struct in the updated struct, replace the original value,
        # otherwise update the subvalues.
        if  not hasattr(updatedStruct, field) or (type(getattr(updateStruct, field)) != struct):
            setattr(updatedStruct, field, getattr(updateStruct, field))
        else:
            if type(getattr(origStruct, field)) == struct:
                setattr(updatedStruct, field, updateStructData(getattr(origStruct, field), getattr(updateStruct, field)))
    
    return updatedStruct