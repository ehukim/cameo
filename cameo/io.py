# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import types
import pickle
from cobra.io import read_sbml_model
from cameo.solver_based_model import OptlangBasedModel, to_solver_based_model

def load_model(path_or_handle):
    """Read a model from a file.

    Arguments
    ---------

    path_or_handle: either a file path or a file handle

    Supported filetypes
    -------------------

    - pickle
    - sbml

    """

    if isinstance(path_or_handle, types.StringType):
        path = path_or_handle
        handle = open(path_or_handle)
    elif hasattr(path_or_handle, 'read'):
        path = path_or_handle.name
        handle = path_or_handle
    else:
        raise ValueError('Provided argument %s has to be either a file path or handle' % path_or_handle)
    logger.debug('Reading file from %s assuming pickled model.' % path)
    try:
        model = pickle.load(handle)
    except Exception:
        logger.debug('Cannot unpickle %s. Assuming sbml model next.' % path)
        try:
            model = read_sbml_model(path)
        except AttributeError:  # TODO: cobrapy doesn't raise a proper exception if a file does not contain an SBML model
            raise ValueError('FIXME')
    
    if not isinstance(model, OptlangBasedModel):
        return to_solver_based_model(model)
    else:
        return model
            