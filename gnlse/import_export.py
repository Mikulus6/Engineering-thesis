"""Import and export \\*.json files.

This module contains functions that enable to read json files (\\*.json)
in python as dictionary, and to export dictionary to \\*.json.

"""

import os
import json
import numpy as np

def read_json(filename):

    with open(filename) as file:
        json_dict = json.loads(file.read())
        assert isinstance(json_dict, dict)
        for key in json_dict.keys():
            json_dict[key] = np.array(unserialize_complex(json_dict[key]))
        return json_dict

def write_json(data: dict, filename):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    json_ready = {key: serialize_complex(value.tolist()) for key, value in data.items()}
    with open(filename, 'w') as file:
        json.dump(json_ready, file)

def serialize_complex(nested: list) -> list:
    for index_value in range(len(nested)):
        match type(nested[index_value]).__name__:
            case "list":  nested[index_value] = serialize_complex(nested[index_value])
            case "complex": nested[index_value] = [nested[index_value].real, nested[index_value].imag]
            case _: pass
    return nested

def unserialize_complex(nested: list) -> list:
    for index_value in range(len(nested)):
        if type(nested[index_value]) == list:
            if len(nested[index_value]) == 2: nested[index_value] = complex(*nested[index_value])
            else:                             nested[index_value] = unserialize_complex(nested[index_value])
    return nested
