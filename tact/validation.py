# -*- coding: utf-8 -*-

"""Various validation functions for inputs."""


def validate_outgroups(ctx, param, value):
    if value is None:
        return
    try:
        value = value.split(",")
    except AttributeError:
        # Tuples and lists shouldn't have the .split method
        pass
    return [x.replace("_", " ") for x in value]
