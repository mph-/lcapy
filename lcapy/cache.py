"""This module provides the cached_property decorator.

Copyright 2021 Michael Hayes, UCECE
"""

try:
    # Requires python3.8
    from functools import cached_property
except:
    from property_cached import cached_property

from functools import lru_cache
