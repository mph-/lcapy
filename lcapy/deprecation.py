"""
Deprecation warning.

This is adapted from sympy/utilities/exceptions.py
"""

import warnings

from sympy.utilities.misc import filldedent


class LcapyDeprecationWarning(DeprecationWarning):
    r"""A warning for deprecated features of Lcapy.
    """

    def __init__(self, feature=None, last_supported_version=None,
                 useinstead=None, issue=None, deprecated_since_version=None):

        self.fullMessage = ""

        if not feature:
            raise ValueError(
                "feature is required argument of LcapyDeprecationWarning")

        if not deprecated_since_version:
            raise ValueError(
                "deprecated_since_version is a required argument of LcapyDeprecationWarning")

        self.fullMessage = "%s has been deprecated since Lcapy %s. " % \
            (feature, deprecated_since_version)

        if last_supported_version:
            self.fullMessage += ("It will be last supported in Lcapy "
                                 "version %s. ") % last_supported_version

        if useinstead:
            self.fullMessage += "Use %s instead. " % useinstead

        if issue is not None:
            self.fullMessage += ("See "
                                 "https://github.com/mph-/lcapy/issues/%d for more "
                                 "info. ") % issue

    def __str__(self):
        return '\n%s\n' % filldedent(self.fullMessage)

    def warn(self, stacklevel=2):
        # the next line is what the user would see after the error is printed
        # if stacklevel was set to 1. If you are writing a wrapper around this,
        # increase the stacklevel accordingly.
        warnings.warn(self, stacklevel=stacklevel)


# Python by default hides DeprecationWarnings, which we do not want.
warnings.simplefilter("once", LcapyDeprecationWarning)
