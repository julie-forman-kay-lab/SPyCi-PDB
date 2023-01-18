"""
SPyCi-PDB Exceptions.

Inspired from:
https://github.com/julie-forman-kay-lab/IDPConformerGenerator/blob/3aef6b085ec09eeebc5812639a5eb6832c0215cd/src/idpconfgen/core/exceptions.py
"""
from spycipdb import count_string_formatters, log


class SPyCiPDBException(Exception):
    r"""
    SPyCi-PDB base exception.

    Parameters
    ----------
    *args
        If the first element of args contains ``'{}'`` it will be
        used as :attr:`errmsg` base string.
        Else, ``args`` are used to feed the sting method ``.format()``
        for the default exception :attr:`errmsg`.

    errmsg : optional
        If given, overrides any previous parameter and the ``str``
        value of ``errmsg`` is used as the Exception message.
        Defaults to ``None``.

    Examples
    --------
    Uses the default errormsg.
    >>> err = SPyCiPDBException(var1, var2)

    >>> err = SPyCiPDBException('An error happened: {}, {}', var1, var2)

    >>> err = SPyCiPDBException('An error happened')

    >>> err = SPyCiPDBException(errmsg='Custom error msg')

    >>> err = SPyCiPDBException(errmsg='')

    """

    errmsg = 'An unknnown error as occurred.'

    def __init__(self, *args, errmsg=None):

        # SPyCiPDBException(errmsg='Custom error msg')
        if errmsg is not None:
            assert isinstance(errmsg, str), f'wrong errmsg type: {type(errmsg)}'
            self.errmsg = errmsg
            self.args = []

        elif len(args) == count_string_formatters(self.errmsg):
            self.args = args

        else:
            assert count_string_formatters(args[0]) == len(args[1:]), \
                'args passed to Exception are not compatible to form a message'
            self.errmsg = args[0]
            self.args = args[1:]

        log.debug(f'Exception errors: {self.errmsg}')
        log.debug(f'Exception args: {self.args}')

        # ensure
        assert isinstance(self.args, (tuple, list)), \
            f'wrong args {type(self.args)}'
        assert count_string_formatters(self.errmsg) == len(self.args), (
            'Bad Exception message:\n'
            f'errmsg: {self.errmsg}\n'
            f'args: {self.args}'
            )

    def __str__(self):
        """Make me a string."""
        return self.errmsg.format(*self.args)

    def __repr__(self):
        return f'{self.__class__.__name__}: {self}'

    def report(self):
        """
        Report error in the form of a string.

        Identifies Error type and error message.

        Returns
        -------
        str
            The formatted string report.
        """
        return f'{self.__class__.__name__} * {self}'


class ReportOnCrashError(SPyCiPDBException):
    """Raised when logger.report_on_crash."""

    errmsg = "Crash reported to {}."""
