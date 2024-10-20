def value_parser(arg):
    """Convert args such as 42p to 42e-12.

    The suffixes are ignored in expressions such as 10 * 42p."""

    suffixes = {'f': 1e-15, 'p': 1e-12, 'n': 1e-9, 'u': 1e-6,
                'm': 1e-3, 'k': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}

    def isfloat(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    if arg is None or len(arg) <  2:
        return arg

    if arg.endswith('Meg'):
        arg = arg[0:-2] + 'M'
    elif arg.endswith('K'):
        arg = arg[0:-1] + 'k'

    if (arg[-1] in suffixes and isfloat(arg[0:-1])):
        return float(arg[0:-1]) * suffixes[arg[-1]]
    else:
        return arg
