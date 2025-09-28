class RcChecker:

    def bool(self, val):

        if val in ('t', 'T', 'true', 'True', True):
            return True
        elif val in ('f', 'F', 'false', 'False', False):
            return False
        raise ValueError(f'Invalid bool: {val}')

    def int(self, val):

        try:
            return int(val)
        except ValueError:
            raise ValueError(f'Invalid int: {val}')

    def float(self, val):

        try:
            return float(val)
        except ValueError:
            raise ValueError(f'Invalid float: {val}')

    def str(self, val):
        return val

    def choice(self, val, options):

        if val not in options:
            s = ', '.join(options)
            return ValueError(f'Invalid value {val}, options: {s}')
        return val
