from .rcdefaults import rcdefaults
from .config import rcparams_user_filename, rcparams_local_filename
from os.path import exists, expanduser
from warnings import warn


class RcParams(dict):

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):

        # Handle aliases here...

        if key not in rcdefaults:
            raise KeyError(f'Unknown param {key}')

        d = rcdefaults[key]

        # Check the value
        if isinstance(d[1], tuple):
            options = d[1]
            if val not in options:
                s = ', '.join(options)
                return ValueError(f'Invalid value {val} for {key}, options: {s}')
        else:
            val = d[1](val)

        dict.__setitem__(self, key, val)

        # Change hook
        if len(d) > 2:
            d[2](val)

    def __getitem__(self, key):

        # Handle aliases here...

        if key not in rcdefaults:
            raise KeyError(f'Unknown param {key}')

        return dict.__getitem__(self, key)

    def load(self, filename):

        rc = {}
        with open(filename) as fd:
            for linenum, line in enumerate(fd, 1):
                sline = line.split('#', 1)[0].strip()
                if not sline:
                    continue
                parts = sline.split(':', 1)
                if len(parts) != 2:
                    warn('Misformed key-val pair in file %r line #%d.' % (filename, linenum))
                    continue
                key, val = parts
                key = key.strip()
                val = val.strip()
                if key in rc:
                    warn('Duplicate key %s in file %r line #%d.' %
                         (key, filename, linenum))
                rc[key] = val
                self[key] = val

    def load_if_exists(self, filename):

        if not exists(filename):
            return

        self.load(filename)

    def load_user(self):

        filename = expanduser(rcparams_user_filename)
        self.load_if_exists(filename)

    def load_local(self):

        self.load_if_exists(rcparams_local_filename)

    def save(self, filename):

        with open(filename, 'w') as fd:

            for key, val in self.items():
                print(key, ':', val, file=fd)

    def set_defaults(self):

        for k, v in rcdefaults.items():
            rcParams[k] = v[0]


rcParams = RcParams()
rcParams.set_defaults()
