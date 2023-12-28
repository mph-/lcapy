from .label import Label


class Labels:

    def __init__(self, flip=False):

        self.flip = flip

        self.current = None
        self.voltage = None
        self.flow = None
        self.inner = None
        self.annotation = None
        self.label = None

    def add(self, key, val):

        if key in ('a', 'v'):
            key = key + ('_' if self.flip else '^')
        elif key in ('l', 'i'):
            key = key + ('^' if self.flip else '_')

        if key[0] == 'a':
            self.annotation = Label(key, val)
        elif key[0] == 'f':
            self.flow = Label(key, val)
        elif key[0] == 'i':
            if key == 'ir':
                key = key + '<'
            self.current = Label(key, val)
        elif key[0] == 'l':
            self.label = Label(key, val)
        elif key[0] == 't':
            self.inner = Label(key, val)
        elif key[0] == 'v':
            if key == 'vr':
                key = key + '>'
            self.voltage = Label(key, val)

    def args(self):

        args = []
        if self.annotation is not None:
            args.append(str(self.annotation))
        if self.flow is not None:
            args.append(str(self.flow))
        if self.current is not None:
            args.append(str(self.current))
        if self.inner is not None:
            args.append(str(self.inner))
        if self.voltage is not None:
            args.append(str(self.voltage))
        if self.label is not None:
            args.append(str(self.label))
        return args

    def __str__(self):

        args = self.args()
        return ', '.join(args)
