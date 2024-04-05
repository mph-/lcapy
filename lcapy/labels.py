from .label import Label
from .latex import latex_format_label


class Labels:

    def __init__(self):

        self.current = None
        self.voltage = None
        self.flow = None
        self.inner = None
        self.annotation = None
        self.label = None

    def add(self, key, val):

        val = latex_format_label(val)

        if key[0] == 'a':
            self.annotation = Label(key, val)
        elif key[0] == 'f':
            self.flow = Label(key, val)
        elif key[0] == 'i':
            if key == 'ir':
                key = 'i' + '<'
            self.current = Label(key, val)
        elif key[0] == 'l':
            self.label = Label(key, val)
        elif key[0] == 't':
            self.inner = Label(key, val)
        elif key[0] == 'v':
            if key == 'vr':
                key = key + '>'
            self.voltage = Label(key, val)

    def args(self, flip=False):

        def fmt2(label, flip=False):
            key, val = label.key, label.val

            if '^' not in key and '_' not in key:
                key = key + ('^' if flip else '_')

            return key + '=' + val

        def fmt(label):
            key, val = label.key, label.val
            return key + '=' + val

        args = []
        if self.annotation is not None:
            args.append(fmt2(self.annotation, not flip))
        if self.voltage is not None:
            args.append(fmt2(self.voltage, not flip))
        if self.flow is not None:
            args.append(fmt(self.flow))
        if self.current is not None:
            args.append(fmt(self.current))
        if self.inner is not None:
            args.append(fmt(self.inner))
        if self.label is not None:
            args.append(fmt2(self.label, flip))
        return args

    def __str__(self):

        args = self.args()
        return ', '.join(args)
