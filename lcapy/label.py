class Label:

    def __init__(self, key, val):

        # Circuitikz does not handle \\ in labels.  As a hack
        # convert l to l2 and a to a2.   Note, l2 and a2 can only
        # support two lines.
        parts = val.split(r'\\', 1)
        if len(parts) == 2:
            if key[0] == 'l':
                key = key.replace('l', 'l2')
                val = ' and '.join(parts)
            elif key[0] == 'a':
                key = key.replace('a', 'a2')
                val = ' and '.join(parts)

        self.key = key
        self.val = val

    def __str__(self):

        return self.key + '=' + self.val
