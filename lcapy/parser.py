from plyplus import Grammar, ParseError
from grammar import grammar

parser = Grammar(grammar)

class Parser(object):

    def __init__(self, cpts):

        self.cpts = cpts

    
    def parse(self, string):

        try:
            thing = parser.parse(string)
        except ParseError as e:
            raise ParseError('Could not parse: %s: due to %s' % (string, e))
        
        classname = str(thing.tail[0].head).capitalize()
        name = str(thing.tail[0].tail[0].tail[0])
        
        try:
            newclass = getattr(self.cpts, classname)
        except:
            newclass = self.cpts.newclasses[classname]

        obj = newclass(name)

        def frepr(self):
            return self.string
        
        obj.__repr__ = frepr
        obj.thing = thing

        # Add attributes.
        obj.nodes = ()
        for field in thing.tail[0].tail[1:]:
            attr, val = field.head, field.tail

            # There must be a cleaner way to do this!
            if attr == 'opts':
                opts = {}
                for val1 in val:
                    if val1.head == 'keyword':
                        opts[val1.tail[0]] = True
                    elif val1.head == 'keypair':
                        opts[val1.tail[0].tail[0]] = val1.tail[1].tail[0][1:]
                    else:
                        raise ValueError('Unknown field %s' % val1.head)
                val = opts

            else:
                try:
                    val = val[0]
                except:
                    pass

            if 'node' in attr:
                obj.nodes += (val, )
            setattr(obj, attr, val)

        obj.string = string

        return obj
