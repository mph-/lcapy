from plyplus import Grammar, ParseError
from grammar import grammar
import schemcpts as cpts

parser = Grammar(grammar)


def parse(string):
    
    try:
        thing = parser.parse(string)
    except ParseError as e:
        raise ParseError('Could not parse: %s: due to %s' % (string, e))
        
    classname = str(thing.tail[0].head).capitalize()
    name = str(thing.tail[0].tail[0].tail[0])

    try:
        newclass = getattr(cpts, classname)
    except:
        newclass = cpts.newclasses[classname]

    obj = newclass(name)

    def frepr(self):
        return self.string
        
    obj.__repr__ = frepr

    # Add attributes.
    obj.nodes = ()
    for field in thing.tail[0].tail[1:]:
        attr, val = field.head, field.tail[0] 
        if 'node' in attr:
            obj.nodes += (val, )
        setattr(obj, attr, val)

    obj.string = string


    return obj
