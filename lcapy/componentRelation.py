from enum import Enum


class ComponentRelation(Enum):
    parallel = 'parallel'
    series = 'series'
    none = None

    def __eq__(self, other):
        if isinstance(other, str):
            return self.value == other
        elif isinstance(other, ComponentRelation):
            return self.value == other.value
        else:
            raise TypeError

    def to_string(self):
        return self.__str__()

    def __str__(self):
        return str(self.value)
