from warnings import warn
from lcapy.parser import Parser
from lcapy import grammar
from lcapy import mnacpts
from typing import Union


class NetlistLine:
    def __init__(self, line: str, validate: bool = True):
        self.line = line.replace('{', '').replace('}', '')

        # parse line
        elementParam, drawParam = self.line.split(';')
        self.drawParam = drawParam.replace(' ', '')

        values = elementParam.split(' ')
        values = [value for value in values if not value == '']

        for i in range(len(values)):
            values[i] = values[i].replace('{', '').replace('}', '')

        if len(values) < 3:
            raise RuntimeError("Cant parse netlist line: %s"
                               "make sure each line has a name startNode endNode; drawing annotation "
                               "e.g. looks like: "
                               "V1 0 1 dc {10}; up --- "
                               "W 2 3; left --- "
                               "C3 6 7 {100}; down", self.line)

        self.parser: Parser = NetlistLineParser()

        self.cpt = self.parser.parse(line)
        self.type = self.cpt.type
        self.typeSuffix = self.cpt.name.replace(self.type, "", 1)
        try:
            self.startNode = int(self.cpt.node_names[0])
            self.endNode = int(self.cpt.node_names[1])
        except ValueError:
            # for sorting in drawWithSchemdraw, Nodes need to be integers
            raise ValueError(f"can't convert {values[1]} or {values[2]} to int. start- and endNode have to be integers")

        if self.type in ["R", "L", "C", "Z"]:
            self.ac_dc, self.value, self.phase, self.omega = self.parseRLCZ()
        elif self.type == "W":
            self.ac_dc, self.value, self.phase, self.omega, self.typeSuffix = self.parseW()
        elif self.type == "V":
            self.ac_dc, self.value, self.phase, self.omega = self.parseV()

    def is_AC_or_DC(self) -> Union[str, None]:
        if not self.type == "V":
            return None
        if self.cpt.has_ac and not self.cpt.has_dc:
            return "ac"
        elif self.cpt.has_dc and not self.cpt.has_ac:
            return "dc"
        else:
            raise RuntimeError(f" if type is V (is: {self.type}) should be ac or dc not both")

    @staticmethod
    def parseW() -> (None, None, None, None, str):
        # self.typeSuffix for wires has to be "", that in drawWithSchemdraw the id-tag in the svg-file is empty
        return None, None, None, None, ""

    def parseRLCZ(self) -> (None, any, None):
        ac_dc = None
        value = self.cpt.args[0]
        phase = None
        omega = None

        return ac_dc, value, phase, omega

    def parseV(self) -> (str, str, str, str):
        value = self.cpt.args[0]
        if self.cpt.has_ac:
            if len(self.cpt.args) > 1:
                phase = self.cpt.args[1]
            else:
                phase = None
            if len(self.cpt.args) > 2:
                omega = self.cpt.args[2]
            else:
                omega = None
        else:
            phase = None
            omega = None

        ac_dc = self.is_AC_or_DC()

        return ac_dc, value, phase, omega

    def label(self):
        if not self.type == "W":
            return self.type + self.typeSuffix
        else:
            return ""

    def reconstruct(self) -> str:
        """
        reconstructs self.line from the parsed elements self.type, self.typeSuffix, self.startNode, self.endNode,
        self.ac_dc, self.value, self.phase, self.omega, self.drawParam
        :return: reconstructed string
        """
        if isinstance(self.value, str):
            self.value = self.value.replace("{", "").replace("}", "")

        if self.type == "W":
            reconstructed =\
                f"{self.type+self.typeSuffix} {self.startNode} {self.endNode}; {self.drawParam}"
        elif self.type in ["R", "L", "C", "Z", "ZR", "ZL", "ZC"]:
            reconstructed =\
                f"{self.type+self.typeSuffix} {self.startNode} {self.endNode} {{{self.value}}}; {self.drawParam}"
        elif self.type == "V" and not self.omega and not self.phase:
            reconstructed =\
                (f"{self.type+self.typeSuffix} {self.startNode} {self.endNode} {self.ac_dc} {self.value};"
                 f" {self.drawParam}")
        elif self.type == "V" and not self.omega and self.phase:
            reconstructed =\
                (f"{self.type+self.typeSuffix} {self.startNode} {self.endNode} {self.ac_dc} {self.value} {self.phase};"
                 f" {self.drawParam}")
        else:
            reconstructed =\
                (f"{self.type+self.typeSuffix} {self.startNode} {self.endNode} {self.ac_dc} {self.value} {self.phase} "
                 f"{self.omega}; {self.drawParam}")

        return reconstructed

    def validate_parsing(self):
        """
        raises an error if parsing fails and warns if it may fail. May fail if the parsed line cant be reconstructed but
        the reconstructed and parsed line without white spaces match.
        :return: void
        """
        # check if the parsing was successful if the line can be reconstructed it should be parsed correctly
        # white space sensitive but without "{"; "}"
        ref = self.line
        # white space insensitive
        ref2 = ref.replace(' ', '')

        if not self.reconstructed == ref and not self.reconstructed == ref2:
            raise RuntimeError(f"Error while parsing {self.line}: reconstructed -> {self.reconstructed}")
        if not ref == self.reconstructed and ref2 == self.reconstructed:
            warn(f"potential error while parsing {self.line}: reconstructed -> {self.reconstructed}")

    def __str__(self):
        return str(self.cpt)


class NetlistLineParser:
    _instance = None

    def __new__(cls, *args, **kwargs) -> Parser:
        if cls._instance is None:
            cls._instance = Parser(mnacpts, grammar, False)
        return cls._instance

    def __init__(self, value=None):
        if not hasattr(self, '_initialized'):  # Ensure `__init__` runs only once
            self.value = value
            self._initialized = True
