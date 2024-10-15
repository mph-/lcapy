from lcapy import state
from lcapy import resistance, inductance, capacitance, voltage, impedance, current
from lcapy.units import ohms, farads, henrys, amperes, volts
from lcapy.mnacpts import R, C, L, Z, V, I
from typing import Union
from lcapy import ConstantFrequencyResponseDomainExpression, ConstantFrequencyResponseDomainImpedance
from lcapy import ConstantDomainExpression


class UnitWorkAround:
    @staticmethod
    def addUnit(val, cptType):
        state.show_units = True
        if cptType == "R":
            returnVal = resistance(val)
        elif cptType == "C":
            returnVal = capacitance(val)
        elif cptType == "L":
            returnVal = inductance(val)
        elif cptType == "Z":
            returnVal = impedance(val)
        elif cptType == "V":
            returnVal = voltage(val)
        elif cptType == "I":
            returnVal = current(val)
        elif cptType == "W":
            return val
        else:
            raise NotImplementedError(f"{cptType} not supported edit UnitWorkAround.addUnit to support")
        return returnVal

    @staticmethod
    def getUnit(element: Union[R, C, L, Z]) -> (
            ConstantFrequencyResponseDomainExpression or ConstantFrequencyResponseDomainImpedance):
        """
        returns the unit of an element
        for R 1*ohm
        for C 1*F
        for L 1*H
        for Z 1*ohm (impedance has unit ohm)
        :param element: element: mnacpts.R | mnacpts.C | mnacpts.L | mnacpts.Z
        :return: for R, C, L ConstantFrequencyResponseDomainExpression; for Z ConstantFrequencyResponseDomainImpedance
        """
        if isinstance(element, R):
            return ohms
        elif isinstance(element, C):
            return farads
        elif isinstance(element, L):
            return henrys
        elif isinstance(element, Z):
            return ohms
        elif isinstance(element, V):
            return amperes
        elif isinstance(element, I):
            return volts
        else:
            raise NotImplementedError(f"{type(element)} not supported edit UnitWorkAround.addUnit to support")

    @staticmethod
    def getElementSpecificValue(element: Union[R, C, L, Z], unit=False) -> ConstantDomainExpression:
        """
        accesses the value resistance, capacitance, inductance, or impedance of an element based on its type
        :param element: mnacpts.R | mnacpts.C | mnacpts.L | mnacpts.Z
        :param unit: if True the Unit (ohm, F, H) are added to the str
        :return: lcapy.ConstantDomainExpression
        """
        if unit:
            return UnitWorkAround.addUnit(UnitWorkAround.getElementSpecificValue(element), element.type)

        state.show_units = False
        if isinstance(element, R):
            returnVal = element.R
        elif isinstance(element, C):
            returnVal = element.C
        elif isinstance(element, L):
            returnVal = element.L
        elif isinstance(element, Z):
            returnVal = element.Z
        elif isinstance(element, V):
            returnVal = element.V
        elif isinstance(element, I):
            returnVal = element.I
        else:
            raise NotImplementedError(f"{type(element)} "
                                      f"not supported edit UnitWorkAround.getElementSpecificValue() to support")

        return returnVal
