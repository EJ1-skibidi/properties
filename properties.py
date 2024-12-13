import clr
from System import String, Double, Array, Reflection
dtlpath = "\\ElliotF\\Desktop\\Article\\DWSIM Thermo\\DTL_v8.3.1\\"
clr.AddReference(dtlpath + "DWSIM.Thermodynamics.StandaloneLibrary.dll")
from DWSIM.Thermodynamics import PropertyPackages, CalculatorInterface
import CapeOpen

class Fluid:
    def __init__(self, components, fractions, mass_flow=1, T=None, P=None, H=None, eos='pr'):
        """
        Constructor for the Fluid class. This class is used to represent a fluid in the DWSIM simulation.
        Two of T, P or H must be specified.

        Parameters:
            T (float): Temperature of the fluid in Kelvin.
            P (float): Pressure of the fluid in Pascal.
            H (float): Enthalpy of the fluid in J/kg.
            mass_flow (float): Mass flow rate of the fluid in kg/s.
            components (list): List of strings representing the components in the fluid.
            fractions (list): List of floats representing the mole or mass fractions of the components in the fluid.
            eos (str): String representing the equation of state to be used in the simulation.
        """    
    
        if len([i for i in [T, P, H] if i is not None]) != 2:
            raise ValueError("Two of T, P, and H must be specified.")
        
        self.dtlc = CalculatorInterface.Calculator()
        self.dtlc.Initialize()
        eos_dict = {"srk": PropertyPackages.SRKPropertyPackage(True), "nrtl": PropertyPackages.NRTLPropertyPackage(True), \
                    "pr": PropertyPackages.PengRobinsonPropertyPackage(True),\
                    "pr-lk": PropertyPackages.PengRobinsonLKPropertyPackage(True),\
                    "uniquac": PropertyPackages.UNIQUACPropertyPackage(True),\
                    "unifac": PropertyPackages.UNIFACPropertyPackage(True)} # NOTE: there are more EoS available in DWSIM. To add to this dictionary, call Fluid.get_valid_eos() and add to this dictionary.
        self.eos = eos_dict[eos.lower()]
        self.dtlc.TransferCompounds(self.eos)
        components = [comp[0].upper() + comp[1:].lower() for comp in components] # lower case all components added except for the first letter.
        self.components = Array[String](components)
        self.fractions = Array[Double](fractions)
        self.T = T
        self.P = P

        # Create a material stream
        self.ms = self.dtlc.CreateMaterialStream(self.components, self.fractions)
        self.ms.SetPropertyPackage(self.eos)
        self.ms.SetPressure(P)
        self.ms.SetTemperature(T)
        self.ms.SetMassFlow(mass_flow) # kg/s
        if T is not None and P is not None:
            self.ms.SetFlashSpec("PT")
        elif T is not None and H is not None:
            self.ms.SetFlashSpec("TH")
        elif P is not None and H is not None:
            self.ms.SetFlashSpec("PH")
        
        self.results = self.ms.Calculate()

    def get_valid_eos(self):
        """
        Method to get the list of available equations of state in the DWSIM simulation.

        Returns:
            list: List of strings representing the available equations of state.
        """
        proppacks = self.dtlc.GetPropPackList()
        return [eos for eos in proppacks]
    
    def get_valid_components(self):
        """
        Method to get the list of available components in the DWSIM simulation.
        Components entered must be exactly matched to the spelling in the list returned
        by this method.

        Returns:
            list: List of strings representing the available components.
        """
        compounds = self.dtlc.GetCompoundList()
        return [comp for comp in compounds]
    
    def get_compound_const_props(self, compound):
        """
        Method to get the constant properties of a compound in the fluid.

        Parameters:
            compound (str): String representing the compound for which properties are to be retrieved.
        Returns:
            dict: Dictionary containing the constant properties of the compound.
        """
        if compound not in self.components:
           raise Exception("Compound not found in the fluid. Components in fluid are: " + ", ".join(self.components))

        compprops = self.dtlc.GetCompoundConstPropList()
        compound_props = {}
        for prop in compprops:
            pval = self.dtlc.GetCompoundConstProp(compound, prop)
            compound_props[prop] = pval
        return compound_props
    
    def get_compound_pdep_props(self, compound):
        """
        Method to get the pressure-dependent properties of a compound in the fluid.

        Parameters:
            compound (str): String representing the compound for which properties are to be retrieved.
        Returns:
            dict: Dictionary containing the pressure-dependent properties of the compound.
        """
        if compound not in self.components:
            raise Exception("Compound not found in the fluid. Components in fluid are: " + ", ".join(self.components))

        compprops = self.dtlc.GetCompoundPDepPropList()
        compound_props = {}
        for prop in compprops:
            pval = self.dtlc.GetCompoundPDepProp(compound, prop, self.P)
            compound_props[prop] = pval
        return compound_props

    def get_compound_tdep_props(self, compound):
        """
        Method to get the temperature-dependent properties of a compound in the fluid.

        Parameters:
            compound (str): String representing the compound for which properties are to be retrieved.
        Returns:
            dict: Dictionary containing the temperature-dependent properties of the compound.
        """
        if compound not in self.components:
            raise Exception("Compound not found in the fluid. Components in fluid are: " + ", ".join(self.components))

        compprops = self.dtlc.GetCompoundTDepPropList()
        compound_props = {}
        for prop in compprops:
            pval = self.dtlc.GetCompoundTDepProp(compound, prop, self.T)
            compound_props[prop] = pval
        return compound_props

    def change_temperature(self, T):
        """
        Method to change the temperature of the fluid.

        Parameters:
            T (float): New temperature of the fluid in Kelvin.
        """
        self.T = T
        self.ms.SetTemperature(T)
        self.results = self.ms.Calculate()
    
    def change_pressure(self, P):
        """
        Method to change the pressure of the fluid.

        Parameters:
            P (float): New pressure of the fluid in Pascal.
        """
        self.P = P
        self.ms.SetPressure(P)
        self.results = self.ms.Calculate()


    def get_phase_composition(self, basis='Mass'):
        """
        Method to get the composition of the phases in the fluid.

        Parameters:
            basis (str): String representing the basis of the composition. Can be 'Mole' or 'Mass'.
        Returns:
            dict: Dictionary containing the overall, vapor, liquid, and liquid2 compositions of the fluid.
        """

        overall_composition, vapor_composition,\
             liquid_composition, liquid2_composition = {}, {}, {}, {}
        
        for comp in self.ms.GetPhase("Overall").Compounds.Values:
            overall_composition[comp.ConstantProperties.Name] = comp.MoleFraction if basis == 'Mole' else comp.MassFraction

        for comp in self.ms.GetPhase("Vapor").Compounds.Values:
            vapor_composition[comp.ConstantProperties.Name] = comp.MoleFraction if basis == 'Mole' else comp.MassFraction

        for comp in self.ms.GetPhase("Liquid1").Compounds.Values:
            liquid_composition[comp.ConstantProperties.Name] = comp.MoleFraction if basis == 'Mole' else comp.MassFraction
        
        for comp in self.ms.GetPhase("Liquid2").Compounds.Values:
            liquid2_composition[comp.ConstantProperties.Name] = comp.MoleFraction if basis == 'Mole' else comp.MassFraction

        return {"TwoPhase": overall_composition, "Vapor": vapor_composition,\
                "Liquid": liquid_composition, "Liquid2": liquid2_composition}


    def get_mixture_properties(self):
            """
            Method to get the mixture properties of the fluid. Performs a flash calculation
            using the specified EoS at T, P and returns properties of the vapour phase and
            the liquid phase.

            Returns:
                dict: Dictionary containing the mixture properties of the fluid.
            """

            vapor_properties = {}
            liquid_properties = {}
            liquid2_properties = {}
            overall_properties = {}

            vapor_phase = self.ms.GetPhase("Vapor").Properties
            liquid_phase = self.ms.GetPhase("Liquid1").Properties
            liquid2_phase = self.ms.GetPhase("Liquid2").Properties
            overall = self.ms.GetPhase("Overall").Properties 

            for attr in dir(vapor_phase):
                if not attr.startswith("_"):
                    try:
                        vvalue = getattr(vapor_phase, attr)
                        lvalue = getattr(liquid_phase, attr)
                        l2value = getattr(liquid2_phase, attr)
                        if callable(vvalue): continue
                        if vvalue is not None:
                            vapor_properties[attr] = vvalue
                        if lvalue is not None:
                            liquid_properties[attr] = lvalue
                        if l2value is not None:
                            liquid2_properties[attr] = l2value
                    except Exception as e:
                        print(f"Could not fetch {attr}: {e}")

            for attr in dir(overall):
                if not attr.startswith("_"):
                    try:
                        value = getattr(overall, attr)
                        if callable(value): continue
                        if value is not None:
                            overall_properties[attr] = value
                    except Exception as e:
                        print(f"Could not fetch {attr}: {e}")

            return {"Vapor": vapor_properties, "Liquid": liquid_properties,\
                    "Liquid2": liquid2_phase, "Overall": overall_properties}