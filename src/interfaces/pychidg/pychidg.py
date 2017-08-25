from __future__ import print_function, absolute_import, division
import _chidg
import f90wrap.runtime
import logging





#class Type_Chidg(f90wrap.runtime.FortranModule):
#    """
#    Module type_chidg
#    
#    
#    Defined at type_chidg.f90 lines 1-1548
#    
#    """
@f90wrap.runtime.register_class("chidg")
class chidg(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=chidg_t)
    
    
    Defined at type_chidg.f90 lines 72-125
    
    """
    def start_up(self, activity):
        """
        start_up(self, activity)
        
        
        Defined at type_chidg.f90 lines 162-244
        
        Parameters
        ----------
        self : Chidg_T
        activity : str
        
        """
        _chidg.f90wrap_start_up(self=self._handle, activity=activity)
    
    def shut_down(self, selection=None):
        """
        shut_down(self[, selection])
        
        
        Defined at type_chidg.f90 lines 261-295
        
        Parameters
        ----------
        self : Chidg_T
        selection : str
        
        """
        _chidg.f90wrap_shut_down(self=self._handle, selection=selection)
    
    
    def set(self, selector, algorithm=None, integer_input=None, real_input=None):
        """
        set(self, selector[, algorithm, integer_input, real_input])
        
        
        Defined at type_chidg.f90 lines 461-562
        
        Parameters
        ----------
        self : Chidg_T
        selector : str
        algorithm : str
        integer_input : int
        real_input : float
        
        """
        _chidg.f90wrap_set(self=self._handle, selector=selector, algorithm=algorithm, \
            integer_input=integer_input, real_input=real_input)

    def project(self, function, ifield):
        """
        project(self, activity)
        
        
        Defined at type_chidg.f90 lines 162-244
        
        Parameters
        ----------
        self : Chidg_T
        activity : str
        
        """
        _chidg.f90wrap_project(chidg=self._handle, func=function._handle, ifield=ifield)
    
    def read_mesh(self, grid_file, interpolation=None, level=None, equation_set=None):
        """
        read_mesh(self, grid_file[, equation_set, bc_wall, bc_inlet, bc_outlet, \
            bc_symmetry, bc_farfield, bc_periodic, interpolation, level])
        
        
        Defined at type_chidg.f90 lines 596-649
        
        Parameters
        ----------
        self : Chidg_T
        grid_file : str
        equation_set : str
        bc_wall : unknown
        bc_inlet : unknown
        bc_outlet : unknown
        bc_symmetry : unknown
        bc_farfield : unknown
        bc_periodic : unknown
        interpolation : str
        level : int
        
        """
        _chidg.f90wrap_read_mesh(self=self._handle, grid_file=grid_file, interpolation=interpolation, level=level, equation_set=equation_set)
    
    def read_fields(self, file_name):
        """
        read_fields(self, file_name)
        
        
        Defined at type_chidg.f90 lines 999-1020
        
        Parameters
        ----------
        self : Chidg_T
        file_name : str
        
        """
        _chidg.f90wrap_read_fields(self=self._handle, file_name=file_name)
    
    def write_mesh(self, file_name):
        """
        write_mesh(self, file_name)
        
        
        Defined at type_chidg.f90 lines 1108-1181
        
        Parameters
        ----------
        self : Chidg_T
        file_name : str
        
        """
        _chidg.f90wrap_write_mesh(self=self._handle, file_name=file_name)
    
    def write_fields(self, file_name):
        """
        write_fields(self, file_name)
        
        
        Defined at type_chidg.f90 lines 1202-1223
        
        Parameters
        ----------
        self : Chidg_T
        file_name : str
        
        """
        _chidg.f90wrap_write_fields(self=self._handle, file_name=file_name)


    def produce_visualization(self, grid_file, solution_file, equation_set=None):
        """
        produce_visualization(self, grid_file, solution_file)
        
        
        Defined at type_chidg.f90 lines 1202-1223
        
        Parameters
        ----------
        self : Chidg_T
        file_name : str
        
        """
        _chidg.f90wrap_produce_visualization(self=self._handle, grid_file=grid_file, solution_file=solution_file, equation_set=equation_set)



    
    def run(self, write_initial=None, write_final=None):
        """
        run(self[, write_initial, write_final])
        
        
        Defined at type_chidg.f90 lines 1310-1444
        
        Parameters
        ----------
        self : Chidg_T
        write_initial : bool
        write_final : bool
        
        """
        _chidg.f90wrap_run(self=self._handle, write_initial=write_initial, \
            write_final=write_final)
    
    def report(self, selection):
        """
        report(self, selection)
        
        
        Defined at type_chidg.f90 lines 1465-1496
        
        Parameters
        ----------
        self : Chidg_T
        selection : str
        
        """
        _chidg.f90wrap_report(self=self._handle, selection=selection)
    
    def __init__(self, handle=None):
        """
        self = Chidg_T()
        
        
        Defined at type_chidg.f90 lines 72-125
        
        
        Returns
        -------
        this : Chidg_T
            Object to be constructed
        
        
        Automatically generated constructor for chidg_t
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _chidg.f90wrap_chidg_t_initialise()
    
    def __del__(self):
        """
        Destructor for class Chidg_T
        
        
        Defined at type_chidg.f90 lines 72-125
        
        Parameters
        ----------
        this : Chidg_T
            Object to be destructed
        
        
        Automatically generated destructor for chidg_t
        """
        if self._alloc:
            _chidg.f90wrap_chidg_t_finalise(this=self._handle)
    
    
    _dt_array_initialisers = []

        
#        def __str__(self):
#            ret = ['<chidg_t>{\n']
#            ret.append('    auxiliary_environment : ')
#            ret.append(repr(self.auxiliary_environment))
#            ret.append(',\n    nterms_s : ')
#            ret.append(repr(self.nterms_s))
#            ret.append(',\n    nterms_s_1d : ')
#            ret.append(repr(self.nterms_s_1d))
#            ret.append(',\n    grid_file : ')
#            ret.append(repr(self.grid_file))
#            ret.append(',\n    solution_file_in : ')
#            ret.append(repr(self.solution_file_in))
#            ret.append(',\n    time_integrator : ')
#            ret.append(repr(self.time_integrator))
#            ret.append(',\n    nonlinear_solver : ')
#            ret.append(repr(self.nonlinear_solver))
#            ret.append(',\n    linear_solver : ')
#            ret.append(repr(self.linear_solver))
#            ret.append(',\n    preconditioner : ')
#            ret.append(repr(self.preconditioner))
#            ret.append(',\n    envinitialized : ')
#            ret.append(repr(self.envinitialized))
#            ret.append('}')
#            return ''.join(ret)
#
#    _dt_array_initialisers = []








@f90wrap.runtime.register_class("function")
class function(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=function_t)
    
    
    Defined at type_function.f90
    
    """
    def __init__(self, handle=None):
        """
        self = function_t()
        
        
        Defined at type_function.f90
        
        
        Returns
        -------
        this : function_t
            Object to be constructed
        
        
        Automatically generated constructor for function_t
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _chidg.f90wrap_function_t_initialise()
    
    def __del__(self):
        """
        Destructor for class function_t
        
        
        Defined at type_chidg.f90 lines 72-125
        
        Parameters
        ----------
        this : Chidg_T
            Object to be destructed
        
        
        Automatically generated destructor for chidg_t
        """
        if self._alloc:
            _chidg.f90wrap_function_t_finalise(this=self._handle)
    


    def create_function(self, fcn_name):
        """
        self = function_t()
        
        
        Defined at type_function.f90
        
        
        Returns
        -------
        this : function_t
            Object to be constructed
        
        
        Automatically generated constructor for function_t
        """
        _chidg.f90wrap_create_function(fcn=self._handle,fcn_name=fcn_name)
    

    def set_option(self, key, val):
        """
        self = function_t()
        
        
        Defined at type_function.f90
        
        
        Returns
        -------
        this : function_t
            Object to be constructed
        
        
        Automatically generated constructor for function_t
        """
        _chidg.f90wrap_set_option(self=self._handle,key=key,val=val)

    _dt_array_initialisers = []



@f90wrap.runtime.register_class("equation_set")
class equation_set(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=equation_set_t)
    
    
    Defined at type_chidg.f90 lines 22-33
    
    """
    def set_name(self, ename):
        """
        set_name(self, ename)
        
        
        Defined at type_chidg.f90 lines 52-57
        
        Parameters
        ----------
        self : Equation_Set_T
        ename : str
        
        """
        _chidg.f90wrap_set_name(self=self._handle, ename=ename)
    
    def add_operator(self, string_bn):
        """
        add_operator(self, string_bn)
        
        
        Defined at type_chidg.f90 lines 99-103
        
        Parameters
        ----------
        self : Equation_Set_T
        string_bn : str
        
        """
        _chidg.f90wrap_add_operator(self=self._handle, string_bn=string_bn)
    
    def add_model(self, string_bn):
        """
        add_model(self, string_bn)
        
        
        Defined at type_chidg.f90 lines 121-125
        
        Parameters
        ----------
        self : Equation_Set_T
        string_bn : str
        
        """
        _chidg.f90wrap_add_model(self=self._handle, string_bn=string_bn)

    def add_io_field(self, field):
        """
        add_io_field(self, field)
        
        
        Defined at type_chidg.f90 lines 121-125
        
        Parameters
        ----------
        self : Equation_Set_T
        string_bn : str
        
        """
        _chidg.f90wrap_add_io_field(self=self._handle, field=field)

    def clear_io_fields(self):
        """
        add_io_field(self, field)
        
        
        Defined at type_chidg.f90 lines 121-125
        
        Parameters
        ----------
        self : Equation_Set_T
        string_bn : str
        
        """
        _chidg.f90wrap_clear_io_fields(self=self._handle)
    
    def __init__(self, handle=None):
        """
        self = Equation_Set_T()
        
        
        Defined at type_chidg.f90 lines 22-33
        
        
        Returns
        -------
        this : Equation_Set_T
            Object to be constructed
        
        
        Automatically generated constructor for equation_set_t
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _chidg.f90wrap_equation_set_t_initialise()
    
    def __del__(self):
        """
        Destructor for class Equation_Set_T
        
        
        Defined at type_chidg.f90 lines 22-33
        
        Parameters
        ----------
        this : Equation_Set_T
            Object to be destructed
        
        
        Automatically generated destructor for equation_set_t
        """
        if self._alloc:
            _chidg.f90wrap_equation_set_t_finalise(this=self._handle)
    
    _dt_array_initialisers = []


class equations(f90wrap.runtime.FortranModule):
    """
    Module mod_equations
    
    
    Defined at mod_equations.f90 lines 9-129
    
    """
    @f90wrap.runtime.register_class("equation_set_factory")
    class equation_set_factory(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=equation_set_factory_t)
        
        
        Defined at mod_equations.f90 lines 19-29
        
        """
        def register(self, eqnset):
            """
            produce_equation_set = produce_equation_set(self, eqnstring, blueprint)
            
            
            Defined at mod_equations.f90 lines 75-80
            
            Parameters
            ----------
            self : Equation_Set_Factory_T
            eqnstring : str
            blueprint : str
            
            Returns
            -------
            produce_equation_set : float
            
            """
            _chidg.f90wrap_register_equation_set(factory_handle=self._handle, eqnset_handle=eqnset._handle)



        def __init__(self, handle=None):
            """
            self = Equation_Set_Factory_T()
            
            
            Defined at mod_equations.f90 lines 19-29
            
            
            Returns
            -------
            this : Equation_Set_Factory_T
            	Object to be constructed
            
            
            Automatically generated constructor for equation_set_factory_t
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _chidg.f90wrap_equation_set_factory_t_initialise()
        
        def __del__(self):
            """
            Destructor for class Equation_Set_Factory_T
            
            
            Defined at mod_equations.f90 lines 19-29
            
            Parameters
            ----------
            this : Equation_Set_Factory_T
            	Object to be destructed
            
            
            Automatically generated destructor for equation_set_factory_t
            """
            if self._alloc:
                _chidg.f90wrap_equation_set_factory_t_finalise(this=self._handle)
        
        _dt_array_initialisers = []
        
    
    @classmethod
    def factory(self):
        """
        Element equation_set_factory ftype=type(equation_set_factory_t) \
            pytype=Equation_Set_Factory_T
        
        
        Defined at mod_equations.f90 line 36
        
        """
        equation_set_factory_handle = _chidg.f90wrap_mod_equations__get__equation_set_factory()
        equation_set_factory = mod_equations.equation_set_factory.from_handle(equation_set_factory_handle)
        return equation_set_factory
    
    
    def __str__(self):
        ret = ['<mod_equations>{\n']
        ret.append('    equation_set_factory : ')
        ret.append(repr(self.equation_set_factory))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []


    

mod_equations = equations()
#type_chidg = Type_Chidg()
#type_function = Type_Function()

