module type_radial_basis_function
    use mod_kinds, only: ik, rk
    implicit none

    type, abstract, public :: radial_basis_function_t

        character(len=:), allocatable, private :: name_

    contains


        procedure                               :: set_name
        procedure                               :: get_name

        procedure(init_interface),      deferred    :: init             !< Initialize function and register options
        procedure(compute_interface), deferred :: compute
        procedure(compute_grad_interface), deferred :: compute_grad

    end type radial_basis_function_t

    abstract interface 
        function compute_interface(self, eval_node, support_node, support_radius)
            use mod_kinds, only: rk
            import radial_basis_function_t

            class(radial_basis_function_t), intent(inout) :: self
            real(rk),                       intent(in)      :: eval_node(3)
            real(rk),                       intent(in)      :: support_node(3)
            real(rk),                       intent(in)      :: support_radius(3)
            real(rk)                                        :: compute_interface
        end function
    end interface

    abstract interface 
        function compute_grad_interface(self, eval_node, support_node, support_radius)
            use mod_kinds, only: rk
            import radial_basis_function_t

            class(radial_basis_function_t), intent(inout) :: self
            real(rk),                       intent(in)      :: eval_node(3)
            real(rk),                       intent(in)      :: support_node(3)
            real(rk),                       intent(in)      :: support_radius(3)
            real(rk)                                        :: compute_grad_interface(3)
        end function
    end interface

    abstract interface
        subroutine init_interface(self)
            import radial_basis_function_t

            class(radial_basis_function_t),  intent(inout)  :: self

        end subroutine
    end interface




contains

    subroutine set_name(self, rbf_name)
        class(radial_basis_function_t), intent(inout) :: self
        character(*),                   intent(in)    :: rbf_name

        self%name_ = rbf_name
        
    end subroutine set_name


    function get_name(self) result(rbf_name)
        class(radial_basis_function_t), intent(inout) :: self
        character(len=:), allocatable                   :: rbf_name

        rbf_name = self%name_

    end function get_name
end module type_radial_basis_function
