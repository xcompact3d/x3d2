Fortran style guide
===================

Basic formatting
----------------

Indentation, whitespaces, and line lengths are enforced by `fprettify`. Project settings are defined in the ``.fprettify.ini`` file located at the root of the repository.

Developers are encouraged to install the associated pre-commit Git hook. Additionally, it is recommended to configure your text editor to run `fprettify` automatically (e.g., on saving a file).

To preview changes without overwriting a file, use the `--stdout` option:

.. code:: bash

   $ fprettify --config .fprettify.ini --stdout

.. warning:: 
   
   CUDA Fortran chevron syntax is not supported by `fprettify`. To handle this, use ``!&`` to deactivate `fprettify` on lines where chevron syntax is used.

   .. code:: fortran
  
      call gpu_kernel<<<blocks, threads>>>(args, ...) !&

Naming conventions
------------------

In the following, "symbol" is a catch-all term for variables, procedures, derived types, type components, and interfaces.

- Use lowercase for all Fortran constructs (e.g., ``do``, ``subroutine``, ``module``).
- Symbols are named using the `snake_case` convention.
- Variables declared with the `parameter` attribute are named using UPPER CASE letters. For all other names use lowercase.
- Symbols should, as much as possible, be named after their meaning rather than their mathematical notation (e.g., ``velocity_field`` instead of ``u``).


Procedure definitions
---------------------

- Procedure prototypes spanning more than 79 characters should be split before the first dummy argument is defined, after the opening parenthesis. If the list of dummy arguments spans more than 79 characters, each argument should be defined on its own line.

  .. code:: fortran

     subroutine my_long_subroutine( &
        argument1, argument2, argument3, argument4, argument5, argument6 &
     )

     subroutine my_very_long_subroutine( &
        argument1, &
        argument2, &
        argument3, &
        argument4, &
        argument5, &
        argument6, &
        argument7 &
     )

- Function prototypes should indicate the return object type unless the return object is an array.

  .. code:: fortran

     pure function compute_area(radius) result(area)
       real :: compute_area
       real, intent(in) :: radius
       real :: area

       area = pi * radius**2
     end function compute_area

- Functions should always be defined with the ``pure`` prefix. A procedure with side effects must be a ``subroutine``.

Modules
-------

- Module names should start with the `m_` suffix. Examples: ``module m_allocator``, ``use m_allocator, only: ...``.
- All module components should be ``private`` by default.
- Always access module components using the ``only`` keyword:

  .. code:: fortran

     module m_a_module
        implicit none
        private
        public :: some_subroutine

        contains

        subroutine some_subroutine()
            ! Subroutine implementation
        end subroutine some_subroutine
     end module m_a_module

  .. code:: fortran

     ! Non-compliant
     use m_stencil

     ! Compliant
     use m_stencil, only: stencil_t

Derived type definitions
------------------------

- Derived type names should end with the ``_t`` suffix. Examples: ``allocator_t``, ``type(stencil_t) :: s``.
- Omit the ``contains`` keyword if the type does not define any type-bound procedures.
- All type components should be ``private`` by default.

  .. code:: fortran

     type :: allocator_t
        private
        ! Type components
     end type allocator_t

     type :: stencil_t
        private
        ! Type components
        contains
        procedure :: some_procedure
     end type stencil_t

Custom structure constructors
-----------------------------

- Name constructors as ``create_<type_root_name>[_<suffix>]``.
- Declare constructors with the ``private`` attribute.
- Define constructors at the top of the module's ``contains`` block.

Example:

.. code:: fortran

   module square_module
      implicit none
      private
      public :: square_t, create_square_from_square, create_square_default_color

      type :: square_t
         real :: size
         character(:), allocatable :: color
      end type square_t

      interface square_t
         module procedure create_square_from_square
         module procedure create_square_default_color
      end interface square_t

   contains

      type(square_t) function create_square_from_square(sq_in)
         type(square_t), intent(in) :: sq_in
         ! Function implementation
         create_square_from_square%size = sq_in%size
         create_square_from_square%color = sq_in%color
      end function create_square_from_square

      type(square_t) function create_square_default_color(sq_size)
         real, intent(in) :: sq_size
         ! Function implementation
         create_square_default_color%size = sq_size
         create_square_default_color%color = 'blue'
      end function create_square_default_color

   end module square_module

.. _in-code-docs:

In-code documentation
---------------------

x3d2 uses `FORD <https://forddocs.readthedocs.io/en/latest/>`_ to extract in-code documentation and generate HTML pages. The syntax for in-code documentation follows FORD's syntax for comments. See the `FORD User Guide <https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html>`_ for more details.

The body of modules, public types, public procedures, and public type-bound methods MUST be preceded by one or more documentation paragraphs. Optionally, the body of private symbols MAY be preceded by a documentation paragraph.

Procedure dummy arguments, interface components, and type-bound procedure declarations MAY be documented using an inline comment either on the same line directly following the statement (using the docmark ``!!``) or on the line directly above the statement (using the predocmark ``!>``).

Including LaTeX in in-code documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can include LaTeX equations in your documentation. For inline math, use ``\( ... \)``. For displayed equations, you can use either ``$$ ... $$`` or ``\[ ... \]``. Note that ``$ ... $`` is not supported for inline math. Displayed equations can be written using ``$$ ... $$``, but this method does not number the equations. To create numbered equations, use the ``\begin{equation} ... \end{equation}`` environment. You can also use ``\label{eq:some_equation}`` to label the equations and ``\eqref{eq:some_equation}`` to reference them within the text.

Example:

.. code:: fortran

   subroutine add(a, b, c)
       !! This is the first paragraph of the procedure's
       !! documentation. Note that it starts with TWO !.
       !! The addition operation is defined in-line as \( c = a + b \).
       !! 
       !! The following operation shows it as a displayed equation:
       !! $$c = a + b$$
       !! 
       real, intent(in) :: a, b !! Optional documentation for dummy argument.
       real, intent(out) :: c !! The result of \( a + b \)

       ! The line below is a regular comment.
       ! Make use of the well-known addition algorithm.
       c = a + b
   end subroutine add