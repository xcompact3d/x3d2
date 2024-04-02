Fortran style guide
===================

Basic formatting
----------------

Identation, whitespaces and line length are enforced by `fprettify`.
Project settings are defined in the `.fprettify.ini` located at the
root of the repository.

In most cases, developers are expected to install the associated
pre-commit git hook.  Additionally, it is recommended to configure
your text editor to run fprettify for you (for instance on saving a
file).

If you want to preview changes without overwriting a file, you can do
so using the `--stdout` option:

.. code:: console

   $ fprettify --config .fprettify.ini --stdout

Note: CUDA Fortran chevron syntax is not supported by `fprettify`.
Thus, we use `!&` to deactivate `fprettify` on lines chevron syntax is
used.

.. code:: fortran

   call gpu_kernel<<<blocks, threads>>>(args, ...) !&

Naming conventions
------------------

In the following "symbol" is a catch all phrase for variables,
procedures, derived types, type components and interfaces.

- Symbols are name using the `snake_case` convention.
- Variable declared with the `parameter` are named using UPPER CASE
  letters.
- Symbols should, as much as possible, be named after they meaning
  rather than their mathematical notation (.e.g `velocity_field`
  instead of `u`).

Procedure definitions
---------------------

- Procedure prototypes spanning more than 79 characters should be split
  before the first dummy argument is defined, after the opening
  parenthese.  If the list of dummy arguments spans more than 79
  characters each argument is defined on it own line.

  .. code:: fortran

     subroutine my_long_subroutine( &
        argument1, argument3, argument4, argument5, argument5, argument6 &
	)

     subroutine my_very_long_subroutine( &
        argument1, &
	argument2, &
	argument3, &
	argument4, &
	argument5, &
	argument6, &
	argument7, &
	)
- Function prototypes indicate the return object type unless the
  return object is an array.
- Functions are always defined with the `pure` prefix.  A procedure
  with side effects must be a `subroutine`.

Modules
-------

- Module names end with the `_m` suffix. Examples:
  `module allocator_m`, `use allocator_m, only: ...`.
- All module components are `private` by default.
- The module components are always acessed through the `only` keyword:

  .. code:: fortran

     module a_module_m
        ! Non compliant
	use stencil_m
	! Compliant
	use stencil_m, only: stencil_t

Derived type definitions
------------------------

- Derived type names end with the `_t` suffix. Examples:
  `allocator_t`, `type(stencil_t) :: s`. 
- The `contains` keyword is omitted is the type does not define any
  type-bound procedures.
- All type components are `private` by default.

Custom structure constructors
-----------------------------

- Are named like `create_<type_root_name>[_<suffix>]`
- Are declared with the `private` attribute
- Are defined at the top of the module's `contains` block.

Example

.. code:: fortran

   module square_module

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
         type(square), intent(in) :: sq_in
         ! ...
      end function create_square_from_square

      type(square_t) function create_square_default_color(sq_size)
         real, intent(in) :: sq_size
         ! ...
      end function create_square_default_color

.. _in-code-docs:

In-code documentation
---------------------

The body of modules, public types, public procedures and public
type-bound methods MUST be preceded of one or more documentation
paragraphs.  Optionally, the body of private symbols MAY be
preceded by documentation paragraph.

Procedure dummy arguments, interface components and type-bound
procedures declarations MAY be documented using an inline comment
either on the same line directly following the statement (using the
docmark `!!`) or on the line directly above the statement (using the
predocmark `!>`).

x3d2 uses `ford <https://forddocs.readthedocs.io/en/latest/>`_ to
extract in-code documentation and generate HTML pages.  The syntax for
in-code documentation follows ford's syntax for comments. See
`(ford)Writing documentation
<https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html>`_

.. code:: fortran

   subroutine add(a, b, c)
       !! This is the first paragraph of the procedures
       !! documentation.  Note that it starts with TWO !.
       real, intent(in) :: a, b !! Optional documentation for dummy argument.
       real, intent(out) :: c !! The result of a + b

       ! The line below is a regular comment.
       ! Make use of the well-known the addition algorithm.
       c = a + b
   end subroutine
