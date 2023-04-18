Style Guide
===========
All new code in the library should be written in a consistent style. The purpose is to make the code easier to understand and maintain. The rules below describe the proper style of the library. All new code for the library should follow these guidelines. If you're editing existing library code that does not match the style guide (especially naming issues), first update it to the new style before making any changes. Add the style updates as a separate commit. *Do not* mix style updates and functionality changes in the same commit (doing this makes it easy to miss important changes).

General Style Guide
-------------------
The style to use for the code is given in the [google c++ style guide](https://google.github.io/styleguide/cppguide.html). All rules from that guide should be followed when writing new library code, unless doing so would be impractical or if there is an exception listed below.

Library Specific Style Guide
----------------------------
There are several additional style items not listed in the google guide that
apply to this library. They are:

* Exception to google guide: We do not enforce the rule against non-const
reference parameters. Function parameters that are references need not be const.

* Exception to google guide: "#pragma once" is allowed in place of the include
guard given in the style guide

* Do not include multiple statements on a single line. This is less clear and
is prone to error. There are two exceptions:
The first is the use of the PETSc functions CHKERRQ, CHKERRV, etc... for
checking error codes. These may be on the same line as the function returning
the error (provided that the maximum line length is not exceeded).
The second is concerning if statements. If blocks containing *only* one statment
may be placed on the same line as the if statement (but do not need to be).
```c++
double d1; int a = 17;  // NO!
int value = calc_value(a); a++;  // NO!
switch (value) {
  case 1: output_string = "one"; break;  // NO!
  case 2: { output_string = "two"; break; }  // still NO!
  case 4: {  // YES. this is the correct way
    output_string = "four";
    break;
  }
}
ierr = VecAssemblyBegin(*vec); CHKERRQ(ierr);  // ok (exception noted above)
if (d1 > 7.5) { update_array(); }  // ok (exception noted above)
if (d1 < 0.0) { update_array(); return; }  // NO! (has two statements in block)
```

* The code following an if statement, else statement, for loop, switch
statement, while loop, etc ... *must* be contained in braces ({ }), even if
there is only a single line. This is for clarity and to avoid future problems
when code is edited.
```c++
// these are wrong
for (int i = 0; i < 10; i++)  // NO! ...
  array[i] = i;  // this should be in braces
for (int i = 0; i < 10; i++) array[i] = i;  // also NO! (missing {})
if (j == 2)  // NO! ...
k *= -1;  // this should be in braces
if (j == 2) k *= -1;  // still NO! (missing {})
// these are right
for (int i = 0; i < 10; i++) {  // OK. this is the proper way
  array[i] = i;
}
for (int i = 0; i < 10; i++) { array[i] = i; }  // OK, note spaces around { }
if (j == 2) {  // OK
  k *= -1;
}
if (j == 2) { k *= -1; }  // also OK, note the spaces around { }
```

Testing Style
=============

Some of the existing code may not match this style, but all new code should. Pay particular attention to the rules for spacing (2 space indent, no tabs) and naming rules. If you followed the installation instructions, you can get a list of some places where your code doesn't match the style guide, by running the following in the library directory:

```
make lint
```
Fix any issues listed before committing your changes. Note that 'make lint' will not catch all of the google style errors and it will not catch any of the library specific style issues.

Documentation
=============
In addition to matching the style guide all classes, functions, and member variables should have comments explaining class. These should be written in [doxygen format](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html) using the style shown below. This will allow automatic generation of documentation for the library, so it is important to get the format correct. Before committing any code, check for documentation errors:

```
make docs
```

(this command only works if the `doxygen` program was installed prior to configuring TalyFEM)

Fix any errors prior to committing the code. Note that currently there are many false errors inside the basis function code due to ridiculous templating - it is safe to ignore these.

In addition to doxygen comments, add comments for anything that isn't obvious. If someone reading the code would not immediately understand what the code does by looking at it, you need to add a comment to explain it. Explain what you are doing and why.

Here is an example of the proper documentation format for a class. In onder to
match the code style, keep comment lines to less than 80 characters in length.

```
/**
 * A test class to show comments.
 *
 * A more detailed description will explain what the class does and how to use
 * it properly. The purpose of this section is to help someone who is unfamiliar
 * with this class quickly understand what it does and how to use it. Include
 * short examples of how the class is used if they would be helpful and do not
 * take up too many lines of code.
 *
 * Use multiple paragraphs if needed. 
 *
 * The format for this block of comments is the same format that will be
 * used for all comment blocks:
 * There is a short *one line* and *one sentence* description of the object,
 * followed by a single blank line, and then a longer multiline description of
 * the object. If a single line description is sufficient to explain the object,
 * just use that and do not include the blank line or the longer description.
 *
 * Markdown format is supported. See details at:
 * https://www.stack.nl/~dimitri/doxygen/manual/markdown.html for more info.
 *
 * Latex formulas are supported, but should only be used if it is not possible
 * to describe it without latex. For details of latex support, see:
 * https://www.stack.nl/~dimitri/doxygen/manual/formulas.html
 *
 * If this class implements a numeric method or something similar, provide
 * references a user can consult for more information.
 */
class Test {
 public:
  /**
   * An enum of type values.
   *
   * More detailed enum description goes here. Again note that we start with a
   * one line summary of the item, leave a blank line, and then have a more
   * detailed description. The detailed description is not needed IF the one
   * line description is sufficient to explain the purpose of the item.
   *
   * Each member of the enum has its own comment in the format shown below.
   * That format is used for enum members and class variables, as shown at the
   * end of this example class.
   */
  enum EnumType {
    TypeVal1,  ///< meaning of value 1 (this goes on one line if possible)
    TypeVal2,  ///< meaning of value 2 (if this is too long to fit on a single
               ///< line, it can be spread out over multiple lines like this).
    TypeVal3,  ///< meaning of value 3.
  }

  /**
   * Constructs the object using the default values.
   *
   * Detailed description of what the constructor does. Note that if the
   * constructor doesn't do anything of interest, there is no need to have any
   * comment at all for the constructor.
   */
  Test();

  /**
   * Destroys the object and all its parts.
   *
   * More details about the destructor go here. Again, if this doesn't do
   * anything of interest, there is no need to have a comment.
   */
  ~Test();

  /**
   * Calculates and returns a useful value.
   *
   * A description of what the function is calculating would go here. Also here
   * would be any special values of the parameters or the return value. For
   * example, if the behaviour is different if an_int is equal to zero, that
   * would be mentioned here.
   *
   * @param an_int description of first value used in the calculation.
   * @param a_double description of second value used in the calculation. If
   *                 this is too long to fit on a single line, wrap it to
   *                 multiple lines. Keep the indention consistent.
   * @return Description of the useful value we just calculated (note that there
   *         is no space between list of param values and the return statement)
   */
  int CalcUsefulValue(int an_int, double a_double);

  /**
   * Prints information about the class.
   *
   * A description of what is printed would go here.
   */
  void PrintInfo();

  /**
   * Does a simple thing that can be fully explained in one line.
   */
  void DoSimpleAction();

  int count_;  ///< description of the member variable for the class.
  double complicated_value_;  ///< Variable requiring more than one line to 
                              ///< explain. In this case, we add additional 
                              ///< lines as needed. Note: if there is more than
                              ///< one sentence, the first one should be a
                              ///< brief, general description of the variable.

  double long_variable_name_that_doesnt_leave_space_to_describe_;
  ///< Variable with extremely long name. Since there is not enough space for a
  ///< comment using the end of line style, we can put the comment on the lines
  ///< following the variable. As before, if there is more than one sentence,
  ///< the first one should be a brief, general description of the variable.
}
```