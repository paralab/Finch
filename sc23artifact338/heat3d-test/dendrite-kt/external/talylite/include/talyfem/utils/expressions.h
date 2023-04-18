#pragma once

#include <petsc.h>
#include <string>
#include <ostream>
#include <limits>

/**
 * This file provides some wrappers for the exprtk library, which can be used to evaluate arbitrary
 * mathematical expressions at run-time (e.g. loaded from a config file).
 * We don't use the exprtk types directly here, as exprtk is completely templated (!) and provided
 * as a single 1.4mb header file (!!), which causes compilation to take upwards of 10 seconds for evey
 * file that includes it (!!!). To get around this long compile time, we only forward-declare the exprtk
 * types here and and only include exprtk.hpp in a single cpp file (utils/expressions.cpp).
 *
 * This also is why we use pointers to exprtk objects (since we only forward-declare the exprtk classes,
 * the compiler doesn't know the size of the objects, so we can't put the objects directly inside our
 * wrapper classes).
 *
 * Note that this only works because we use only 1 specialization of exprtk (on PetscScalar), otherwise
 * we would have to template our wrapper classes, which would require us to bring exprtk.hpp in the header,
 * solving nothing. :)
 *
 * These wrapper classes are very bare-bones, so feel free to add extra wrapper methods as you need them.
 */

// forward declarations for exprtk classes (so we can declare pointers)
namespace exprtk {

template <typename T>
class symbol_table;

template <typename T>
class expression;

}

namespace TALYFEMLIB {

// forward declaration for friend statement
class Expression;

/**
 * Wrapper for exprtk's symbol_table class. Maps symbols (variables, constants, functions, ...) to values.
 * The current interface only supports variables (you can't register custom functions).
 *
 * There are two options for setting variables:
 *
 *   1. Call add_variable("name", &val). As you change val, the variable will automatically update.
 *      This is almost ideal, but requires you to have declared val before your Expressions are created
 *      (e.g. before you parse your InputData).
 *
 *   2. Call create_variable("name"), then call set_variable("name", val) every time you need to change
 *      the variable. This allows you to "declare" variables to the Expression system without actually
 *      carrying a C++ variable all the way from config initialization to where you want to set the variables.
 */
class SymbolTable {
 public:
  SymbolTable();
  virtual ~SymbolTable();

  SymbolTable(const SymbolTable&);
  SymbolTable& operator=(const SymbolTable&);

  /**
   * Register a fixed value with the symbol table. Use of "name" will be automatically substituted with val
   * (for expressions registered with this symbol table).
   * @param name
   * @param val
   */
  void add_constant(const std::string& name, PetscScalar val);

  /**
   * Register a variable with the symbol table. Use of "name" will be automatically be substituted with
   * the current value of val (for expressions registered with this symbol table).
   * @param name variable name, as it should be referred to in expressions
   * @param val value for name
   */
  void add_variable(const std::string& name, PetscScalar* val);

  /**
   * Create a variable managed by the SymbolTable.
   * You can later call "table.variable_ref(name) = 123" to change the variable.
   * @param name variable name, as it should be referred to in expressions
   * @param initial_value initial value, defaults to NaN
   */
  void create_variable(const std::string& name,
      PetscScalar initial_value = std::numeric_limits<PetscScalar>::signaling_NaN());

  /**
   * Get a reference to the variable registered as name. If you change this reference, the variable will change.
   * Commonly used with create_variable.
   * NOTE: If the name is invalid (e.g. for a variable that has not been add_variable'd or create_variable'd),
   *       this function returns a reference to a statically allocated "null" variable (and changing it does nothing).
   * @param name variable name
   * @return mutable reference to variable
   */
  PetscScalar& variable_ref(const std::string& name);

  /**
   * Helper method for changing the value of a variable by name.
   * Particularly useful for variables created with create_variable().
   * @param name variable name
   * @param new_val value to set the variable to
   */
  inline void set_variable(const std::string& name, PetscScalar new_val) {
    variable_ref(name) = new_val;
  }

 protected:
  friend Expression;
  exprtk::symbol_table<PetscScalar>* table_;  ///< exprtk symbol_table we are wrapping
};

/**
 * Wrapper for exprtk's expression class.
 * To use:
 * 1. Create an expression object with the default constructor.
 * 2. Optionally, call register_symbol_table() to register additional symbol tables.
 * 3. Make sure you have filled in the symbol table (i.e. call table.add_variable() or table.create_variable(),
 *    see the SymbolTable class comment for more info on the difference between the two).
 * 4. Call set_expression() to specify what expression to evaluate. This will compile the expression.
 * 5. Call value (probably many times).
 */
class Expression {
 public:
  /**
   * Create a "blank" expression, pre-registered with the SymbolTable global_symbol_table().
   */
  Expression();

  /**
   * Helper for fully constructing an expression that's ready to be evaluated immediately.
   * @param expr expression string
   */
  explicit Expression(const std::string& expr);

  Expression(const Expression&);  // copy constructor
  Expression& operator= (const Expression&);  // assignment operator
  virtual ~Expression();

  /**
   * This returns a "global" symbol table that is automatically registered for every expression.
   * You don't have to use this, but it may be more convenient.
   * @return global symbol table (pre-registered with every equation)
   */
  static SymbolTable& global_symbol_table();

  /**
   * Register a new symbol table with the expression object. This should be done before calling set_expression().
   * It is possible to register multiple symbol tables at once. Conflicts are resolved based on order
   * of registration (see exprtk docs).
   * @param table symbol table to register
   */
  void register_symbol_table(const SymbolTable& table);

  /**
   * Set or change the expression to evaluate. This should be called after register_symbol_table.
   * @param expr_str new expression
   * @throws TALYException on parse error
   */
  void set_expression(const std::string& expr_str);

  /**
   * Evaluate the expression.
   * @return value value of expression
   */
  PetscScalar value() const;

  friend std::ostream& operator<<(std::ostream& stream, const Expression& expr);

 protected:
  static SymbolTable g_global_symbol_table_;
  exprtk::expression<PetscScalar>* expr_;  ///< exprtk expression we are wrapping
};

std::ostream& operator<<(std::ostream& stream, const Expression& expr);

}  // namespace TALYFEMLIB

