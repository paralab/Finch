#include <talyfem/utils/expressions.h>
#include <talyfem/utils/utils.h>  // for PrintError
#include <assert.h>

#include <exprtk.hpp>
#include <talyfem/common/exceptions.h>

namespace TALYFEMLIB {

SymbolTable::SymbolTable() {
  table_ = new exprtk::symbol_table<PetscScalar>();
  table_->add_constants();  // pi, etc
}

SymbolTable::~SymbolTable() {
  delete table_;
}

void SymbolTable::add_variable(const std::string& name, PetscScalar* val) {
  if (!table_->add_variable(name, *val))
    throw TALYException() << "Failed to add variable '" << name << "' (invalid symbol? registering twice?)";
}

void SymbolTable::add_constant(const std::string& name, PetscScalar val) {
  table_->add_constant(name, val);
}

SymbolTable::SymbolTable(const SymbolTable& rhs) {
  table_ = new exprtk::symbol_table<PetscScalar>(*rhs.table_);
}

SymbolTable& SymbolTable::operator=(const SymbolTable& rhs) {
  delete table_;
  table_ = new exprtk::symbol_table<PetscScalar>(*rhs.table_);
  return *this;
}

void SymbolTable::create_variable(const std::string &name, PetscScalar initial_value) {
  table_->create_variable(name, initial_value);
}

PetscScalar& SymbolTable::variable_ref(const std::string &name) {
  return table_->variable_ref(name);
}


SymbolTable Expression::g_global_symbol_table_;

Expression::Expression() {
  expr_ = new exprtk::expression<PetscScalar>();
  register_symbol_table(global_symbol_table());
}

Expression::Expression(const std::string &expr_str) : Expression() {
  set_expression(expr_str);
}

Expression::Expression(const Expression& rhs) {
  expr_ = new exprtk::expression<PetscScalar>(*rhs.expr_);
}

Expression &Expression::operator=(const Expression& rhs) {
  delete expr_;
  expr_ = new exprtk::expression<PetscScalar>(*rhs.expr_);
  return *this;
}

Expression::~Expression() {
  delete expr_;
}

PetscScalar Expression::value() const {
  return expr_->value();
}

void Expression::set_expression(const std::string &expr_str) {
  exprtk::parser<PetscScalar> parser;
  if (!parser.compile(expr_str, *expr_)) {
    PrintError(parser.error());

    // print detailed error info as per exprtk readme
    for (std::size_t i = 0; i < parser.error_count(); i++) {
      typedef exprtk::parser_error::type error_t;

      error_t error = parser.get_error(i);
      PrintError("Error[", i, "] Position: ", error.token.position,
          " Type: [", exprtk::parser_error::to_str(error.mode),
          "] Msg: ", error.diagnostic);
    }
    throw TALYException() << "Error parsing expression";
  }
}

void Expression::register_symbol_table(const SymbolTable &table) {
  expr_->register_symbol_table(*table.table_);
}

SymbolTable &Expression::global_symbol_table() {
  return g_global_symbol_table_;
}

std::ostream& operator<<(std::ostream &stream, const Expression& expr) {
  return stream << "[expression]";
}

}  // namespace TALYFEMLIB
