import csv
from itertools import combinations
import os
import pandas as pd
import sympy as sp
from sympy import Eq, simplify, Add, S, Symbol
import re
from difflib import SequenceMatcher
from functools import lru_cache
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
from abc import ABC, abstractmethod
from typing import List, Tuple, Set, Dict
import logging
from sympy.parsing.latex import parse_latex

#def parse_latex(latex_str: str) -> List[Eq]:
#    """Handle LaTeX expressions with logical operators"""
#    from sympy.parsing.latex import parse_latex
#    import re
#    
#    # Split on logical operators (\vee, \wedge, etc.)
#    split_expr = re.split(r'\\vee|\\wedge', latex_str)
#    
#    equations = []
#    for expr in split_expr:
#        try:
#            eq = parse_latex(expr.strip())
#            equations.append(eq)
#        except Exception as e:
#            print(f"Could not parse: {expr}")
#    
#    return equations

class LinearEquation:
    """Find linear equations"""
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self._variables = sp.symbols('s r b c w n y m d v o p q x a f t')
    
    def is_linear(self, expr: sp.Expr) -> bool:
        """Check if an expression is linear while also ensuring no fractions or square roots"""
        try:
            poly = expr.as_poly()
            if poly is None or poly.degree() != 1:
                return False
            
            # Check for fractions
            if any(term.is_Pow and term.exp.is_Number and term.exp < 0 for term in expr.atoms(sp.Pow)):
                return False

            # Check for square roots
            if any(term.is_Pow and term.exp == sp.Rational(1, 2) for term in expr.atoms(sp.Pow)):
                return False

            return True
        except Exception as e:
            self.logger.error(f"Linear check failed: {str(e)}")
            return False

class EquationPreprocessor:
    """Handles preprocessing of equation strings"""
    _MULTIPLY_PATTERNS = [
    (re.compile(r'(\d)([a-zA-Z(])'), r'\1*\2'),
    (re.compile(r'([a-zA-Z])(\d)'), r'\1*\2'),
    (re.compile(r'([a-zA-Z])([a-zA-Z])'), r'\1*\2'),
    (re.compile(r'\)([a-zA-Z(])'), r')*\1'),
    (re.compile(r'([a-zA-Z)])(\()'), r'\1*\2'),
]

    def preprocess(self, expr_str: str) -> str:
        """Normalize and insert implicit multiplication operators"""
        expr_str = expr_str.replace('−', '-').strip()
        for pattern, replacement in self._MULTIPLY_PATTERNS:
            expr_str = pattern.sub(replacement, expr_str)
        return expr_str
        

class EquationParser:
    """Handles parsing and caching of equations"""
    def __init__(self, preprocessor: EquationPreprocessor):
        self.preprocessor = preprocessor

    @lru_cache(maxsize=128)
    def parse(self, eq_str: str) -> Tuple[Eq, Eq]:
        """Parse equation string into evaluated and raw forms"""
        try:
            lhs, rhs = map(str.strip, eq_str.split('=', 1))
        except ValueError:
            raise ValueError("Equation must contain exactly one '='")

        lhs_expr = self._parse_side(lhs)
        rhs_expr = self._parse_side(rhs)

        return (
            Eq(simplify(lhs_expr), simplify(rhs_expr)),
            Eq(lhs_expr, rhs_expr,  evaluate=False)
        )

    def _parse_side(self, side_str: str) -> sp.Expr:
        """Parse a single side of an equation"""
        return sp.parse_expr(
            self.preprocessor.preprocess(side_str),
            evaluate=False
        )

class TransformationsGenerator:
    """Generates equation transformations with transformation tracking"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def generate(self, equation: Eq) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Generate all possible equation transformations with transformation tracking"""
        transformations = set()
        try:
            # Add original equation with empty transformation history
            canonical = simplify(equation.lhs - equation.rhs)
             # Extract terms from LHS and subtract RHS terms manually
            lhs_terms = equation.lhs.args if equation.lhs.func == sp.Add else (equation.lhs,)
            rhs_terms = equation.rhs.args if equation.rhs.func == sp.Add else (equation.rhs,)

            # Combine terms: LHS - RHS = LHS + (-RHS_terms)
            unevaluated_terms = list(lhs_terms) + [-term for term in rhs_terms]
            canonical_expr = sp.Add(*unevaluated_terms, evaluate=False)

            transformations.add((equation, tuple(), tuple()))
            transformations.update(self._generate_term_combinations(canonical))
            transformations.update(self._generate_term_combinations(canonical_expr))
            transformations.update(self.generate_all_transformations(equation))

            transformations.update(self._divid_singel_term_by_coff(transformations))
            
            # Get all terms in canonical form (lhs - rhs = 0)
            terms = self._get_terms(equation)
            
            # Generate moves for each term
            for term in terms:
                moved_equations = self._move_term(equation, term, terms)
                for eq, actions, _terms in moved_equations:
                    transformations.add((eq, actions, _terms))
                #    equation = eq  #
            
            # Generate divisions for each variable
            variables = equation.free_symbols
            new_set = set()
            for eq, actions, terms in transformations:
                if len(Add.make_args(eq.lhs)) ==1:
                    for var in variables:
                        if eq.lhs.has(var):
                            coeff = eq.lhs.as_coefficients_dict().get(var, None)
                            if coeff and coeff != 1:
                                new_set.update(self._divide_factor(eq, var, actions, terms))
            
            transformations.update(new_set)



            for var in variables:
            #    transformations.update(self._divide_factor(equation, var))
                transformations.update(self._solve_for_variable(equation, var))

            # Convert to decimal representations
            new_set_decimal = set()
            for eq, actions, terms in transformations:
                new_set_decimal.add(self._convert_to_decimals(eq, actions, terms))
            transformations.update(new_set_decimal)

            
            
            # Filter out invalid equations
            # Remove non equivalent equations
            valid_transformations = {
                (eq, tuple(actions), tuple(terms))
                for eq, actions, terms in transformations
                if self._is_valid_equation(equation, eq)
            }
            
            return transformations
            
        except Exception as e:
            self.logger.error(f"Generation failed: {str(e)}")
            return set()
    
    def _get_terms(self, equation: Eq) -> List[sp.Expr]:
        """Extract terms from both sides with proper sign handling"""
        try:
            lhs_terms = Add.make_args(equation.lhs)
            rhs_terms =  Add.make_args(equation.rhs)
            return list(lhs_terms) + list(rhs_terms)
        except Exception as e:
            self.logger.error(f"Term extraction failed: {str(e)}")
            return []
    
    def _move_term(self, equation: Eq, expr: sp.Expr, terms: List[sp.Expr], 
               actions: List[str] = None, moved_terms: List[sp.Expr] = None) -> Set[Tuple[Eq, Tuple[str, ...], Tuple[sp.Expr, ...]]]:  
        """Move a term between sides with transformation tracking"""
        results = set()

        # Initialize actions and moved_terms lists if not provided
        if actions is None:
            actions = []
        if moved_terms is None:
            moved_terms = []

        # Remove term from terms (creating a new list)
        new_terms = [t for t in terms if t != expr]

        try:
            # Create both simplified and unsimplified forms
            str_lhs = str(equation.lhs)
            str_rhs = str(equation.rhs)
            str_expr = str(expr)

            # Append action and moved term
            new_actions = actions + [f"move term '{expr}'"]
            new_moved_terms = moved_terms + [expr]

            # Convert to tuples for hashing
            new_actions_tuple = tuple(new_actions)
            new_moved_terms_tuple = tuple(new_moved_terms)

            new_lhs = sp.parse_expr(f"{str_lhs} - {str_expr}", evaluate=False)
            new_rhs = sp.parse_expr(f"{str_rhs} - {str_expr}", evaluate=False)

            lhs_list = [xpr.expand(mul=True, evaluate=False) for xpr in sp.Add.make_args(new_lhs)] 
            rhs_list = [xpr.expand(mul=True, evaluate=False) for xpr in sp.Add.make_args(new_rhs)] 
            seme_simpl_lhs = sp.Add(*lhs_list, evaluate=False)
            seme_simpl_rhs = sp.Add(*rhs_list, evaluate=False)
            # Generate different forms
            unsimplified = Eq(new_lhs, new_rhs, evaluate=False)  # Unsimplified form
            simplified = Eq(sp.simplify(new_lhs), sp.simplify(new_rhs))  # Simplified form
            sem_simplified = Eq(sp.simplify(new_lhs), new_rhs) if expr in Add.make_args(equation.lhs) else Eq(new_lhs, sp.simplify(new_rhs))
            sem_simplified_2 = Eq(seme_simpl_lhs, seme_simpl_rhs, evaluate=False)
            # **Check if equation already exists in results before adding**
            for eq in [unsimplified, simplified, sem_simplified, sem_simplified_2]:
                if not any(existing_eq == eq for existing_eq, _, _ in results):
                    results.add((eq, new_actions_tuple, new_moved_terms_tuple))

            # Recursively call _move_term for remaining terms
            if new_terms:
                new_expr = new_terms[0]  # First term in new_terms
                results.update(self._move_term(
                    unsimplified, new_expr, new_terms, new_actions, new_moved_terms
                ))

        except Exception as e:
            self.logger.error(f"Term move failed: {str(e)}")

        return results

    def _divide_factor(self, equation: Eq, variable: Symbol, actions: List[str] = None,  moved_terms: List[sp.Expr] = None) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Divide equation by variable's coefficient with tracking"""
        results = set()
        try:
            coeff = equation.lhs.coeff(variable)
            
            if coeff == 0:
                return results
                
            # Create divided forms
            new_lhs = equation.lhs / coeff
            new_rhs = sp.parse_expr(f"({equation.rhs}) / ({coeff})", evaluate=False)

            rhs_list = [(xpr/coeff).expand(mul=True, evaluate=False) for xpr in sp.Add.make_args(equation.rhs)] 
            seme_simpl_lhs = sp.Add(*rhs_list, evaluate=False)

            # Ensure actions and moved_terms are initialized
            if actions is None:
                actions = []
            if moved_terms is None:
                moved_terms = []

            # Append action and moved term
            new_actions = actions + (f"divide_factor '{equation.lhs}'",)
            new_moved_terms = moved_terms + (coeff,)

            # Convert to tuples for hashing
            new_actions_tuple = tuple(new_actions)
            new_moved_terms_tuple = tuple(new_moved_terms)
            
            seme_simplified = Eq(new_lhs, seme_simpl_lhs)
            # Unsimpified form
            unsimplified = Eq(new_lhs, new_rhs, evaluate=False)
            simplified = Eq(new_lhs, equation.rhs / coeff)
            for eq in [unsimplified, simplified, seme_simplified]:
                results.add((eq, new_actions_tuple, new_moved_terms_tuple ))
            
          

            
            
        except Exception as e:
            self.logger.error(f"Division failed: {str(e)}")
        return results    
    
    def _solve_for_variable(self, equation: Eq, variable: Symbol) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Solve equation for variable with transformation tracking"""
        results = set()
        try:
            
            for sol in sp.solve(equation, variable):
                results.add((Eq(variable, sol), ("solve_variable",), (variable,)))
                
        except Exception as e:
            self.logger.error(f"Solving failed: {str(e)}")
        return results
    
    def _is_valid_equation(self, equation: Eq, transformation:Eq) -> bool:
        """Validate transformed equation"""
        try:
            # check if the transformed equation is not equivilant to the original equation
            # ??????

            return not equation.is_Identity and \
                   equation.lhs != equation.rhs and \
                   not any(arg.is_Boolean for arg in equation.args)
        except Exception:
            return False

    def _convert_to_decimals(self, equation: Eq, actions: List[str] = None,  terms: List[sp.Expr] = None) -> Tuple[Eq, List[str], List[sp.Expr]]:
        """Create decimal representations of equations with exactly two decimals."""
        decimal_equations = set()
        try:
                def format_term(term):
                    """Format numbers in the term but keep symbols unchanged."""
                    if term.is_Number:
                        return int(term) if term == int(term) else round(term, 2)
                    return term  # Keep symbols unchanged
                
                lhs_rounded = sp.simplify(equation.lhs.replace(lambda expr: expr.is_Number, format_term))
                rhs_rounded = sp.simplify(equation.rhs.replace(lambda expr: expr.is_Number, format_term))
                new_eq = Eq(sp.sympify(lhs_rounded), sp.sympify(rhs_rounded))

                # Append action and moved term
                new_actions = actions + (f"convert_to_decimals ",)

                # Convert to tuples for hashing
                new_actions_tuple = tuple(new_actions)
                decimal_equations.add((new_eq, new_actions_tuple ,terms))
                return (new_eq, new_actions, terms)
        except Exception as e:
                # If conversion fails, simply skip this equation.
                self.logger.error(f"Converting to decimal: {str(e)}")
        return  (new_eq, new_actions, terms)

    def _generate_term_combinations(self, expr: sp.Expr) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Generate term combinations without combining like terms."""

        rearrangements = set()
        terms = sp.Add.make_args(expr)
        for k in range(len(terms) + 1):
            for combo in combinations(terms, k):
                # Split into selected (combo) and remaining terms
                remaining = [t for t in terms if t not in combo]
                # Construct left and right without evaluation
                left = sp.Add(*combo, evaluate=False)
                right = sp.Add(*[-term for term in remaining], evaluate=False)
                canonical = simplify(left - right)
                # check if the canonical form equals the original expression
                if canonical - simplify(expr) == 0:
                    term_tuple = tuple(combo)


                # Add all equation permutations
                    rearrangements.add((Eq(left, right), ("term_combination",), term_tuple))
                    rearrangements.add((Eq(right, left), ("term_combination",), term_tuple))
                    rearrangements.add((Eq(-left, -right), ("term_combination",), term_tuple))
                    rearrangements.add((Eq(-right, -left), ("term_combination",), term_tuple))

                    rearrangements.add((Eq(left, right, evaluate=False), ("term_combination",), term_tuple))
                    rearrangements.add((Eq(right, left, evaluate=False), ("term_combination",), term_tuple))
                    rearrangements.add((Eq(-left, -right, evaluate=False), ("term_combination",), term_tuple))
                    rearrangements.add((Eq(-right, -left, evaluate=False), ("term_combination",), term_tuple))
        return rearrangements
    

    def divide_equation_by_term(self, eq: Eq) -> Eq:
        """Safe division with zero-check"""
        try:
            variable = next(iter(eq.lhs.free_symbols), None)
            if variable is None:
                return eq, eq
                
            coff = eq.lhs.coeff(variable)
            if coff == 0:
                raise ValueError("Zero coefficient")
            if coff == 1:
                return eq, eq
            eq_divided = Eq(eq.lhs/coff, eq.rhs/coff)
            lhs_str = str(eq.lhs/coff)
            rhs_str = f"({eq.rhs})/({coff})"
            new_eq = sp.Eq(sp.sympify(lhs_str), sp.parse_expr(rhs_str, evaluate=False), evaluate=False)
   
            return eq_divided, new_eq
        except Exception as e:
            self.logger.debug(f"Division skipped: {str(e)}")
            return eq, eq

    def _divid_singel_term_by_coff(self, triples: Set[Tuple[Eq, List[str], List[sp.Expr]]]) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Robust transformation sequence with error isolation"""
        transformations = set()
        # extract the equations from the triples
        equations = {eq for eq, _, _ in triples}

        try:
            for eq, actions, terms in triples:
                if len(eq.lhs.as_ordered_terms()) == 1 and not eq.lhs.as_ordered_terms()[0].is_integer:
                # check if the term is an integer then skip the division
                    
                    eq1_div, new_eq = self.divide_equation_by_term(eq)
                    simnew_eq = Eq(simplify(new_eq.lhs), simplify(new_eq.rhs))

                    simp_eq1_div = simplify(eq1_div)
                    flip_eq1_div = Eq(eq1_div.rhs, eq1_div.lhs)
                    simp_flip_eq1_div = simplify(flip_eq1_div)


                    if eq1_div not in equations:
                        transformations.add((eq1_div, (actions ,"Divide by coefficient",), (terms , eq.lhs,)))
                    if flip_eq1_div not in equations:
                        transformations.add((flip_eq1_div, (actions ,"Divide by coefficient",), (terms , eq.lhs,)))
                    if simp_eq1_div not in equations:
                        transformations.add((simp_eq1_div, (actions ,"Divide by coefficient",), (terms , eq.lhs,)))
                    if simp_flip_eq1_div not in equations:
                        transformations.add((simp_flip_eq1_div, (actions ,"Divide by coefficient",), (terms , eq.lhs,)))
                    if new_eq not in equations:
                        transformations.add((new_eq, (actions , "Divide by coefficient"),(terms, new_eq)))
                        transformations.add((simnew_eq, (actions , "Divide by coefficient"),(terms, new_eq)))
                 
            return transformations
      


        except Exception as e:
            self.logger.error(f"divid_singel_term_by_coff: {str(e)}")
            return transformations

    def _convert_to_decimals_2(self, triples: Set[Tuple[Eq, List[str], List[sp.Expr]]]) -> Set[Eq]:
        """Create decimal representations of equations with exactly two decimals."""
        decimal_equations = set()
        eqs = {eq for eq, _, _ in triples}
        
        for eq, actions, terms in triples:
            try:
                simplified_eq = Eq(sp.simplify(eq.lhs) ,sp.simplify(eq.rhs))
                lhs_simplified_dec = sp.N(simplified_eq.lhs,2) 
                rhs_simplified_dec = sp.N(simplified_eq.rhs,2) 
                new_simplified_dec = Eq(lhs_simplified_dec, rhs_simplified_dec)
                decimal_equations.add((simplified_eq, actions + ("Simplify",), terms ))
                decimal_equations.add((new_simplified_dec, actions + ("Simplify",), terms ))

                lhs_val = sp.N(eq.lhs,2)
                rhs_val = sp.N(eq.rhs,2 )

                # Convert back to a sympy expression (Float)
                new_eq = Eq(sp.sympify(lhs_val), sp.sympify(rhs_val))
                decimal_eq = Eq(lhs_val, rhs_val, evaluate=False)
                if new_eq not in eqs:
                    decimal_equations.add((new_eq, actions + ("convert_to_decimals",), terms))

                if decimal_eq not in eqs:
                    decimal_equations.add((decimal_eq, actions + ("convert_to_decimals",), terms))
            except Exception as e:
                # If conversion fails, simply skip this equation.
                continue
        return decimal_equations
    
    def subtract_term(self, eq: Eq, term: sp.Expr) -> Eq:
        """Subtract a given term from both sides of an equation."""
        try:
            lhs_str =""
            rhs_str = ""
            term_str = str(term)
            first_char_in_term = term_str[0]
            term_sign = '+'
            # remove the first character if it is a negative or positive sign
            if first_char_in_term == '-' or first_char_in_term == '+':
                term_str = term_str[1:]
            
            if first_char_in_term == '-':
                term_sign = '-'


            # initialize the equation string
            eq_str = Eq(eq.lhs, eq.rhs)
            eq_ter_evaluated = Eq(eq.lhs, eq.rhs)
            # check if term is scalar 
            
            if not term.is_Number :
                if term in eq.lhs.as_ordered_terms() :


                    lhs_str = str(eq.lhs - term)
                    rhs_str = str(eq.rhs)
                    if  term_sign == '+':
                        
                        rhs_str = rhs_str.__add__(' - ').__add__(term_str)
                    else:
                        rhs_str = rhs_str.__add__(' + ').__add__(term_str)
                    
                
                else:

                    rhs_str = str(eq.rhs - term)
                    lhs_str = str(eq.lhs)

                    if term_sign == '+':

                        lhs_str = lhs_str.__add__(' - ').__add__(term_str)
                    
                    else:
                        lhs_str = lhs_str.__add__(' + ').__add__(term_str)

                eq_str = Eq(sp.parse_expr(lhs_str, evaluate=False), sp.parse_expr(rhs_str, evaluate=False))
                eq_ter_evaluated = Eq(sp.parse_expr(lhs_str), sp.parse_expr(rhs_str))

            # remove term_str from the lhs_str
            # Add - term_str to the rhs_str
            #print(f"lhs_str: {lhs_str}")
            
            eq_original = Eq(eq.lhs - term, eq.rhs - term)
            # return both the string and the equation
            return eq_original, eq_str, eq_ter_evaluated
        except Exception as e:
            self.logger.debug(f"Subtraction skipped: {str(e)}")
            return eq, eq, eq

    def transformation_sequence(self, eq: Eq, terms: List[sp.Expr], index: int) -> Set[Eq]:
        """Robust transformation sequence with error isolation"""
        transformations = set()
        try:
            if index >= len(terms):
                return transformations
                
            term = terms[index]
            eq_original, eq1, eq2 = self.subtract_term(eq, term)
            transformations.add((eq1, ("Sbutract term",), (term,)))
            transformations.add((eq2, ("Sbutract term",), (term,)))
            transformations.add((eq_original, ("Sbutract term",), (term,)))
            
            # Add simplified version
            try:
                transformations.add((simplify(eq1), ("Sbutract term and simplify the result",), (term,)))
            except Exception as e:
                self.logger.debug(f"Simplification failed: {str(e)}")

            # Attempt division
            try:
                if len(eq1.lhs.as_ordered_terms()) == 1:
                # check if the term is an integer then skip the division
                    if not eq1.lhs.as_ordered_terms()[0].is_integer:
                        eq1_div , new_eq= self.divide_equation_by_term(eq1)
                        #print(f"eq1_div: {eq1_div}")
                        eq1_div_simpified = simplify(eq1_div)
                        #print(f"eq1_div_simpified: {eq1_div_simpified}")

                        transformations.add((eq1_div, ("Sbutract term", "Divide by coefficient",), (term, eq1.lhs,)))   
                        transformations.add((eq1_div_simpified, ("Sbutract term", "Divide by coefficient",), (term,eq1.lhs,))) 
                    

                # Now subtract each of the other terms from eq1.
                for j, _term in enumerate(terms):
                    if j == index:
                        continue
                    eq_sub1, eq_sub2, eq_3 = self.subtract_term(eq1, _term)
                    #print(f"eq_sub: {eq_sub1}")
                    #print(f"eq_sub2: {eq_sub2}")
                    #print(f"simplify eq_sub: {simplify(eq_sub)}")
                    transformations.add((eq_sub1, ("Sbutract term","Sbutract term",), (term, _term,)))
                    transformations.add((eq_sub2, ("Sbutract term","Sbutract term",), (term, _term,)))
                    transformations.add((simplify(eq_sub1), ("Sbutract term","Sbutract term",), (term, _term,)))
                    transformations.add((eq_3, ("Sbutract term","Sbutract term",), (term, _term,)))
                    # Check if there is only one term on the left side
                    if len(eq_sub1.lhs.as_ordered_terms()) == 1 and not eq_sub1.lhs.as_ordered_terms()[0].is_integer:
                        eq1_div, new_eq = self.divide_equation_by_term(eq_sub1)
                        transformations.add((eq1_div, ("Sbutract term","Sbutract term", "Divide by coefficient",), (term, _term, eq1.lhs,))) 
                        transformations.add((simplify(eq1_div), ("Sbutract term","Sbutract term", "Divide by coefficient",), (term, _term, eq1.lhs,))) 
                        transformations.add((new_eq, ("Sbutract term","Sbutract term", "Divide by coefficient",), (term, _term, eq1.lhs,)))
                        #print(f"eq1_div: {eq1_div}")

                    if len(eq_sub2.lhs.as_ordered_terms()) == 1 and not eq_sub2.lhs.as_ordered_terms()[0].is_integer:
                        eq2_div, new_eq = self.divide_equation_by_term(eq_sub2)
                        transformations.add((eq2_div, ("Sbutract term","Sbutract term", "Divide by coefficient",), (term, _term, eq1.lhs,))) 
                        transformations.add((simplify(eq2_div), ("Sbutract term","Sbutract term", "Divide by coefficient",), (term, _term, eq1.lhs,))) 
                        transformations.add((new_eq, ("Sbutract term","Sbutract term", "Divide by coefficient",), (term, _term, eq1.lhs,)))

                        #print(f"eq2_div: {eq2_div}")
                return transformations


            except Exception as e:
                self.logger.warning(f"Transformation sequence aborted: {str(e)}")
        except Exception as e:
            self.logger.error(f"Transformation sequence failed: {str(e)}")
            return transformations
            
        return transformations

    def generate_all_transformations(self, eq: Eq) -> Set[Eq]:
        """Fault-tolerant transformation generator"""
        transformations = set()
        try:
            terms = self.get_terms(eq)
            for i in range(len(terms)):
                try:
                    transforms = self.transformation_sequence(eq, terms, i)
                    transformations.update(transforms)
                except Exception as e:
                    self.logger.debug(f"Skipping term {i}: {str(e)}")
                    continue
        except Exception as e:
            self.logger.error(f"Transformation generation failed: {str(e)}")
            
        return transformations

    def get_terms(self, eq: Eq) -> List[sp.Expr]:
        """Safe term extraction"""
        try:
            return Add.make_args(eq.lhs) + Add.make_args(eq.rhs)
        except Exception as e:
            self.logger.error(f"Term extraction failed: {str(e)}")
            return []
    
class SimilarityStrategy(ABC):
    """Abstract base class for similarity strategies"""
    @abstractmethod
    def calculate(self, eq1: Eq, eq2: Eq) -> float:
        pass

class StructuralSimilarity(SimilarityStrategy):
    """Structural similarity using sequence matching"""
    @lru_cache(maxsize=1024)
    def _structural_ratio(self, a: str, b: str) -> float:
        return SequenceMatcher(None, a, b).ratio()

    def calculate(self, eq1: Eq, eq2: Eq) -> float:
        try:
            canon1 = simplify(eq1.lhs - eq1.rhs)
            canon2 = simplify(eq2.lhs - eq2.rhs)
            canon3 = simplify(-eq2.lhs + eq2.rhs)
             # Handle cases where simplify() returns True or False
            if isinstance(canon1, bool):
                canon1 = sp.S.Zero
            if isinstance(canon2, bool):
                canon2 = sp.S.Zero
            if isinstance(canon3, bool):
                canon3 = sp.S.Zero

            base_score1 = np.mean([
                self._structural_ratio(sp.srepr(eq1.lhs), sp.srepr(eq2.lhs)),
                self._structural_ratio(sp.srepr(eq1.rhs), sp.srepr(eq2.rhs)),
                self._structural_ratio(sp.srepr(canon1), sp.srepr(canon2))
            ])

            base_score11 = np.mean([
                self._structural_ratio(sp.srepr(eq1.lhs), sp.srepr(eq2.lhs)),
                self._structural_ratio(sp.srepr(eq1.rhs), sp.srepr(eq2.rhs)),
                self._structural_ratio(sp.srepr(canon1), sp.srepr(canon3))
            ])

            base_score2 = np.mean([
                self._structural_ratio(sp.srepr(eq1.rhs), sp.srepr(eq2.lhs)),
                self._structural_ratio(sp.srepr(eq1.lhs), sp.srepr(eq2.rhs)),
                self._structural_ratio(sp.srepr(canon1), sp.srepr(canon2))
            ])

            base_score22 = np.mean([
                self._structural_ratio(sp.srepr(eq1.rhs), sp.srepr(eq2.lhs)),
                self._structural_ratio(sp.srepr(eq1.lhs), sp.srepr(eq2.rhs)),
                self._structural_ratio(sp.srepr(canon1), sp.srepr(canon3))
            ])
            
        
        
            return max(base_score1, base_score2, base_score11, base_score22)
        except Exception as e:
            logging.error(f"Structural similarity failed: {str(e)}")
            return 0

class JaccardSimilarity(SimilarityStrategy):
    """Jaccard similarity based on token overlap"""
    def _tokenize(self, expr: sp.Expr) -> Set[str]:
        return set(re.findall(r'[a-zA-Z]+|\d+|[\+\-\*/\^()]', sp.srepr(expr)))
    
    def calculate(self, eq1: Eq, eq2: Eq) -> float:
        try:
            tokens1 = self._tokenize(eq1.lhs) | self._tokenize(eq1.rhs)
            tokens2 = self._tokenize(eq2.lhs) | self._tokenize(eq2.rhs)
            # Prevent boolean issues
            if isinstance(tokens1, bool) or isinstance(tokens2, bool):
                return 0
            
            intersection = len(tokens1 & tokens2)
            union = len(tokens1 | tokens2)
            return intersection / union if union else 0
        except Exception as e:
            logging.error(f"Jaccard similarity failed: {str(e)}")
            return 0

class CosineSimilarity(SimilarityStrategy):
    """Cosine similarity using token frequency vectors"""
    def _tokenize(self, expr: sp.Expr) -> List[str]:
        """ Tokenize by extracting full terms instead of just individual tokens """
        terms = expr.as_ordered_terms()  # Extract ordered terms like ['4*x', '-7/x']
        return [str(term) for term in terms]  # Convert terms to string format

    def calculate(self, eq1: Eq, eq2: Eq) -> float:
        try:
            # Normalize expressions to canonical form
            expr1 = sp.simplify(eq1.lhs - eq1.rhs)
            expr2 = sp.simplify(eq2.lhs - eq2.rhs)

            # Tokenize the full terms
            tokens1 = self._tokenize(expr1)
            tokens2 = self._tokenize(expr2)

            # Prevent boolean issues
            if isinstance(tokens1, bool) or isinstance(tokens2, bool):
                return 0
            all_tokens = list(set(tokens1 + tokens2))
            vec1 = np.array([tokens1.count(t) for t in all_tokens])
            vec2 = np.array([tokens2.count(t) for t in all_tokens])
            sim = cosine_similarity([vec1], [vec2])[0][0]

            return sim
        except Exception as e:
            logging.error(f"Cosine similarity failed: {str(e)}")
            return 0


#
class EquationComparator:
    """Main comparison handler with strategy pattern"""
    def __init__(self):
        self.preprocessor = EquationPreprocessor()
        self.parser = EquationParser(self.preprocessor)
        self.transformer = TransformationsGenerator()
        self.strategies = {
            'structural': StructuralSimilarity(),
        #    'jaccard': JaccardSimilarity(),
        #    'cosine': CosineSimilarity()
        }

    def compare(self, eq1_str: str, eq2_str: str) -> Dict:
        """Compare two equations using all strategies"""
        try:
            eq1_eval,  eq1_raw  = self.parser.parse(eq1_str)
            eq2_eval,  eq2_raw  = self.parser.parse(eq2_str)
        except Exception as e:
            logging.error(f"Equation parsing failed: {str(e)}")
            return {}

        transformations = self.transformer.generate(eq1_eval)
        transformations.update(self.transformer.generate(eq1_raw))

        equations = {eq for eq, _, _ in transformations}

    #    if eq2_raw in equations:
    #        equations.remove(eq2_raw)

        try:

            results = {
                'total_transformations': len(equations),
            #    'exact_matches': self._find_exact_matches(equations, eq2_raw),
            #    'similarities': self._calculate_similarities(transformations, eq2_raw, eq1_eval),
            #    'similarity_parsed_target': self._calculate_similarities(transformations, eq2_eval, eq1_eval),
                'transformations':  equations
                
            }
            results.update(self._calculate_similarities(transformations, eq2_raw, eq1_eval)['structural'])
        except Exception as e:
            logging.error(f"Comparison failed: {str(e)}")
            return {'exact_matches': [], 'similarities': {}, 'transformations': [], 'total_transformations': 0}

        return results

    def _find_exact_matches(self, transformations: Set[Eq], target: Eq) -> List[str]:
        return [
            str(eq) for eq in transformations
            if simplify(eq.lhs + eq.rhs - (target.lhs + target.rhs)) == 0 # Check the difference between all terms.
        ]

    def _calculate_similarities(self, transformations:Set[Tuple[Eq, List[str], List[sp.Expr]]], target: Eq, original_eq: Eq) -> Dict:
        try:
            results = {}
            target_expr = simplify(target.lhs - target.rhs)
            equations = {eq for eq, _, _ in transformations}

            
            for strategy_name, strategy in self.strategies.items():
                best_match, beset_score = max(
                        ((eq, strategy.calculate(eq, target)) for eq in equations),
                        key = lambda x:x[1],
                        default=(None, -1)
                        )

                actions =[]
                terms = []
                for eq, _actions, _terms in transformations:
                    if simplify(eq.lhs - eq.rhs - (best_match.lhs - best_match.rhs)) == 0:
                        actions = _actions
                        terms  = _terms
                        break

                if best_match is None:
                    results[strategy_name] = {'error':'No valid compersion'}
                    continue

                # calculate symbolic difference
                
                diff1 = simplify((best_match.lhs - best_match.rhs) - (target.lhs - target.rhs))
                diff2 = simplify((best_match.lhs - best_match.rhs) - (target.rhs - target.lhs))
                diff = diff2
                # diff is the less terms have between diff1 and diff2
                diff  = diff1 if len(diff1.args) < len(diff2.args) else diff2


                # Analyze difference
                # finde 
                result = {
                    'best_match': str(f"{best_match.lhs} = {best_match.rhs}"),
                    'max_score': beset_score,
                    'symbolic_difference': str(diff),
                    'variable_factor': [], 
                    'typo': None,
                    'suggested_errors': [],
                    'actions':actions,
                    'terms':terms

                }
                # check if neumeric diff
                if diff ==0:
                    result['suggested_errors'].append( f"Correct")


                elif diff.is_constant():
                    result['suggested_errors'].append( f"Arithmetic error: Constant difference of {diff}")

                # Check for variable factor
                variables_in_diff = diff.free_symbols
                variables_in_target = target_expr.free_symbols
                variables_in_original = original_eq.free_symbols
                variables_in_best_match = best_match.free_symbols
                for var in variables_in_diff:
                    if var not in variables_in_original:
                        result['typo'] = {'variable': str(var)}
                        result['suggested_errors'].append( f"There is no such a vriable  {var} in the original equatoin.")

                    if var not in variables_in_target:
                        result['variable_factor'].append(f"variable: {str(var)}")
                        result['suggested_errors'].append(f"Sign error: where did {var} disappear?")
                        
                    if var not in variables_in_best_match:
                        result['variable_factor'].append(f"variable: {str(var)}")
                        result['suggested_errors'].append( f"Sign error: factor of {var} sould be zero?" )
                                            
                    else:
                        result['variable_factor'].append(f"variable: {str(var)}")
                        result['suggested_errors'].append( f"There is an error calculatin the factor of {var}" )

                    
                results[strategy_name] = result
                

             
        except Exception as e:
            logging.error(f"Similarity calculation failed: {str(e)}")
            results[strategy_name] = {'error': str(e)}
        return results
    

def print_similarity_results(similarities, similarities_parsed, transformations, total_transformations):


    print("\n=== Similarity Analysis ===\n")
    
    for method, details in similarities.items():
        parsed_details = similarities_parsed.get(method, {})  # Get parsed details for the same method
        file.write(f"\n🔍 Similarity Method: {method.capitalize()}\n")
        file.write(f"   ✅ Best Matching Equation: {details['best_match']}\n")
        file.write(f"   ✅ Best Matching to simplified target: {parsed_details['best_match']}\n")

        file.write(f"   📊 Max Similarity Score: {details['max_score']:.4f}\n")
        file.write(f"   📊 Max Similarity Score simplified: {parsed_details['max_score']:.4f}\n")
        file.write(f"   🔀 Symbolic Difference: {details['symbolic_difference']}\n")
        file.write(f"   Actions: {details['actions']}\n")
        file.write(f"   Terms: {details['terms']}\n")
            
        print(f"\n🔍 Similarity Method: {method.capitalize()}")
        print(f"   ✅ Best Matching Equation: {details['best_match']}")
        print(f"   📊 Max Similarity Score: {details['max_score']:.4f}")
        print(f"   🔀 Symbolic Difference: {details['symbolic_difference']}")

        
        if details['neumeric_differnce'] is not None:
            print(f"   🔢 Numeric Difference: {details['neumeric_differnce']}")
            file.write(f"   🔢 Numeric Difference: {details['neumeric_differnce']}\n")
        
        if 'variable_factor' in details and details['variable_factor']:
            print(f"   🔎 Variable Factor Error: {details['variable_factor'].get('variable')}")
            file.write(f"   🔎 Variable Factor Error: {details['variable_factor'].get('variable')}\n")
        
        if 'typo' in details and details['typo']:
            print(f"   ✍️ Possible Typo: {details['typo'].get('variable')}")
            file.write(f"   ✍️ Possible Typo: {details['typo'].get('variable')}\n")
        
        print("   ❌ Suggested Errors:")
        file.write("   ❌ Suggested Errors:\n")
        for error in details.get('suggested_errors', []):
            file.write(f"      - {error}\n")
            print(f"      - {error}")
    
    print("\n=== Rearranged Equations ===\n")
    file.write("\n=== Rearranged Equations ===\n")
    for i, eq in enumerate(transformations, start=1):
        file.write(f"   🔄 {eq.lhs} = {eq.rhs} \n")
        print(f"{i}. {eq}")

    print(f"\nTotal Rearrangements: {total_transformations}\n")
    file.write(f"\n📌 Total Rearrangements: {total_transformations}\n")
    file.write("--------------------------------------\n")
    file.write("--------------------------------------\n")
    file.write("--------------------------------------\n")

# Example usage
if __name__ == "__main__":
    try:
        comparator = EquationComparator()
        df_data = pd.read_json('0_9999_v7.json')
        i = 1
        linear_eq = LinearEquation()  # Linearity checker
        expression_pairs = [
            ("x*(1-a)","-7","x","-a-7/1")
#    ("-4*s + 2", "4*s - r", "2", "-r"),
#  ("-2 - 3", "b + 3", "-3*b", "6"),
#  ("-8*s", "-r - 2", "-8*s", "-2*r"),
#   ("-4*c", "-w/2 - 2", "c", "w/2 - 1/2"),
#   ("-4*n - 3", "3*n - y", "-3", "-1*n - y"),
#   ("3*c + 1", "-c + m", "c", "(m + 1)/4"),
#   ("3*c + 1", "-c + m", "4 + 1", "m"),
#   ("4", "3*d + 4*v", "o", "3*d - 4 + 4*v"),
#   ("o", "3*d + 4*v - 4", "-4*v", "3*d - 4"),
#   #("d*m - d*x - x", "m", "m", "((d + 1)*x)/d - 1"),
#   ("r - 4", "12", "r - 4", "-sqrt(12)"),
#   ("x + 4", "-x - 5*o", "0", "4 + 5*o"),
#   ("8*c", "5 - p", "c", "1.6 - p"),
#   ("-5*c + 5", "3*c + p", "-8*c", "-5*p"),
#   ("b - 3/2", "3/2", "b - 3/2", "-sqrt(3/2)"),
#   ("2*p", "-5 - q", "p", "-2"),
#   ("2*p", "-9*r - 5", "p", "-9/2 - 5/2"),
#   ("-8*s", "-r - 2", "4*s", "r/2"),
#    ("r*(2 - 3)", "0", "r", "0"),
#    ("3*r + 2", "2*x", "1.5*r", "x"),
#    ("2*m + 5", "0", "3*m + 5", "0"),
#    ("2*x", "1", "x", "0"),
#    ("-4*o - 4*y", "2*y + 4", "0", "-1"),
#    ("-f - 1", "0", "-5*f - 4", "0"),
#    ("x*(1 + 2*a)", "-4", "x", "-4/(1 + 2*a)"),
#      ("-4*x", "7 + a*x", "x", "-7/4 + a*x"),
#      ("-9*t - 2", "5*o", "(-9*t)/5 (- 2/5)", "o"),
#       ("x", "3/2*x", "1", "3/2")
    ]

        for index, row in df_data[df_data['category'] == 2].iterrows():
     #   for exp1_lhs, exp1_rhs, exp2_lhs, exp2_rhs in expression_pairs:
            try:
                exp1 = row['t0']
                exp2 = row['t1']
                eq1_str = parse_latex(exp1)
                eq2_str = parse_latex(exp2)
                eq1_str = f"{eq1_str.lhs}={eq1_str.rhs}"
                eq2_str = f"{eq2_str.lhs}={eq2_str.rhs}"
        #        eq1_str = f"{exp1_lhs}={exp1_rhs}"
        #        eq2_str = f"{exp2_lhs}={exp2_rhs}"

                eq1_eval, eq1_raw = comparator.parser.parse(eq1_str)
                eq2_eval, eq2_raw = comparator.parser.parse(eq2_str)

                eq1_eval = eq1_raw.lhs - eq1_raw.rhs
                if not linear_eq.is_linear(eq1_eval):
                    print(f"Skipping: Non-linear equation - {eq1_str}")
                    continue

                eq1_symbols = eq1_eval.free_symbols
                eq2_symbols = eq2_eval.free_symbols
                typo = False
                simis_list = {'o': '0', 'p': 'P', 'q': 'g', 's': 'S', 'x': 'X', 'y': 'Y', '0': 'o', 'P': 'p', 'g': 'q', 'S': 's', 'X': 'x', 'Y': 'y', 'L': 'l', 'l': 'L'}

                if not eq2_symbols.issubset(eq1_symbols):
                    
                    print(f"Typo detected in {eq1_str} or {eq2_str}")
                    new_symbols = eq2_symbols - eq1_symbols

                    for var in new_symbols:
                        var = str(var)
                        if var in simis_list:
                            
                            if simis_list[var] in eq1_str:
                                typo = True
                                print(f"Replacing {var} by {simis_list[var]}")
                                eq2_str = eq2_str.replace(var, simis_list[var])

                print(f"{i} - Comparing:\n{eq1_str} \n{eq2_str}")
                i += 1

                result = comparator.compare(eq1_str, eq2_str)
                if typo:
                    result['typo'] = f"{result.get('typo', '')}?  ,Typo: replacing {var} by {simis_list[var]}"

                result['expression 1'] = eq1_str
                result['expression 2'] = eq2_str

                csv_file = "similarity_results_full_Set.csv"
                fieldnames = list(result.keys())
                #print(result)

                write_header = not os.path.exists(csv_file)
                with open(csv_file, mode='a', newline='', encoding='utf-8') as file:
                    writer = csv.DictWriter(file, fieldnames=fieldnames)
                    if write_header:
                        writer.writeheader()
                    writer.writerow(result)

            except SyntaxError as e:
                print(f"Skipping due to SyntaxError: {e}")
                continue
            except Exception as e:
                print(f"Skipping due to error: {e}")
                continue

    except Exception as e:
        logging.error(f"Critical Error: {e}")

            


