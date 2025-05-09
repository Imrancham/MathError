from collections import defaultdict
import csv
from itertools import combinations
import json
import os
import pandas as pd
import sympy as sp
from sympy import Eq, SympifyError, factor, simplify, Add, S, Symbol, sympify
import re
from difflib import SequenceMatcher
from functools import lru_cache
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
from abc import ABC, abstractmethod
from typing import List, Tuple, Set, Dict
import logging
from sympy import Mul
from sympy import Abs
from sympy.parsing.latex import parse_latex
from sympy.core.expr import Expr

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
        self._variables = sp.symbols('s r b c w n y m d v o p q x a f t h j k l i e g z u X T U V W A B C D E F G H I J K L M N O P Q R S T U V W X Y Z')
    
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
    

    @lru_cache(maxsize=128)
    def parse(self, eq_str: str) -> Tuple[Eq, Eq]:
        """Parse equation string into evaluated and raw forms"""
        try:
            lhs, rhs = map(str.strip, eq_str.split('=', 1))
        except ValueError:
            raise ValueError("Equation must contain exactly one '='")

        lhs_expr = sp.parse_expr(lhs, evaluate=False)
        rhs_expr = sp.parse_expr(rhs, evaluate=False)

        simplified = Eq(parse_latex(lhs), parse_latex(rhs)) 
        not_simplified = Eq(parse_latex(lhs), parse_latex(rhs), evaluate=False)
        equation = Eq(lhs_expr, rhs_expr)
      #  eq2 = Eq(lhs, rhs)
      #  eq3 = Eq(lhs, rhs, evaluate=False)


        return (
            Eq(simplify(lhs_expr), simplify(rhs_expr)),
            Eq(lhs_expr, rhs_expr,  evaluate=False)
        )


class TransformationsGenerator:
    """Generates equation transformations with transformation tracking"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def generate(self, equation: Eq) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Generate all possible equation transformations with transformation tracking"""
        transformations = set()
        try:


            transformations.add((equation, ("Original",), tuple()))

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
                new_set_decimal.update(self._convert_to_decimals(eq, actions, terms))
            for eq , actions, terms in new_set_decimal:
                # check in eq in eqsuations of the transformations tuble

                if any(eq == t[0] for t in transformations):
                    continue
            
                transformations.add((eq, actions, terms))
            
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

        # Initialize mutable arguments if None
        actions = actions or []
        moved_terms = moved_terms or []

        # Create a new term list without 'expr'
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
            unsimplified = Eq(new_lhs, new_rhs, evaluate=False)  # Unsimplified form

            lhs_simplified = Eq(simplify(new_lhs), new_rhs)
            rhs_simplified = Eq(new_lhs, simplify(new_rhs))

             

            simplified = Eq(sp.parse_expr(f"{str_lhs} - {str_expr}", evaluate=True), sp.parse_expr(f"{str_rhs} - {str_expr}", evaluate=True))  # Simplified form
            sem_simplified = Eq(sp.parse_expr(f"{str_lhs} - {str_expr}", evaluate=True), new_rhs) if expr in Add.make_args(equation.lhs) else Eq(new_lhs, sp.parse_expr(f"{str_rhs} - {str_expr}", evaluate=True))
            
            lhs_list = [xpr.expand(mul=True, evaluate=False) for xpr in sp.Add.make_args(new_lhs)] 
            rhs_list = [xpr.expand(mul=True, evaluate=False) for xpr in sp.Add.make_args(new_rhs)]
            seme_simpl_lhs = sp.Add(*lhs_list, evaluate=False)
            seme_simpl_rhs = sp.Add(*rhs_list, evaluate=False)
            sem_simplified_2 = Eq(seme_simpl_lhs, seme_simpl_rhs, evaluate=False)
            
            
            # **Check if equation already exists in results before adding**
            equations = {unsimplified, simplified, sem_simplified, sem_simplified_2, lhs_simplified, rhs_simplified}
            results.update((eq, new_actions_tuple, new_moved_terms_tuple) for eq in equations if eq not in {r[0] for r in results})

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
            # Ensure actions and moved_terms are initialized
            actions = actions or ()
            moved_terms = moved_terms or ()

            # Append action and moved term
            new_actions_tuple = actions + (f"divide_factor '{equation.lhs}'",)
            new_moved_terms_tuple = moved_terms + (coeff,)
                
            if coeff == 0:
                return results
                
            # Create divided forms
            new_lhs = sp.parse_expr(f"({equation.lhs}) / ({coeff})", evaluate=False)
            new_rhs = sp.parse_expr(f"({equation.rhs}) / ({coeff})", evaluate=False)
            # Unsimpified form
            unsimplified = Eq(new_lhs, new_rhs, evaluate=False)

            lhs_simplified = Eq(simplify(new_lhs), new_rhs, evaluate=False)
            rhs_simplified = Eq(new_lhs, simplify(new_rhs), evaluate=False)

            rhs_list = [(xpr/coeff).expand(mul=True, evaluate=False) for xpr in sp.Add.make_args(equation.rhs)] 
            seme_simpl_lhs = sp.Add(*rhs_list, evaluate=False)
            seme_simplified = Eq(new_lhs, seme_simpl_lhs)
            seme_simplified2 = Eq(simplify(new_lhs), seme_simpl_lhs)


            rhs_simplified2 = Eq(new_lhs, simplify(new_rhs), evaluate=True) 
            simplified = Eq(simplify(new_lhs), simplify(new_rhs), evaluate=True) 
            equations = {unsimplified, simplified, lhs_simplified, rhs_simplified, seme_simplified, rhs_simplified2, seme_simplified2}
            for eq in equations:
                results.add((eq, new_actions_tuple, new_moved_terms_tuple ))
            
            
        except Exception as e:
            self.logger.error(f"Division failed: {str(e)}")
        return results    
    
    def _solve_for_variable(self, equation: Eq, variable: Symbol) -> Set[Tuple[Eq, List[str], List[sp.Expr]]]:
        """Solve equation for variable with transformation tracking"""
        results = set()
        try:
            
            for sol in sp.solve(equation, variable):
                results.add((Eq(variable, sol), (f"solve for variable: {variable}",), (variable,)))
                
        except Exception as e:
            self.logger.error(f"Solving failed: {str(e)}")
        return results
    
    def _is_valid_equation(self, equation: Eq, transformed:Eq) -> bool:
        """Validate transformed equation"""
        try:
            # check if the transformed equation is  equivilant to the original equation
            if transformed.lhs == transformed.rhs:
                return False
            # check if the transformed equation is equivalent to the original equation

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
                lhs_simplified = Eq(sp.simplify(equation.lhs), equation.rhs)
                rhs_simplified = Eq(equation.lhs, sp.simplify(equation.rhs))    
                simplified = Eq(sp.simplify(lhs_rounded), sp.simplify(rhs_rounded))
                eqs = {lhs_simplified, rhs_simplified, simplified,new_eq }
                for eq in eqs:
                
                    decimal_equations.add((eq, new_actions_tuple ,terms))
                return decimal_equations
        except Exception as e:
                # If conversion fails, simply skip this equation.
                self.logger.error(f"Converting to decimal: {str(e)}")
        return  decimal_equations


#   
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
            eq1_lhs = str(eq1.lhs).replace(" ", "")
            eq1_rhs = str(eq1.rhs).replace(" ", "")
            eq2_lhs = str(eq2.lhs).replace(" ", "")
            eq2_rhs = str(eq2.rhs).replace(" ", "")

            if( eq1_lhs == eq2_lhs and eq1_rhs == eq2_rhs) or( eq1_lhs == eq2_rhs and eq1_rhs == eq2_lhs):
                return 1.0
          #  if sp.srepr(eq1) == sp.srepr(eq2):
           #         return 1.0

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



class EquationComparator:
    """Main comparison handler with strategy pattern"""
    def __init__(self):
        #self.preprocessor = EquationPreprocessor()
        self.parser = EquationParser()
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


        try:

            results = {
                'total_transformations': len(equations),
            #    'exact_matches': self._find_exact_matches(equations, eq2_raw),
            #    'similarities': self._calculate_similarities(transformations, eq2_raw, eq1_eval),
            #    'similarity_parsed_target': self._calculate_similarities(transformations, eq2_eval, eq1_eval),
                'transformations':  equations
                
            }
            
            # check if eq2_eval,  eq2_raw as strings are the same
            for strategy_name, strategy in self.strategies.items():
                if str(f"{eq2_eval.lhs}={eq2_eval.rhs}").replace(" ", "") == str(f"{eq2_raw.lhs}={eq2_raw.rhs}").replace(" ", ""):
                    results.update(self._calculate_similarities(transformations, eq2_eval, eq1_raw)[strategy_name])

                else:

                    results.update(self._calculate_similarities(transformations, eq2_raw, eq1_raw)[strategy_name])
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

            
            for strategy_name, strategy in self.strategies.items():
                best_match, best_score, best_actions, terms_ = max(
                    ((eq, strategy.calculate(eq, target), actions, terms_) for eq, actions, terms_ in transformations),
                    key=lambda x: x[1],
                    default=(None, -1, [])
                )

                if best_match is None:
                    results[strategy_name] = {'error': 'No valid comparison'}
                    continue

                # calculate symbolic difference
                def similarity_score(expr1: Expr, expr2: Expr) -> int:
                    """Rough similarity score based on matching variables/structure."""
                    return len(expr1.free_symbols & expr2.free_symbols)
                
                lhs_score = similarity_score(target.lhs, best_match.lhs)
                rhs_score = similarity_score(target.lhs, best_match.rhs)

                # Choose the most similar side
                if lhs_score >= rhs_score:
                    eq2_match = best_match.lhs
                    eq2_other = best_match.rhs
                else:
                    eq2_match = best_match.rhs
                    eq2_other = best_match.lhs


                lhs_diff = simplify(target.lhs - eq2_match)
                rhs_diff = simplify(target.rhs - eq2_other)
                diff = simplify(lhs_diff - rhs_diff)



                # Analyze difference
                # finde 
                Missing_variable, variable_mismatch = self._find_typo_mismach(diff)
                result = {
                    f"{strategy_name}_best_match": str(f"{best_match.lhs} = {best_match.rhs}"),
                    f"{strategy_name}_max_score": best_score,
                    f"{strategy_name}_symbolic_difference": str(diff),
                    f"{strategy_name}_variable_factor": [], 
                    "typo": None,
                    f"{strategy_name}_suggested_errors": [],
                    f"{strategy_name}_actions":best_actions,
                    f"{strategy_name}_terms":terms_
                }
                target_terms = target.lhs.as_ordered_terms() + target.rhs.as_ordered_terms()
                original_terms = original_eq.lhs.as_ordered_terms() + original_eq.rhs.as_ordered_terms()

                tolerance = 1e-3
                if result[f"{strategy_name}_max_score"] == 1.0 or ( diff.is_number and Abs(diff) < tolerance ):
                    if self._is_num_over_linear(target.lhs) or self._is_num_over_linear(target.rhs):
                        result[f"{strategy_name}_suggested_errors"].append( f"Correct")
                        result[f"{strategy_name}_suggested_errors"].append( f"Check assumption. Denominator might be zero.")
                    
                    else:
                        result[f"{strategy_name}_suggested_errors"].append( f"Correct")
                    
                
                        

                else:
                    if diff.is_number and result[f"{strategy_name}_suggested_errors"] == []:
                        half_diff = diff/2
                        matched_term = next((term for term in original_terms if simplify(term - half_diff) == 0), None)
                        if matched_term is not None:
                            result[f"{strategy_name}_suggested_errors"].append(f"Sign error 1: Be aware when adding {matched_term}")

                        else:
                            matched_term = next((term for term in original_terms if simplify(term + half_diff) == 0), None)
                            if matched_term is not None:
                                result[f"{strategy_name}_suggested_errors"].append(f"Sign error 2: Be aware when adding {matched_term}")

                            else:
                                matched_term = next((term for term in target_terms if simplify(term + half_diff) == 0), None)
                                if matched_term is not None:
                                    result[f"{strategy_name}_suggested_errors"].append(f"Sign error 3: Be aware of {matched_term} sign")

                                else:
                                    matched_term = next((term for term in target_terms if simplify(term - half_diff) == 0), None)
                                    if matched_term is not None:
                                        result[f"{strategy_name}_suggested_errors"].append(f"Sign error 4: Be aware of {matched_term} sign")



                                    else:
                                            result[f"{strategy_name}_suggested_errors"].append( f"Arithmetic error")

                    # check if the target_eq and best_mach of the form number*sybmbol = number
                    # if yes check if there is a sign error
                    elif Missing_variable:
                        result['typo'] = {f'Missing variable {Missing_variable[0][1]}'}
                        result[f"{strategy_name}_suggested_errors"].append( f"Typo: A missing variable {Missing_variable[0][1]} of the cofficient {{Missing_variable[0][0]}}.")
                    
                    elif variable_mismatch:
                        result['typo'] = {f'variable mismatch {variable_mismatch[0][1]}'}
                        result[f"{strategy_name}_suggested_errors"].append( f"Typo: A variable mismatch {variable_mismatch[0][1]} of the cofficient {variable_mismatch[0][0]}.")

                    elif self._is_eq_of_form(best_match) and self._is_eq_of_form(target): # form a*x * b
                        # Check for sign error
                        bool_sing, sign_type = self._is_sign_error(best_match, target)
                        if bool_sing and sign_type == 'value':
                            result[f"{strategy_name}_suggested_errors"].append( f"Sign error 5: Be awair of {diff/2} sign")
                        elif bool_sing and sign_type == 'variable':
                            result[f"{strategy_name}_suggested_errors"].append( f"Sign error 6: Be awair of {diff/2} sign when moving terms.")
                    #   else:
                    #      result[f"{strategy_name}_suggested_errors"].append( f"Check your calculatoin")

                    if result[f"{strategy_name}_suggested_errors"] == []:
                        if self._is_num_over_linear(target.lhs) or self._is_num_over_linear(target.rhs):
                          result[f"{strategy_name}_suggested_errors"].append( f"Check assumption. Denominator might be zero.")
                        for term in diff.as_ordered_terms():
                                half_diff = term/2

                                matched_term = next((term for term in original_terms if simplify(term - half_diff) == 0), None)
                                if matched_term is not None:
                                    result[f"{strategy_name}_suggested_errors"].append(f"Sign error 7: Be aware when adding {matched_term}")

                                else:
                                    matched_term = next((term for term in original_terms if simplify(term + half_diff) == 0), None)
                                    if matched_term is not None:
                                        result[f"{strategy_name}_suggested_errors"].append(f"Sign error 8: Be aware when adding {matched_term}")

                                    else:
                                        matched_term = next((term for term in target_terms if simplify(term + half_diff) == 0), None)
                                        if matched_term is not None:
                                            result[f"{strategy_name}_suggested_errors"].append(f"Sign error 9: Be aware of {matched_term} sign")

                                        else:
                                            matched_term = next((term for term in target_terms if simplify(term - half_diff) == 0), None)
                                            if matched_term is not None:
                                                result[f"{strategy_name}_suggested_errors"].append(f"Sign error 10: Be aware of {matched_term} sign")

                            #    else:
                            #        result[f"{strategy_name}_suggested_errors"].append( f"Arithmetic error")

                
                        # Check for variable factor
                        variables_in_diff = diff.free_symbols
                        variables_in_target = target_expr.free_symbols
                        variables_in_original = original_eq.free_symbols
                        variables_in_best_match = best_match.free_symbols
                        for var in variables_in_diff:
                            if var not in variables_in_original:
                                result['typo'] = {'variable': str(var)}
                                result[f"{strategy_name}_suggested_errors"].append( f"Typo: There is no such a vriable  {var} in the original equatoin.")

                            elif var not in variables_in_target:
                                result[f"{strategy_name}_variable_factor"].append(f"variable: {str(var)}")
                                result[f"{strategy_name}_suggested_errors"].append(f"Missing variable: where did {var} disappear?")
                                
                            elif var not in variables_in_best_match:
                                result[f"{strategy_name}_variable_factor"].append(f"variable: {str(var)}")
                                result[f"{strategy_name}_suggested_errors"].append( f"Extra varialbe: factor of {var} sould be zero?" )
                                                    
                        #    else  :
                        #        result[f"{strategy_name}_variable_factor"].append(f" {str(var)}")
                        #        result[f"{strategy_name}_suggested_errors"].append( f"There is an error calculatin the factor of {var}" )

                        
                results[strategy_name] = result
                

             
        except Exception as e:
            logging.error(f"Similarity calculation failed: {str(e)}")
            results[strategy_name] = {'error': str(e)}
        return results
    
    
    def _is_num_over_linear(self, expr: Expr) -> bool:
                num, denom = expr.as_numer_denom()

                # Check numerator is a number
                if not num.is_number:
                    return False

                # Denominator must be an addition of two terms
                if not isinstance(denom, Add):
                    return False

                terms = denom.as_ordered_terms()

                if len(terms) != 2:
                    return False

                has_linear_term = False
                has_constant_term = False

                for term in terms:
                    if term.is_number:
                        has_constant_term = True
                    elif isinstance(term, Mul):
                        coeff, vars_ = term.as_coeff_mul()
                        # Check only one symbol, and it has power 1
                        if len(vars_) == 1 and isinstance(vars_[0], Symbol):
                            has_linear_term = True
                    elif isinstance(term, Symbol):
                        has_linear_term = True
                    else:
                        return False

                return has_linear_term and has_constant_term
    def _is_eq_of_form(self, eq: Eq) -> bool:
        """
        Check if the equation is of the form:
        symbol = number,
        number = symbol,
        number*symbol = number,
        number = number*symbol
        including all equivalent orientations.
        """
        lhs, rhs = eq.lhs, eq.rhs

        def is_symbol_or_number_mul_symbol(expr):
            # Case 1: single variable
            if isinstance(expr, Symbol):
                return True
            # Case 2: multiplication like 2*x or x*2
            if expr.is_Mul:
                has_symbol = any(isinstance(arg, Symbol) for arg in expr.args)
                has_number = any(arg.is_number for arg in expr.args)
                return has_symbol and has_number
            return False

        def is_symbolic_form(expr1, expr2):
            return (
                expr1.is_number and is_symbol_or_number_mul_symbol(expr2)
            ) or (
                expr2.is_number and is_symbol_or_number_mul_symbol(expr1)
            )

        return is_symbolic_form(lhs, rhs)


    def _is_sign_error(self, eq_ref: Eq, eq_cmp: Eq) -> Tuple[bool, str]:
        """
        Check if eq_cmp differs from eq_ref only by a sign error.
        Returns a tuple (True/False, 'variable'/'value'/'both'/None)
        """
        def get_var_side(eq: Eq) -> Tuple[Expr, Expr]:
            # Return (var_expr, value_expr) such that var_expr contains the symbol
            if eq.rhs.is_number:
                return eq.lhs, eq.rhs
            return eq.rhs, eq.lhs

        ref_var, ref_val = get_var_side(eq_ref)
        cmp_var, cmp_val = get_var_side(eq_cmp)

        # Determine which component has the sign error
        var_diff = simplify(ref_var - cmp_var)
        var_sign_flip = simplify(ref_var + cmp_var)

        val_diff = simplify(ref_val - cmp_val)
        val_sign_flip = simplify(ref_val + cmp_val)

        if var_diff == 0 and val_sign_flip == 0:
            return True, 'variable'
        elif val_diff == 0 and var_sign_flip == 0:
            return True, 'value'

        else:
            return False, None

    def _extract_coeff_and_var(self, term):
        if isinstance(term, Mul):
            coeff, vars_ = term.as_coeff_mul()
            return coeff, vars_
        elif term.is_number:
            return term, ()
        else:
            return 1, (term,)

    

    def _find_typo_mismach(self, expr: Expr):
        """Check for missing variables and variable mismatch in the expression."""

        try:
            terms = expr.as_ordered_terms()
        except Exception as e:
            print(f"[Error] Failed to get terms from expression: {e}")
            return [], []

        grouped = defaultdict(list)
        missing_var_details = []  # (coefficient, [vars])
        mismatch_var_details = []  # (coefficient, [[vars], [vars], ...])
        new_expr_parts = []

        for term in terms:
            try:
                coeff, vars_ = self._extract_coeff_and_var(term)
                grouped_key = abs(coeff)
                grouped[grouped_key].append((term, coeff, vars_))
            except Exception as e:
                print(f"[Warning] Failed to extract coefficient/vars from term '{term}': {e}")
                continue

        for coeff_abs, term_data in grouped.items():
            try:
                current_terms = [term for term, _, _ in term_data]
                has_constant = any(not vars_ for _, _, vars_ in term_data)
                vars_list = [tuple(sorted(vars_)) for _, _, vars_ in term_data if vars_]
                unique_vars = set(vars_list)

                if has_constant and len(term_data) >= 2:
                    flat_vars = sorted({str(v) for var_tuple in vars_list for v in var_tuple})
                    missing_var_details.append((coeff_abs, flat_vars))

                if len(unique_vars) > 1:
                    mismatch_groups = [[str(v) for v in tup] for tup in unique_vars]
                    mismatch_var_details.append((coeff_abs, mismatch_groups))

                combined = Add(*current_terms)
                if not has_constant and len(unique_vars) == 1:
                    try:
                        combined = factor(combined)
                    except Exception as e:
                        print(f"[Info] Could not factor terms with coeff {coeff_abs}: {e}")

                new_expr_parts.append(combined)

            except Exception as e:
                print(f"[Error] Processing group with coeff {coeff_abs} failed: {e}")
                continue

        return missing_var_details, mismatch_var_details


if __name__ == "__main__":
    try:
        comparator = EquationComparator()
        df_data = pd.read_json('0_9999_v7.json')
        i = 1
        linear_eq = LinearEquation()  # Linearity checker
        expression_pairs = [
            ("-20*x","64", "x","-3.20000000000000"),
            ("x(5*a - 4)","3","x","3/(5*a - 4)"), 
            ("5*h","7*s - 1*2", "h","(7/5)*h - 2/5"),
            ("-4*x - 3","7","4*x + 3","-7"),
            ("8","-36*x","-8/36","x"),
          #  ("-5*k - 5","-k + 4*m", "-4*k","4*m + 2"), #Check move terms operations.
          #  ("-2*r","2*z + 1","r","-z + 1/2"), # Sould get 'Sign error'
          #  ("-2*d","2*q + 5", "d","-q + 5/2"), # Sould get 'Sign error'
          #  ("5*x + 3","3", "5*x","-3"), # Sould get 'Sign error'
          #  ("x - 3","-3","x","3"), # Sould get 'Sign error'
          #  ("-4*c","3*p - 1","c","-3/4*p + 1/4"),
          #  ("x*(a + 1)","0","x","0/0"),
          #  ("-3*l - 3","d - 5*l","-8*l","4*d"), # Too many errors, need to reduce the similarity by contiputing the symbolic_difference.
          #  ("-2*o","3*y - 5","o","(3/2)*y + 5/2"),  #Check differenace functions
          #  ("5*c + 4*w","3*c + 2","4*w","2*c + 2"), 
          #  ("-3*y","4*v + 3","y","-4/3*v - 7"),
          ("3*z - 2","0", "2*z + 5","0"),
            ("x + 1","25*(a*x) - 5","-5*x + x","-6"), # Check differenace functions
            ("1 - 5*z","2*z - 5","7*z","-6"), # Check differenace functions
            ("4 - 2*x","5*x - 5","7*x","-9"), # Check differenace functions
            ("-6*(x + 6)","7","x + 6","-7"), # Handel cases where '()' are exists
            ("3*x - 4","-2*u","1.5*x + 2","u"), # Sould get 'Sign error'
        
            ("0","x(3*a + 1) + 6","1/x","(3*a + 1) + 6"),
            ("5*x - 4","-3","x","-1/5"),
            ("-o + (5*x + (-5*x - 15))","0","-10*x - 25","0"),
            ("2 - 3*q","2*n + 4*q", "2","2*n + 7*p"),
            ("(2*x+(5*x+5))-6","0","7*x-1","0"),
             ("3-4*x","7","4*x+3","-7"),
            ("8*x+7","1","x","-3/4"),	
            ("-6*x","1","x","-1/6"),
#            ("-1","x(a-2)","x","-1/(a-2)")
#            ("-22*x","-38","x","38/22"),
#            ("4-9*x","6","9*x-4","-6")
#            ("6*x-2","-5","x","-1/3")
#            ("9*x-4","-1","x","6/9")
#            ("x*(1-a)","-7","x","-a-7/1")
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
            
            try:
                exp1 = row['t0']
                exp2 = row['t1']
                eq1_str = parse_latex(exp1)
                eq2_str = parse_latex(exp2)
                eq1_str = f"{eq1_str.lhs}={eq1_str.rhs}"
                eq2_str = f"{eq2_str.lhs}={eq2_str.rhs}"
    #    for exp1_lhs, exp1_rhs, exp2_lhs, exp2_rhs in expression_pairs:
    #        try:
    #            
    #            eq1_str = f"{exp1_lhs}={exp1_rhs}"
    #            eq2_str = f"{exp2_lhs}={exp2_rhs}"

                #eq1_str = parse_latex(eq1_str)
                #eq2_str = parse_latex(eq2_str)

                eq1_eval, eq1_raw = comparator.parser.parse(eq1_str)
                eq2_eval, eq2_raw = comparator.parser.parse(eq2_str)

                eq1_eval = eq1_raw.lhs - eq1_raw.rhs
                if not linear_eq.is_linear(eq1_eval):
                    print(f"Skipping: Non-linear equation - {eq1_str}")
                    continue

                eq1_symbols = eq1_eval.free_symbols
                eq2_symbols = eq2_eval.free_symbols
                typo = False
                simis_list = {'o': '0', 'p': 'P', 'q': 'g', 's': 'S', 'x': 'X', 'y': 'Y', '0': 'o', 'P': 'p', 'g': 'q', 'S': 's', 'X': 'x', 'Y': 'y', 'L': 'l', 'l': 'L', 'w':'W', 'W':'w',  'c':'C', 'C':'c', 'f':'F', 'F':'f'}

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
                print(result)

        #        write_header = not os.path.exists(csv_file)
        #        with open(csv_file, mode='a', newline='', encoding='utf-8') as file:
        #               writer = csv.DictWriter(file, fieldnames=fieldnames)
        #               if write_header:
        #                   writer.writeheader()
        #               writer.writerow(result)
#
                
            except Exception as e:
                print(f"Skipping due to error: {e}")
                continue

    except Exception as e:
        logging.error(f"Critical Error: {e}")

            


