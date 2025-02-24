import pandas as pd
import sympy as sp
from sympy import Eq, simplify, symbols, Add, S
from itertools import combinations
import re
from difflib import SequenceMatcher
from functools import lru_cache
from collections import Counter
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
from abc import ABC, abstractmethod
from typing import List, Tuple, Set, Dict
from sympy.parsing.latex import parse_latex

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
            Eq(lhs_expr, rhs_expr)
        )

    def _parse_side(self, side_str: str) -> sp.Expr:
        """Parse a single side of an equation"""
        return sp.parse_expr(
            self.preprocessor.preprocess(side_str),
            evaluate=False
        )

class RearrangementGenerator:
    """Generates algebraic rearrangements of equations"""
    def __init__(self, max_terms: int = 4):
        self.max_terms = max_terms

    def generate(self, equation: Eq) -> Set[Eq]:
        """Generate all valid rearrangements of an equation"""
        canonical = simplify(equation.lhs - equation.rhs)
        rearrangements = set()

        rearrangements.update(self._generate_term_combinations(canonical))
        rearrangements.update(self._solve_for_variables(canonical))
        rearrangements.update(self._create_canonical_forms(canonical))
        rearrangements.update(self._convert_to_decimals(rearrangements))

        all_transforms = set()
        try:
            terms = get_terms(equation)
            for i in range(len(terms)):
                transforms = transformation_sequence(equation, terms, i)
                all_transforms.update(transforms)
        except Exception as e:
            raise Exception(f"Error in generate_all_transformations: {e}")
        rearrangements.update(all_transforms)

        return rearrangements

    def _generate_term_combinations(self, expr: sp.Expr) -> Set[Eq]:
        """Generate term combinations with limited size"""
        terms = tuple(Add.make_args(expr))
        rearrangements = set()

        for k in range(min(len(terms), self.max_terms) + 1):
            for combo in combinations(terms, k):
                remaining = tuple(t for t in terms if t not in combo)
                left, right = Add(*combo), -Add(*remaining)
                rearrangements.update({
                    Eq(left, right), Eq(right, left),
                    Eq(-left, -right), Eq(-right, -left)
                })

        return rearrangements

    def _solve_for_variables(self, expr: sp.Expr) -> Set[Eq]:
        """Solve equation for each variable"""
        solutions = set()
        for var in expr.free_symbols:
            try:
                solution = sp.solve(expr, var, dict=True)
                if solution:
                    solutions.add(Eq(var, solution[0][var]))
            except NotImplementedError:
                continue
        return solutions

    def _create_canonical_forms(self, expr: sp.Expr) -> Set[Eq]:
        """Create standard canonical forms"""
        return {Eq(expr, 0), Eq(0, expr), Eq(-expr, 0), Eq(0, -expr)}

    def _convert_to_decimals(self, equations: Set[Eq]) -> Set[Eq]:
        """Create decimal representations of equations with exactly two decimals."""
        decimal_equations = set()
        for eq in equations:
            try:
                lhs_val = sp.N(eq.lhs)
                rhs_val = sp.N(eq.rhs)
                # Format the evaluated expressions to two decimals
                lhs_str = f"{lhs_val:.2f}"
                rhs_str = f"{rhs_val:.2f}"
                # Convert back to a sympy expression (Float)
                new_eq = Eq(sp.sympify(lhs_str), sp.sympify(rhs_str))
                decimal_equations.add(new_eq)
            except Exception as e:
                # If conversion fails, simply skip this equation.
                continue
        return decimal_equations


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
        canon1 = simplify(eq1.lhs - eq1.rhs)
        canon2 = simplify(eq2.lhs - eq2.rhs)
        return np.mean([
            self._structural_ratio(sp.srepr(eq1.lhs), sp.srepr(eq2.lhs)),
            self._structural_ratio(sp.srepr(eq1.rhs), sp.srepr(eq2.rhs)),
            self._structural_ratio(sp.srepr(canon1), sp.srepr(canon2))
        ])

class JaccardSimilarity(SimilarityStrategy):
    """Jaccard similarity based on token overlap"""
    def _tokenize(self, expr: sp.Expr) -> Set[str]:
        return set(re.findall(r'[a-zA-Z]+|\d+|[\+\-\*/\^()]', sp.srepr(expr)))

    def calculate(self, eq1: Eq, eq2: Eq) -> float:
        tokens1 = self._tokenize(eq1.lhs) | self._tokenize(eq1.rhs)
        tokens2 = self._tokenize(eq2.lhs) | self._tokenize(eq2.rhs)
        
        intersection = len(tokens1 & tokens2)
        union = len(tokens1 | tokens2)
        return intersection / union if union else 0

class CosineSimilarity(SimilarityStrategy):
    """Cosine similarity using token frequency vectors"""
    def _tokenize(self, expr: sp.Expr) -> List[str]:
        return re.findall(r'[a-zA-Z]+|\d+|[\+\-\*/\^()]', sp.srepr(expr))

    def calculate(self, eq1: Eq, eq2: Eq) -> float:
        tokens1 = self._tokenize(eq1.lhs) + self._tokenize(eq1.rhs)
        tokens2 = self._tokenize(eq2.lhs) + self._tokenize(eq2.rhs)

        all_tokens = list(set(tokens1 + tokens2))
        vec1 = np.array([tokens1.count(t) for t in all_tokens])
        vec2 = np.array([tokens2.count(t) for t in all_tokens])

        return cosine_similarity([vec1], [vec2])[0][0]

class EquationComparator:
    """Main comparison handler with strategy pattern"""
    def __init__(self):
        self.preprocessor = EquationPreprocessor()
        self.parser = EquationParser(self.preprocessor)
        self.rearranger = RearrangementGenerator()
        self.strategies = {
            'structural': StructuralSimilarity(),
            'jaccard': JaccardSimilarity(),
            'cosine': CosineSimilarity()
        }

    def compare(self, eq1_str: str, eq2_str: str) -> Dict:
        """Compare two equations using all strategies"""
        try:
            eq1_eval,  eq1_raw = self.parser.parse(eq1_str)
            eq2_eval,  _       = self.parser.parse(eq2_str)
        except Exception as e:
            return {'error': str(e)}

        rearrangements = self.rearranger.generate(eq1_eval)
        rearrangements.update(self.rearranger.generate(eq1_raw))

        results = {
            'exact_matches': self._find_exact_matches(rearrangements, eq2_eval),
            'similarities': self._calculate_similarities(
                rearrangements, eq2_eval
            ),
            'rearrangements': list(map(str, rearrangements)),
            'total_rearrangements': len(rearrangements)
        }

        return results

    def _find_exact_matches(self, rearrangements: Set[Eq], target: Eq) -> List[str]:
        return [
            str(eq) for eq in rearrangements
            if simplify(eq.lhs - eq.rhs - (target.lhs - target.rhs)) == 0
        ]

    def _calculate_similarities(self, rearrangements: Set[Eq], target: Eq) -> Dict:
        return {
            strategy_name: {
                'best_match': str(max(
                    [(eq, strategy.calculate(eq, target)) for eq in rearrangements],
                    key=lambda x: x[1],
                    default=(Eq(S.Zero, S.Zero), -1)
                )[0]),
                'max_score': max(
                    [(eq, strategy.calculate(eq, target)) for eq in rearrangements],
                    key=lambda x: x[1],
                    default=(None, -1)
                )[1]
            }
            for strategy_name, strategy in self.strategies.items()
        }
    
def subtract_term(eq: Eq, term: sp.Expr) -> Eq:
    """Subtract a given term from both sides of an equation."""
    lhs_str =""
    rhs_str = ""
    term_str = str(term)
    variable = term_str.__getitem__(-1)
    # initialize the equation string
    eq_str = Eq(eq.lhs, eq.rhs)
    
    # check if term is scalar 
    
    if not term.is_Number :
        if term in eq.lhs.as_ordered_terms() :
            # check the sign of the term
            #print("on the left side")
            term_sign = eq.lhs.coeff(variable)
            #print(f"term: {term}")
            #print(f"term_sign: {term_sign}")

            # Subtract the term from the lhs
            lhs_str = str(eq.lhs - term)
            rhs_str = str(eq.rhs)
            if term_sign > 0:
                
                rhs_str = rhs_str.__add__(' - ').__add__(term_str)
            else:
                rhs_str = rhs_str.__add__(' + ').__add__(term_str)
        
        else:
            #print("on the right side")
            # check the sign of the term
            term_sign = eq.rhs.coeff(term)
            #print(f"term: {term}")

            #print(f"term_sign: {term_sign}")
            # Subtract the term from the rhs
            rhs_str = str(eq.rhs - term)
            lhs_str = str(eq.lhs)

            if term_sign > 0:

                lhs_str = lhs_str.__add__(' - ').__add__(term_str)
            else:
                lhs_str = lhs_str.__add__(' + ').__add__(term_str)

        eq_str = Eq(sp.parse_expr(lhs_str, evaluate=False), sp.parse_expr(rhs_str, evaluate=False))

    # remove term_str from the lhs_str
    # Add - term_str to the rhs_str
    #print(f"lhs_str: {lhs_str}")
    
    eq = Eq(eq.lhs - term, eq.rhs - term)
    # return both the string and the equation
    return eq, eq_str

def divide_equation_by_term(eq: Eq) -> Eq:
    """
    Divide both sides of the equation by the coefficient of `term` in the left-hand side.
    If the coefficient is zero, return the original equation.
    """
    # check if the tern is an integer then skip the division
    lhs_str = str(eq.lhs)
    variable = lhs_str.__getitem__(-1)
    coeff = eq.lhs.coeff(variable)
    if coeff == 0:
        return eq
    return Eq(eq.lhs / coeff, eq.rhs / coeff)

def transformation_sequence(eq: Eq, terms: List[sp.Expr], index: int) -> Set[Eq]:
    """
    Generate a set of transformations starting by subtracting the term at position `index`
    from the equation, then applying simplification, division, and subtracting the other terms.
    """
    transformations = set()
    T = terms[index]
    # Step 1: Subtract the chosen term T.
    eq1, eq2 = subtract_term(eq, T)
    #print(f"eq1: {eq1}")
    #print(f"simplify eq1: {simplify(eq1)}")
    transformations.add(eq1)
    transformations.add(eq2)
    # Step 2: Add its simplified version.
    transformations.add(simplify(eq1))
    # Step 3: Divide both sides by the factor of T.
    # Check if there is only one term on the left side

    if len(eq1.lhs.as_ordered_terms()) == 1:
        # check if the term is an integer then skip the division
        if not eq1.lhs.as_ordered_terms()[0].is_integer:
            eq1_div = divide_equation_by_term(eq1)
            #print(f"eq1_div: {eq1_div}")
            eq1_div_simpified = simplify(eq1_div)
            #print(f"eq1_div_simpified: {eq1_div_simpified}")

            transformations.add(eq1_div)   
            transformations.add(eq1_div_simpified)
        

    # Now subtract each of the other terms from eq1.
    for j, term in enumerate(terms):
        if j == index:
            continue
        eq_sub1, eq_sub2 = subtract_term(eq1, term)
        #print(f"eq_sub: {eq_sub1}")
        #print(f"eq_sub2: {eq_sub2}")
        #print(f"simplify eq_sub: {simplify(eq_sub)}")
        transformations.add(eq_sub1)
        transformations.add(eq_sub2)
        transformations.add(simplify(eq_sub1))
        # Check if there is only one term on the left side
        if len(eq_sub1.lhs.as_ordered_terms()) == 1 and not eq_sub1.lhs.as_ordered_terms()[0].is_integer:
            eq1_div = divide_equation_by_term(eq_sub1)
            transformations.add(eq1_div)
            transformations.add(simplify(eq1_div))
            #print(f"eq1_div: {eq1_div}")

        if len(eq_sub2.lhs.as_ordered_terms()) == 1 and not eq_sub2.lhs.as_ordered_terms()[0].is_integer:
            eq2_div = divide_equation_by_term(eq_sub2)
            transformations.add(eq2_div)
            transformations.add(simplify(eq2_div))
            #print(f"eq2_div: {eq2_div}")
    return transformations

def get_terms(eq: Eq) -> List[sp.Expr]:
    """
    Return a list of terms from the equation.
    For an equation of the form e1 + e2 = e3 + e4, we return [e1, e2, e3, e4].
    """
    try:
        lhs_terms = eq.lhs.as_ordered_terms()
        rhs_terms = eq.rhs.as_ordered_terms()
        return lhs_terms + rhs_terms
    except Exception as e:
        raise Exception(f"Error in get_terms: {e}")

def generate_all_transformations(eq: Eq) -> Set[Eq]:
    """
    Starting from an equation, generate many rearrangements by:
      - Iterating over each term in the equation.
      - For each term, subtracting it and then applying a sequence of transformations.
    """
    all_transforms = set()
    try:
        terms = get_terms(eq)
        for i in range(len(terms)):
            transforms = transformation_sequence(eq, terms, i)
            all_transforms.update(transforms)
        return all_transforms
    except Exception as e:
        raise Exception(f"Error in generate_all_transformations: {e}")


# Example usage
if __name__ == "__main__":
    try:
        comparator = EquationComparator()
        # Example list of expression pairs for comparison.
        expression_pairs = [
        ("-4*s + 2", "4*s - r", "2", "-r"),
        ("-2 - 3", "b + 3", "-3*b", "6"),
        ("-8*s", "-r - 2", "-8*s", "-2*r"),
        ("-4*c", "-w/2 - 2", "c", "w/2 - 1/2"),
        ("-4*n - 3", "3*n - y", "-3", "-1*n - y"),
        ("3*c + 1", "-c + m", "c", "(m + 1)/4"),
        ("3*c + 1", "-c + m", "4 + 1", "m"),
        ("4", "3*d + 4*v", "o", "3*d - 4 + 4*v"),
        ("o", "3*d + 4*v - 4", "-4*v", "3*d - 4"),
        #("d*m - d*x - x", "m", "m", "((d + 1)*x)/d - 1"),
        ("r - 4", "12", "r - 4", "-sqrt(12)"),
        ("x + 4", "-x - 5*o", "0", "4 + 5*o"),
        ("8*c", "5 - p", "c", "1.6 - p"),
        ("-5*c + 5", "3*c + p", "-8*c", "-5*p"),
        ("b - 3/2", "3/2", "b - 3/2", "-sqrt(3/2)"),
        ("2*p", "-5 - q", "p", "-2"),
        ("2*p", "-9*r - 5", "p", "-9/2 - 5/2"),
        ("-8*s", "-r - 2", "4*s", "r/2"),
        ("r*(2 - 3)", "0", "r", "0"),
        ("3*r + 2", "2*x", "1.5*r", "x"),
        ("2*m + 5", "0", "3*m + 5", "0"),
        ("2*x", "1", "x", "0"),
        #("-4*o - 4*y", "2*y + 4", "0", "-1"),
        ("-f - 1", "0", "-5*f - 4", "0"),
        ("x*(1 + 2*a)", "-4", "x", "-4/(1 + 2*a)"),
        ("-4*x", "7 + a*x", "x", "-7/4 + a*x"),
        ("-9*t - 2", "5*o", "(-9*t)/5 - 2/5", "o"),
        ("x", "3/2*x", "1", "3/2")
    ]
    #    df_data = pd.read_json('0_9999_v7.json')
        i = 1
    #    for index, row in df_data.head(10).iterrows():
    #        exp1 = row['t0']
    #        exp2 = row['t1']
    #        eq1_str = parse_latex(exp1)
    #        eq2_str = parse_latex(exp2)
    #        #eq1_str = f"{parsed_exp1.lhs}={parsed_exp1.rhs}"
    #        #eq2_str = f"{parsed_exp2.lhs}={parsed_exp2.rhs}"

        for exp1, exp2, exp3, exp4 in expression_pairs:
            
            eq1_str = f"{exp1}={exp2}"
            eq2_str = f"{exp3}={exp4}"
        
            # Compare two sample equations.
            result = comparator.compare(eq1_str, eq2_str)
            print(f" {i} -Comparison Results for: \n {eq1_str} and \n {eq2_str}:")
            i = i + 1
            print("Comparison Results:")
            for key, value in result.items():
                print(f"{key}: {value}")
            print("--------------------------------------")
    except Exception as e:
        print(f"Error during comparison: {e}")