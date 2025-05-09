import sympy as sp
from sympy import Eq, simplify, symbols
from itertools import combinations
import re

def preprocess_expression(expr_str):
    """Add explicit multiplication operators where needed"""
    expr_str = expr_str.replace('−', '-')
    expr_str = re.sub(r'[\u200B\u2060\uFEFF]', '', expr_str)  # Removes U+200B, U+2060, U+FEFF

    # Add * between: digit-letter, letter-digit, letter-letter
    expr_str = re.sub(r'(\d)([a-zA-Z(])', r'\1*\2', expr_str)
    expr_str = re.sub(r'([a-zA-Z])(\d)', r'\1*\2', expr_str)
    expr_str = re.sub(r'([a-zA-Z])([a-zA-Z])', r'\1*\2', expr_str)
    return expr_str

def parse_equation(eq_str):
    """Parse a string equation into a SymPy Eq object"""
    try:
        lhs_str, rhs_str = eq_str.split('=', 1)
    except ValueError:
        raise ValueError("Equation must contain exactly one '='")

    lhs = sp.sympify(preprocess_expression(lhs_str.strip()))
    rhs = sp.sympify(preprocess_expression(rhs_str.strip()))
        # Parse both sides without automatic simplification
    lhs2 = sp.parse_expr(preprocess_expression(lhs_str.strip()), evaluate=False)
    rhs2 = sp.parse_expr(preprocess_expression(rhs_str.strip()), evaluate=False)


    return Eq(lhs, rhs), Eq(lhs2, rhs2)

def generate_rearrangements(equation, raw_eq):
    """Generate all algebraic rearrangements of an equation"""
    # Move all terms to left side
    all_terms = sp.expand(equation.lhs - equation.rhs)
    
    # Handle special case of 0 = 0
    if all_terms == 0:
        return {equation}
    
    # Get individual terms
    terms = sp.Add.make_args(all_terms)
    n = len(terms)
    
    rearrangements = set()
    
    # Generate all possible splits into left/right sides
    for k in range(1, n):  # Avoid empty splits
        for combo in combinations(terms, k):
            left = sum(combo)
            right = -sum(t for t in terms if t not in combo)
            rearrangements.add(Eq(left, right))
            rearrangements.add(Eq(right, left))  # Commutative

    eq = Eq(all_terms, 0)
    rearrangements.add(eq)
    rearrangements.add(Eq(0, all_terms))  # Commutative
     # Identify all variables dynamically
    variables = set(re.findall(r'[a-zA-Z]+', raw_eq))
    symbols_map = {var: symbols(var) for var in variables}  # Define variables in SymPy
    eq1 = equation
    # Solve the equation for each symbol
    for var in variables:
        res = sp.solve(eq1, symbols_map[var])
        if res:
            rearrangements.add(Eq(symbols_map[var], res[0]))
    
    return rearrangements

def equation_similarity(eq1, eq2):
    """Compute similarity score between two equations based on term overlap."""
    # Extract terms from both equations
    terms1 = set(eq1.lhs.as_ordered_terms() + eq1.rhs.as_ordered_terms())
    terms2 = set(eq2.lhs.as_ordered_terms() + eq2.rhs.as_ordered_terms())

    # Compute Jaccard similarity (|A ∩ B| / |A ∪ B|)
    intersection = len(terms1 & terms2)
    union = len(terms1 | terms2)

    return intersection / union if union != 0 else 0

def find_most_similar_match(rearrangements, target_eq):
    """Find the most similar equation in the rearrangements set to target_eq."""
    best_match = None
    best_score = 0

    for eq in rearrangements:
        score = equation_similarity(eq, target_eq)
        if score > best_score:
            best_score = score
            best_match = eq

    return best_match, best_score

def compare_expressions(exp1_str, exp2_str):
    """Compare two equations and find equivalent forms"""
    try:
        eq1, eq12 = parse_equation(exp1_str)
        eq2, eq22  = parse_equation(exp2_str)
    except Exception as e:
        return f"Error parsing equations: {str(e)}"

    # Generate all rearrangements of exp1
    rearrangements = generate_rearrangements(eq1, exp1_str)
    rearrangements2 = generate_rearrangements(eq12, exp1_str)
    rearrangements.update(rearrangements2)
    most_similar, similarity_score = find_most_similar_match(rearrangements, eq2)

    
    
    # Find matches with exp2
    matches = []
    for eq in rearrangements:
        print(f"Equation:  {eq} - {eq22}")
        print(f"subractions: {simplify(eq.rewrite(sp.Add) - eq22.rewrite(sp.Add))}")
        if simplify(eq.rewrite(sp.Add) == eq22.rewrite(sp.Add)):
            matches.append(eq)
    
    return matches, most_similar, similarity_score

# Example usage
if __name__ == "__main__":
    print("Enter equations using variables and coefficients like '3x + 4 = y'")
    exp1_input = "4=3d+4v"
    exp2_input = "0=3d−4+4v"

    matches, most_similar, similarity_score = compare_expressions(exp1_input, exp2_input)
    print(f"Most similar equation: {most_similar}")
    print(f"Similarity score: {similarity_score:.2f}")
    print(f"\nAll matching rearrangements: {matches}")
    