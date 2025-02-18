import sympy as sp
from sympy import Eq, simplify
from itertools import combinations

def generate_rearrangements(exp1):
    """Generate all unique algebraic rearrangements of exp1."""
    # Extract terms from both sides of the equation
    lhs = exp1.lhs
    rhs = exp1.rhs
    
    # Move all terms to the left side and expand
    all_terms = sp.expand(lhs - rhs)
    terms = all_terms.as_ordered_terms()
    n_terms = len(terms)
    
    rearrangements = set()
    
    # Generate all possible splits of terms into left and right
    for k in range(1, n_terms):  # Avoid empty splits
        for left_terms in combinations(terms, k):
            left_sum = sum(left_terms)
            right_sum = all_terms - left_sum
            rearranged_eq = Eq(left_sum, -right_sum)
            rearrangements.add(rearranged_eq)
    
    return rearrangements

def compare_expressions(exp1_str, exp2_str, symbols):
    """Compare exp2 to all rearrangements of exp1 using SymPy."""
    # Parse input strings into SymPy equations
    try:
        exp1 = sp.sympify(exp1_str, locals=symbols)
        exp2 = sp.sympify(exp2_str, locals=symbols)
    except Exception as e:
        raise ValueError(f"Error parsing input expressions: {str(e)}")
    
    if not isinstance(exp1, Eq) or not isinstance(exp2, Eq):
        raise ValueError("Inputs must be equations")
    
    # Generate all rearrangements of exp1
    all_rearrangements = generate_rearrangements(exp1)
    
    # Check if exp2 matches any rearrangement
    matches = []
    for eq in all_rearrangements:
        if simplify(eq.rewrite(sp.Add) - exp2.rewrite(sp.Add)) == 0:
            matches.append(eq)
    
    return matches

# Example usage
if __name__ == "__main__":
    # Define possible symbols (could be extended based on the user's input)
    symbols = sp.symbols('x y a b c d')

    # Input expressions (can include variables and constants)
    exp1_input = input("Enter exp1 (e.g., 'a + b = c + d'): ")
    exp2_input = input("Enter exp2 (e.g., 'c + d - a = b'): ")
    
    try:
        matches = compare_expressions(exp1_input, exp2_input, symbols)
        if matches:
            print("\nMatching rearrangements found:")
            for i, match in enumerate(matches, 1):
                print(f"{i}. {match}")
        else:
            print("\nNo matching rearrangements found.")
    
    except Exception as e:
        print(f"Error: {str(e)}")
