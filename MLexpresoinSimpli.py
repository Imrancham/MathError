
from sympy import symbols, Eq, simplify, expand, sympify, Add
import re
from itertools import combinations

def format_equation(equation_str):
    """Ensure multiplication is explicit in the equation (e.g., 5x → 5*x, 2(x + y) → 2*(x + y))."""
    equation_str = re.sub(r'(\d)([a-zA-Z(])', r'\1*\2', equation_str)  # Add * between numbers and variables
    return equation_str

def parse_equation(equation_str):
    """Extract expressions from a linear equation string."""
    equation_str = format_equation(equation_str)  # Ensure correct formatting
    lhs_str, rhs_str = equation_str.split("=")
    
    # Identify all variables dynamically
    variables = set(re.findall(r'[a-zA-Z]+', equation_str))
    symbols_map = {var: symbols(var) for var in variables}  # Define variables in SymPy
    
    # Convert string expressions to sympy expressions
    lhs_expr = sympify(lhs_str.strip(), locals=symbols_map)  
    rhs_expr = sympify(rhs_str.strip(), locals=symbols_map)  
    
    # Expand both sides to simplify expressions
    lhs_expr = expand(lhs_expr)
    rhs_expr = expand(rhs_expr)
    
    # Extract individual terms from both sides
    lhs_terms = lhs_expr.as_ordered_terms()
    rhs_terms = rhs_expr.as_ordered_terms()
    
    # Create an equation
    equation = Eq(lhs_expr, rhs_expr)

    # Move all terms to the left side and expand
    all_terms = expand(lhs_expr - rhs_expr)
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
    
    
    return rearrangements, equation, list(lhs_terms), list(rhs_terms), symbols_map

def compare_expressions(eq1_str, eq2_str):
    """Compare all expression pairs from eq1 with the left and right side of eq2."""
    
    # Parse both equations
    _, eq1, lhs_terms1, rhs_terms1, _ = parse_equation(eq1_str)
    _, eq2, lhs_terms2, rhs_terms2, _ = parse_equation(eq2_str)
    
    print(f"Equation 1: {eq1}")
    print(f"Equation 2: {eq2}\n")
    
    # Create all pairs from eq1 expressions
    expressions1 = lhs_terms1 + rhs_terms1
    pairs1 = []
    
    num_expressions = len(expressions1)
    for i in range(num_expressions):
        for j in range(i + 1, num_expressions):
            pair = simplify(expressions1[i] + expressions1[j])  # Simplify sum of two expressions
            pairs1.append(pair)
    
    # Expressions from eq2
    expressions2 = lhs_terms2 + rhs_terms2
    
    # Compare each pair from eq1 with expressions from eq2
    matching_pairs = []
    for pair in pairs1:
        for expr in expressions2:
            if simplify(pair - expr) == 0:  # Check if pair is mathematically equal to an expression in eq2
                matching_pairs.append((pair, expr))
    
    # Print results
    if matching_pairs:
        print("Matching Expressions Found:")
        for pair, expr in matching_pairs:
            print(f"{pair} == {expr}")
    else:
        print("No matching expressions found.")



def compare_expressions2(exp1_str, exp2_str):
    """Compare exp2 to all rearrangements of exp1 using SymPy."""
    # Parse input strings into SymPy equations
    try:
        equation_str1, _  = parse_equation(exp1_str)  # Ensure correct formatting
        equation_str2, _  = parse_equation(exp2_str)  # Ensure correct formatting
        
        # Identify all variables dynamically
        variables = set(re.findall(r'[a-zA-Z]+', equation_str1))
        symbols_map = {var: symbols(var) for var in variables}  # Define variables in SymPy
        
        # Convert string expressions to sympy expressions
        exp1 = sympify(equation_str1, locals=symbols_map)  
        exp2 = sympify(equation_str2, locals=symbols_map)  
        
 
    except Exception as e:
        raise ValueError(f"Error parsing input expressions: {str(e)}")
    
    if not isinstance(exp1, Eq) or not isinstance(exp2, Eq):
        raise ValueError("Inputs must be equations")
    
    # Generate all rearrangements of exp1 takle only the first output of the pqrse_equation function
    all_rearrangements = parse_equation(exp1_str)[0]
    
    # Check if exp2 matches any rearrangement
    matches = []
    for eq in all_rearrangements:
        if simplify(eq.rewrite(Add) - exp2.rewrite(Add)) == 0:
            matches.append(eq)
    
    return matches


# Example usage
if __name__ == "__main__":
    # Define possible symbols (could be extended based on the user's input)
    symbols = symbols('x y a b c d')

    # Input expressions (can include variables and constants)
    exp1_input = input("Enter exp1 (e.g., 'a + b = c + d'): ")
    exp2_input = input("Enter exp2 (e.g., 'c + d - a = b'): ")
    
    try:
        matches = compare_expressions2(exp1_input, exp2_input)
        if matches:
            print("\nMatching rearrangements found:")
            for i, match in enumerate(matches, 1):
                print(f"{i}. {match}")
        else:
            print("\nNo matching rearrangements found.")
        
        compare_expressions(exp1_input, exp2_input)
    
    except Exception as e:
        print(f"Error: {str(e)}")
