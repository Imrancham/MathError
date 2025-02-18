from sympy import symbols, Eq, simplify, expand, sympify
import re

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
    
    return equation, list(lhs_terms), list(rhs_terms), symbols_map

def compare_expressions(eq1_str, eq2_str):
    """Compare all expression pairs from eq1 with the left and right side of eq2."""
    
    # Parse both equations
    eq1, lhs_terms1, rhs_terms1, _ = parse_equation(eq1_str)
    eq2, lhs_terms2, rhs_terms2, _ = parse_equation(eq2_str)
    
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

# Example input equations
eq1 = "5x + 9y = 2(x + 3y)"
eq2 = "5x - 3x = 6y + 9y"

compare_expressions(eq1, eq2)
