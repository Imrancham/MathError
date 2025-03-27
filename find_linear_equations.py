import sympy as sp
from sympy import Eq, solve, simplify

# Define symbols
s, r, b, c, w, n, y, m, d, v, o, p, q, x, a, f, t = sp.symbols('s r b c w n y m d v o p q x a f t')

# List of expression pairs
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
    ("d*m - d*x - x", "m", "m", "((d + 1)*x)/d - 1"),
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
    ("3*x + 12", "3*(x + 4)", "3*x + 12", "33*x + 12"),
    ("2*m + 5", "0", "3*m + 5", "0"),
    ("2*x", "1", "x", "0"),
    ("2*o", "-3*c + 4", "o", "-1"),
    ("-4*o - 4*y", "2*y + 4", "0", "-1"),
    ("-f - 1", "0", "-5*f - 4", "0"),
    ("x*(1 + 2*a)", "-4", "x", "-4/(1 + 2*a)"),
    ("-4*x", "7 + a*x", "x", "-7/4 + a*x"),
    ("-9*t - 2", "5*o", "(-9*t)/5 - 2/5", "o"),
    ("x", "3/2*x", "1", "3/2")
]

# Function to check if an expression is linear
def is_linear(expr, variables):
    """
    Check if an expression is linear with respect to the given variables.
    """
    return expr.as_poly().degree() == 1




# Function to extract features
def extract_features(expr1, expr2, expr3, expr4):
    # Parse expressions
    eq1 = Eq(sp.sympify(expr1), sp.sympify(expr2))
    eq2 = Eq(sp.sympify(expr3), sp.sympify(expr4))
    
    # Syntactic features
    syntactic_features = {
        "depth_expr1": sp.srepr(eq1).count("("),  # Depth of expression 1
        "depth_expr2": sp.srepr(eq2).count("("),  # Depth of expression 2
        "num_operators_expr1": str(eq1).count("+") + str(eq1).count("-") + str(eq1).count("*") + str(eq1).count("/"),
        "num_operators_expr2": str(eq2).count("+") + str(eq2).count("-") + str(eq2).count("*") + str(eq2).count("/"),
    }
    
    # Semantic features
    eq1_diff = simplify(eq1.lhs - eq1.rhs)
    eq2_diff = simplify(eq2.lhs - eq2.rhs)
    is_equivalent = False
    if not isinstance(eq1_diff, bool) and not isinstance(eq2_diff, bool):
        is_equivalent = eq1_diff == eq2_diff

    semantic_features = {
        "is_equivalent": is_equivalent,  # Check equivalence
        "is_simplified": simplify(eq1) == simplify(eq2),  # Check simplification
    }
    
    # Contextual features
    contextual_features = {
        "is_linear": is_linear(eq1.lhs - eq1.rhs) and is_linear(eq2.lhs - eq2.rhs),  # Check if equations are linear
        "is_solution_correct": solve(eq1) == solve(eq2),  # Check if solutions match
    }
    
    return {**syntactic_features, **semantic_features, **contextual_features}

# Extract features for all pairs
for idx, (expr1, expr2, expr3, expr4) in enumerate(expression_pairs):
    try:
        print(f"Pair {idx + 1}:")
        features = extract_features(expr1, expr2, expr3, expr4)
        for key, value in features.items():
            print(f"  {key}: {value}")
        print("\n")
    except Exception as e:
        print(f"Error processing pair {idx + 1}: {e}")