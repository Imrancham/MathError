import sympy as sp
import difflib
import pandas as pd
import re
#from sympy.parsing.latex import parse_latex
from pylatexenc.latex2text import LatexNodes2Text
from sympy.parsing.sympy_parser import standard_transformations, implicit_multiplication_application


def parse_expression(expr):
    try:
        #print("trying to parse: ", expr)
        expr = LatexNodes2Text().latex_to_text(expr)
        expr = expr.replace(" ", "")  # Remove spaces
        expr = expr.replace("^", "**")  # Fix exponentiation issue
        expr = expr.replace(" (", "*(")  # Fix square root issue
        #expr = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', expr)  # Convert 2x to 2*x
        #expr = re.sub(r'([a-zA-Z])\(', r'\1*(', expr)  # Fix h (1-q) → h*(1-q)
        return sp.parse_expr(expr, transformations=(standard_transformations + (implicit_multiplication_application,)))
    except Exception      as e:
        print("error: ", e)
        return None


def compare_expressions(expr1, expr2):

   

    
    parsed_expr1 = parse_expression(expr1)
    parsed_expr2 = parse_expression(expr2)
    

    #print("ex1: ", expr1,"+", parsed_expr1)
    #print("ex2: ", expr2,"+", parsed_expr2)
    
    if parsed_expr1 is None or parsed_expr2 is None:
        return "Invalid Expression"
    
    try:
        ex  = sp.simplify(parsed_expr1)
        if sp.simplify(ex - parsed_expr2) == 0:
            return "No Error"
    except TypeError:
        return "Invalid Expression"
    
    diffs = list(difflib.ndiff(expr1, expr2))
    print("diffs: ", diffs)
    if '-' in ''.join(diffs) or '+' in ''.join(diffs):
        return "Possible Term Manipulation Error"
    
    return "Unknown Error"

def annotate_errors(df):
    df['Error Type'] = df.apply(lambda row: compare_expressions(row['t0'], row['t1']), axis=1)
    return df

# Load dataset (replace 'dataset.csv' with actual filename)
df = pd.read_json("0_9999_v7.json") 

# Annotate errors
df = annotate_errors(df.head(10))

# Save annotated dataset
df.to_csv("annotated_dataset.csv", index=False)

print(df)
