import os
import json
import pandas as pd
import sympy as sp
import logging
import re
from sympy import Eq, simplify
from sympy.parsing.sympy_parser import parse_expr
from functools import lru_cache
from typing import Tuple
from sympy.parsing.latex import parse_latex

# ====== EquationPreprocessor, EquationParser, and LinearEquation (same as before) ======



class EquationParser:
    @lru_cache(maxsize=128)
    def parse(self, eq_str: str) -> Tuple[Eq, Eq]:
        try:
            lhs, rhs = map(str.strip, eq_str.split('=', 1))
        except ValueError:
            raise ValueError("Equation must contain exactly one '='")

        lhs_expr = parse_expr(lhs, evaluate=False)
        rhs_expr = parse_expr(rhs, evaluate=False)

        return (
            Eq(simplify(lhs_expr), simplify(rhs_expr)),  # simplified
            Eq(lhs_expr, rhs_expr, evaluate=False)       # original
        )


class LinearEquation:
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self._variables = sp.symbols('s r b c w n y m d v o p q x a f t h j k l i e g z u X T U V W A B C D E F G H I J K L M N O P Q R S T U V W X Y Z')

    def is_linear(self, expr: sp.Expr) -> bool:
        try:
            poly = expr.as_poly()
            if poly is None or poly.degree() != 1:
                return False

            if any(term.is_Pow and term.exp.is_Number and term.exp < 0 for term in expr.atoms(sp.Pow)):
                return False

            if any(term.is_Pow and term.exp == sp.Rational(1, 2) for term in expr.atoms(sp.Pow)):
                return False

            return True
        except Exception as e:
            self.logger.error(f"Linear check failed: {str(e)}")
            return False


# ====== Main loop to process all JSON files in a directory ======

def process_equations_from_directory(directory_path: str):
    parser = EquationParser()
    linear_checker = LinearEquation()

    all_linear_equations = []
    i = 0
    all_inequalities = []
    for filename in os.listdir(directory_path):
        if filename.endswith(".json"):
            df_data = pd.read_json(os.path.join(directory_path, filename))


            parser = EquationParser()
            linear_checker = LinearEquation()
            if 'category' not in df_data.columns or 't0' not in df_data.columns:
                logging.warning(f"Skipping {filename} because it lacks required columns.")
                continue


            for index, row in df_data[df_data['category'] == 2].iterrows():
                
                try:
                    eq_str = row['t0']
                    # if it contain an inequality, skip it
                    if re.search(r'[<>]', eq_str):
                        all_inequalities.append(row)
                        continue

                    # split the equation at the first '='

                    eq1_str = parse_latex(eq_str)
                    eq1_str = f"{eq1_str.lhs}={eq1_str.rhs}"


                    eq1_eval, eq1_raw = parser.parse(eq1_str)
                

                
                    canonical_expr = eq1_eval.lhs - eq1_eval.rhs

                    if linear_checker.is_linear(canonical_expr):
                        all_linear_equations.append({
                            "source_file": "first_file.json",
                        
                            "original": eq_str,
                            "s1": row['t1'],
                            "assumption": row['assumptions'],
                            "preprocessed": eq1_str,
                            "canonical": str(canonical_expr)
                        })

                except Exception as e:
                    logging.error(f"Error processing equation: {eq_str} ")
            continue

    df = pd.DataFrame(all_linear_equations)
    df_inequalities = pd.DataFrame(all_inequalities)
    df_inequalities.to_csv(f"all_inequalities.csv", index=False)
    print(f"\nTotal number of linear equations across all files: {len(df)}")
    print(f"\nTotal number of inequalities across all files: {len(df_inequalities)}")

    # Save results
    df.to_csv(f"all_linear_equations_no_inequalities.csv", index=False)
    i+=1
    return df


# ====== Run ======
if __name__ == "__main__":
    folder = r"C:\Users\imran.chamieh\Desktop\data\json"
    df = process_equations_from_directory(folder)
