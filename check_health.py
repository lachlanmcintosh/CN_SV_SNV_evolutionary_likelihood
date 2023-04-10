import ast
from collections import defaultdict

class FunctionAnalyzer(ast.NodeVisitor):
    def __init__(self):
        self.global_vars = defaultdict(set)
        self.global_scope_vars = set()

    def visit_FunctionDef(self, node):
        defined_vars = set()
        used_vars = set()

        for child in ast.walk(node):
            if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Load):
                used_vars.add(child.id)
            elif isinstance(child, (ast.Assign, ast.AnnAssign, ast.arg)):
                targets = [child.target] if isinstance(child, ast.AnnAssign) else child.targets if isinstance(child, ast.Assign) else [child]
                for target in targets:
                    if isinstance(target, ast.Name):
                        defined_vars.add(target.id)
                    elif isinstance(target, ast.arg):
                        defined_vars.add(target.arg)

        globals_in_function = (used_vars - defined_vars) & self.global_scope_vars
        if globals_in_function:
            self.global_vars[node.name] = globals_in_function

        self.generic_visit(node)

    def visit_Assign(self, node):
        for target in node.targets:
            if isinstance(target, ast.Name) and isinstance(target.ctx, ast.Store):
                self.global_scope_vars.add(target.id)
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        if isinstance(node.target, ast.Name) and isinstance(node.target.ctx, ast.Store):
            self.global_scope_vars.add(node.target.id)
        self.generic_visit(node)

def find_unused_functions(source_code):
    tree = ast.parse(source_code)

    defined_functions = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            defined_functions.add(node.name)

    called_functions = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                called_functions.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                called_functions.add(node.func.attr)

    return defined_functions - called_functions

if __name__ == "__main__":
    file_path = "run_simulation_and_analysis2.py"

    with open(file_path, "r") as file:
        code = file.read()

    # Analyze global variable usage
    tree = ast.parse(code)
    analyzer = FunctionAnalyzer()
    analyzer.visit(tree)

    for func_name, globals_used in analyzer.global_vars.items():
        print(f"Function '{func_name}' uses true global variables: {', '.join(globals_used)}")

    # Find unused functions
    unused_functions = find_unused_functions(code)
    print(f"Unused functions: {', '.join(unused_functions)}")

