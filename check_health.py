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

import ast
import sys



class CallTreeVisitor(ast.NodeVisitor):
    def __init__(self, call_tree, defined_functions):
        self.call_tree = call_tree
        self.defined_functions = defined_functions

    def visit_Call(self, node):
        if hasattr(node, 'parent_function'):
            if isinstance(node.func, ast.Name):
                called_function = node.func.id
            elif isinstance(node.func, ast.Attribute):
                called_function = node.func.attr
            else:
                called_function = None

            if called_function and called_function in self.defined_functions:
                if node.parent_function not in self.call_tree:
                    self.call_tree[node.parent_function] = set()
                self.call_tree[node.parent_function].add(called_function)
        self.generic_visit(node)


def print_call_tree(call_tree, current_function=None, depth=0, visited=None):
    if visited is None:
        visited = set()

    if current_function is None:
        for function in call_tree:
            if function not in visited:
                print_call_tree(call_tree, function, depth, visited)
    else:
        print("  " * depth + current_function)
        visited.add(current_function)
        if current_function in call_tree:
            for child_function in call_tree[current_function]:
                if child_function not in visited:
                    print_call_tree(call_tree, child_function, depth + 1, visited)


def create_call_tree(file_path):
    with open(file_path, "r") as file:
        code = file.read()

    tree = ast.parse(code)
    call_tree = {}
    defined_functions = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            node.parent_function = node.name
            defined_functions.add(node.name)

        for child_node in ast.iter_child_nodes(node):
            if hasattr(node, 'parent_function'):
                child_node.parent_function = node.parent_function

    CallTreeVisitor(call_tree, defined_functions).visit(tree)
    print_call_tree(call_tree)



import re


class FunctionSummaryVisitor(ast.NodeVisitor):
    def __init__(self):
        self.function_ranges = []

    def visit_FunctionDef(self, node):
        self.function_ranges.append((node.lineno, node.end_lineno))
        self.generic_visit(node)


def remove_comments(line):
    line = re.sub(r"#.*", "", line)
    return line


def summarize_code(file_path):
    with open(file_path, "r") as file:
        code_lines = file.readlines()

    tree = ast.parse(''.join(code_lines))
    visitor = FunctionSummaryVisitor()
    visitor.visit(tree)

    in_function = False
    for i, line in enumerate(code_lines, start=1):
        line = remove_comments(line)
        for start, end in visitor.function_ranges:
            if i == start:
                in_function = True
            elif i == end:
                in_function = False
                break

        if in_function and i != start:
            if i == start + 1:
                print("    ...")
        else:
            if line.strip():
                print(line, end='')


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

    create_call_tree(file_path)

    summarize_code(file_path)
