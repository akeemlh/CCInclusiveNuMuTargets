import sys
import ROOT
import math

# --- Configuration ---
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Colors for terminal output
GREY = "\033[90m"
RED = "\033[91m"
GREEN = "\033[92m"
RESET = "\033[0m"


def get_color(is_match):
    return GREEN if is_match else RED

def format_line(name, status, comment=""):
    color = get_color(status)
    msg = f"{color}{name}{RESET}"
    if comment:
        msg += f" {GREY}({comment}){RESET}"
    return msg

def get_obj_safe(directory, name):
    if directory and directory.GetListOfKeys().Contains(name):
        return directory.Get(name)
    return None

def compare_values(val1, val2):
    try:
        if isinstance(val1, float) and isinstance(val2, float):
            return math.isclose(val1, val2, rel_tol=1e-5)
        return val1 == val2
    except:
        return False


def print_progress(iteration, total, branch_name, suffix='', bar_length=30):

    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(bar_length * iteration // total)
    bar = '=' * filled_length + '-' * (bar_length - filled_length)
    
    # Shorten branch name if it's too long
    display_name = (branch_name[:25] + '..') if len(branch_name) > 25 else branch_name
    
    output_line = f"\rProgress: |{bar}| {percent}% - Checking: {display_name} {suffix}"
    
    sys.stdout.write(output_line)
    sys.stdout.flush()

def check_branch_contents(tree1, tree2, branch_name):
    """
    Checks entries and prints the final PASS/FAIL status immediately 
    to the console before returning the data for the final summary.
    """
    n1 = tree1.GetEntries()
    n2 = tree2.GetEntries()
    
    # Immediate Fail for entry mismatch
    if n1 != n2:
        msg = f"{n1} vs {n2} entries"
        # Print a static failure line immediately
        sys.stdout.write(f"\rProgress: |{'='*30}| 100.0% - Checking: {branch_name} - {RED}FAILED{RESET} ({msg})\n")
        sys.stdout.flush()
        return False, f"{msg} mismatch"

    total_checked = n1
    mismatch_count = 0
    
    # Initialize progress bar
    print_progress(0, total_checked, branch_name, suffix='')

    for i in range(total_checked):
        tree1.GetEntry(i)
        tree2.GetEntry(i)
        try:
            leaf1 = tree1.GetLeaf(branch_name)
            leaf2 = tree2.GetLeaf(branch_name)
            if not leaf1 or not leaf2:
                sys.stdout.write(f"\rProgress: |{'='*30}| 100.0% - Checking: {branch_name} - {RED}ERROR{RESET} (Complex object)\n")
                return False, "Complex object or leaf not found"
            v1 = leaf1.GetValue()
            v2 = leaf2.GetValue()
            """ if hasattr(tree1, branch_name):
                v1 = getattr(tree1, branch_name)
                v2 = getattr(tree2, branch_name)
            else:
                leaf1 = tree1.GetLeaf(branch_name)
                leaf2 = tree2.GetLeaf(branch_name)
                if not leaf1 or not leaf2:
                    sys.stdout.write(f"\rProgress: |{'='*30}| 100.0% - Checking: {branch_name} - {RED}ERROR{RESET} (Complex object)\n")
                    return False, "Complex object or leaf not found"
                v1 = leaf1.GetValue()
                v2 = leaf2.GetValue() """

            if not compare_values(v1, v2):
                mismatch_count += 1
                
        except Exception as e:
            sys.stdout.write(f"\rProgress: |{'='*30}| 100.0% - Checking: {branch_name} - {RED}ERROR{RESET} ({e})\n")
            return False, f"Error: {e}"
        
        # Update progress bar every 200 entries
        if (i + 1) % 200 == 0:
            print_progress(i + 1, total_checked, branch_name, suffix='')

    # -- FINAL STATUS PRINTING --
    if mismatch_count == 0:
        # Green Passed Message
        final_suffix = f"- {GREEN}PASSED{RESET}"
        comment = ""
        match = True
    else:
        # Red Failed Message
        final_suffix = f"- {RED}FAILED{RESET} ({mismatch_count} mismatches)"
        comment = f"{mismatch_count}/{total_checked} entries have value mismatch"
        match = False

    print_progress(total_checked, total_checked, branch_name, suffix=final_suffix)
    sys.stdout.write('\n')
    sys.stdout.flush()

    return match, comment
# ----------------------------------------------------------------------

def recurse_structure(obj1, obj2, name, depth=0):
    indent = "    " * depth
    children_lines = []
    children_match = True

    keys1 = set()
    keys2 = set()
    
    is_tree = isinstance(obj1, ROOT.TTree) or isinstance(obj2, ROOT.TTree)

    if obj1:
        if is_tree:
            keys1 = set(b.GetName() for b in obj1.GetListOfBranches())
        else:
            keys1 = set(k.GetName() for k in obj1.GetListOfKeys())
            
    if obj2:
        if is_tree:
            keys2 = set(b.GetName() for b in obj2.GetListOfBranches())
        else:
            keys2 = set(k.GetName() for k in obj2.GetListOfKeys())

    all_keys = sorted(list(keys1.union(keys2)))

    for key in all_keys:
        sub_obj1 = None
        sub_obj2 = None
        
        if is_tree:
            if key in keys1: sub_obj1 = obj1.GetBranch(key)
            if key in keys2: sub_obj2 = obj2.GetBranch(key)
        else:
            sub_obj1 = get_obj_safe(obj1, key)
            sub_obj2 = get_obj_safe(obj2, key)

        if sub_obj1 and not sub_obj2:
            children_match = False
            children_lines.append(f"{indent}├──{format_line(key, False, 'Found in file1.root, absent in file2.root')}")

        elif not sub_obj1 and sub_obj2:
            children_match = False
            children_lines.append(f"{indent}├──{format_line(key, False, 'Absent in file1.root, found in file2.root')}")

        else:
            if isinstance(sub_obj1, ROOT.TBranch):
                match, comment = check_branch_contents(obj1, obj2, key)
                if not match:
                    children_match = False
                children_lines.append(f"{indent}├──{format_line(key, match, comment)}")
            
            elif isinstance(sub_obj1, (ROOT.TDirectory, ROOT.TTree)):
                sub_match, sub_lines = recurse_structure(sub_obj1, sub_obj2, key, depth + 1)
                if not sub_match:
                    children_match = False
                children_lines.extend(sub_lines)
            
            else:
                children_lines.append(f"{indent}├──{format_line(key, True, 'Object present')}")

    final_output = []
    
    if depth > 0:
        parent_indent = indent[:-4] if depth > 0 else ""
        final_output.append(f"{parent_indent}{format_line(name, children_match, '')}")
    
    final_output.extend(children_lines)

    return children_match, final_output

def main():
    if len(sys.argv) != 3:
        print("Usage: python root_diff_live.py <file1.root> <file2.root>")
        sys.exit(1)

    path1 = sys.argv[1]
    path2 = sys.argv[2]

    f1 = ROOT.TFile.Open(path1)
    f2 = ROOT.TFile.Open(path2)

    if not f1 or f1.IsZombie():
        print(f"Error opening file 1: {path1}")
        sys.exit(1)
    if not f2 or f2.IsZombie():
        print(f"Error opening file 2: {path2}")
        sys.exit(1)

    print(f"Comparing:\n 1: {path1}\n 2: {path2}\n")
    print("--- Live Check Log ---")

    match, lines = recurse_structure(f1, f2, "ROOT")

    print("\n" + "="*40)
    print("--- Final Tree Summary ---")
    print("="*40)
    for line in lines:
        print(line)
        
    f1.Close()
    f2.Close()

if __name__ == "__main__":
    main()