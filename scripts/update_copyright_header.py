#!/usr/bin/env python3
"""
20260321 Sat   update_copyright_header.py
Updates copyright headers in specified C/C++ and Python files
based on copyright.tmpl template.

Target files:
    - include/* (all C/C++ files)
    - src/* (all C/C++ files)
    - scripts/run.py
    - setup.py

Excluded files:
    - include/kseq.h
    - include/sse2neon.h
    - include/ssw.h
    - src/sortmerna/ssw_example.c
    - src/sortmerna/ssw.c

Usage:
    python3 scripts/update_copyright_header.py                    # Update all target files
    python3 scripts/update_copyright_header.py --dry-run          # Preview changes
    python3 scripts/update_copyright_header.py --list             # List target files
    python3 scripts/update_copyright_header.py include/options.hpp  # Specific file
"""

import re
import sys
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Set
from subprocess import run as subprocess_run

DEFAULT_TEMPLATE = "copyright.tmpl"

# Target directories to scan recursively
TARGET_DIRS = [
    "include",
    "src",
]

# Target specific files (relative to repo root)
TARGET_FILES = [
    "scripts/run.py",
    "setup.py",
]

EXCLUDE_FILES = [
    "include/kseq.h",
    "include/sse2neon.h",
    "include/ssw.h",
    "src/sortmerna/ssw_example.c",
    "src/sortmerna/ssw.c"
]

# File extensions and their comment styles
FILE_CONFIGS = {
    # C/C++ files: /* */ block comments
    ".h":   {"style": "c_block", "remove_parblock": False},
    ".hpp": {"style": "c_block", "remove_parblock": False},
    ".hxx": {"style": "c_block", "remove_parblock": False},
    ".c":   {"style": "c_block", "remove_parblock": False},
    ".cpp": {"style": "c_block", "remove_parblock": False},
    ".cxx": {"style": "c_block", "remove_parblock": False},
    ".cc":  {"style": "c_block", "remove_parblock": False},
    # Python files: # line comments (remove @parblock markers)
    ".py":  {"style": "python", "remove_parblock": True},
}

EXCLUDE_DIRS = {".git", "build", "cmake-build", "venv", "__pycache__", 
                "node_modules", "bin", "obj", "out", "dist", ".eggs"}
# ─────────────────────────────────────────────────────────────────────────────


def get_repo_root() -> Path:
    """Get the root directory of the git repository."""
    result = subprocess_run(
        ["git", "rev-parse", "--show-toplevel"],
        capture_output=True, text=True, check=False
    )
    if result.returncode == 0:
        return Path(result.stdout.strip())
    return Path.cwd()


def load_template(template_path: str, repo_root: Path) -> str:
    """Load and return the copyright template content."""
    path = Path(template_path)
    if not path.is_absolute():
        path = repo_root / template_path
    
    if not path.exists():
        raise FileNotFoundError(f"Template file not found: {path}")
    
    with open(path, "r", encoding="utf-8") as f:
        return f.read().strip()


def extract_existing_header_c(content: str) -> Optional[str]:
    """Extract existing /* */ header block from C/C++ file content."""
    match = re.match(r'^\s*/\*.*?\*/', content, re.DOTALL)
    return match.group(0) if match else None


def extract_existing_header_python(content: str) -> Optional[str]:
    """
    Extract existing # comment header block from Python file content.
    Handles consecutive # lines at the start of file.
    """
    lines = content.split('\n')
    header_lines = []
    in_header = True
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        if in_header:
            if stripped.startswith('#') and not stripped.startswith('#!'):
                header_lines.append(line)
            elif stripped == '' and header_lines:
                header_lines.append(line)
            else:
                in_header = False
                if i == 0 and lines[0].startswith('#!'):
                    header_lines = []
                break
        else:
            break
    
    if header_lines:
        while header_lines and header_lines[-1].strip() == '':
            header_lines.pop()
        return '\n'.join(header_lines)
    
    return None


def clean_template_for_python(template: str) -> str:
    """
    Remove @parblock and @endparblock markers from template for Python files.
    Also cleans up extra blank lines that may result.
    """
    # Remove @parblock and @endparblock lines
    cleaned = re.sub(r'^\s*@parblock\s*$', '', template, flags=re.MULTILINE)
    cleaned = re.sub(r'^\s*@endparblock\s*$', '', cleaned, flags=re.MULTILINE)
    
    # Clean up multiple consecutive blank lines (max 1 blank line)
    cleaned = re.sub(r'\n{3,}', '\n\n', cleaned)
    
    return cleaned.strip()


def generate_header_c(template: str, current_year: int, contributors: List[str]) -> str:
    """Generate C/C++ style header with /* */ block comments."""
    header = template.replace("${current_year}", str(current_year))
    header = header.replace("${contributors}", "\n".join(c for c in contributors))
    
    header = header.strip()
    if not header.startswith("/*"):
        header = f"/*\n{header}"
    if not header.endswith("*/"):
        header = f"{header}\n*/"
    
    return header


def generate_header_python(template: str, current_year: int, contributors: List[str]) -> str:
    """Generate Python style header with # line comments."""
    # Clean template for Python (remove @parblock markers)
    cleaned_template = clean_template_for_python(template)
    
    header = cleaned_template.replace("${current_year}", str(current_year))
    header = header.replace("${contributors}", "\n".join(c for c in contributors))
    
    lines = header.split('\n')
    commented_lines = []
    for line in lines:
        if line.strip():
            commented_lines.append(f"# {line}")
        else:
            commented_lines.append("#")
    
    return '\n'.join(commented_lines)


def update_file_header(
    filepath: Path,
    template: str,
    style: str,
    remove_parblock: bool,
    current_year: int,
    contributors: List[str],
    dry_run: bool = False
) -> bool:
    """Update the copyright header in a single source file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            original_content = f.read()
    except Exception as e:
        print(f"✗ Error reading {filepath}: {e}", file=sys.stderr)
        return False
    
    if style == "c_block":
        existing_header = extract_existing_header_c(original_content)
        new_header = generate_header_c(template, current_year, contributors)
    elif style == "python":
        existing_header = extract_existing_header_python(original_content)
        new_header = generate_header_python(template, current_year, contributors)
    else:
        print(f"✗ Unknown style '{style}' for {filepath}", file=sys.stderr)
        return False
    
    # Replace or prepend header
    if existing_header:
        new_content = original_content.replace(existing_header, new_header, 1)
    else:
        lines = original_content.split("\n", 1)
        # Preserve shebang for Python files
        if lines[0].startswith("#!"):
            new_content = f"{lines[0]}\n{new_header}\n\n{lines[1] if len(lines) > 1 else ''}"
        else:
            new_content = f"{new_header}\n\n{original_content}"
    
    if new_content == original_content:
        return False
    
    if not dry_run:
        try:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(new_content)
            print(f"✓ Updated: {filepath}")
        except Exception as e:
            print(f"✗ Error writing {filepath}: {e}", file=sys.stderr)
            return False
    else:
        print(f"⊘ Would update: {filepath}")
    
    return True


def normalize_path(filepath: Path, repo_root: Path) -> str:
    """Normalize a file path to a relative path string for comparison."""
    try:
        return str(filepath.relative_to(repo_root))
    except ValueError:
        return str(filepath)


def is_excluded(filepath: Path, repo_root: Path, exclude_set: Set[str]) -> bool:
    """Check if a file is in the exclusion list."""
    rel_path = normalize_path(filepath, repo_root)
    return rel_path in exclude_set


def find_target_files(repo_root: Path, 
                      specific_files: Optional[List[str]] = None) -> Dict[str, List[Path]]:
    """
    Find target files from specified directories and files.
    Excludes files listed in EXCLUDE_FILES.
    
    Returns dict: {"c_block": [files], "python": [files]}
    """
    result = {"c_block": [], "python": []}
    
    # Build exclusion set with normalized paths
    exclude_set = set()
    for excl in EXCLUDE_FILES:
        excl_path = repo_root / excl
        if excl_path.exists():
            exclude_set.add(normalize_path(excl_path, repo_root))
        else:
            exclude_set.add(excl)
    
    if specific_files:
        # Use provided specific file list only
        for f in specific_files:
            path = Path(f)
            if not path.is_absolute():
                path = repo_root / path
            if path.exists():
                # Check exclusion
                if is_excluded(path, repo_root, exclude_set):
                    print(f"⊘ Excluded: {path.relative_to(repo_root)}")
                    continue
                
                suffix = path.suffix.lower()
                if suffix in FILE_CONFIGS:
                    style = FILE_CONFIGS[suffix]["style"]
                    result[style].append(path)
        return result
    
    # Scan target directories recursively
    for target_dir in TARGET_DIRS:
        dir_path = repo_root / target_dir
        if not dir_path.exists():
            print(f"⚠️  Directory not found: {dir_path}")
            continue
        
        for filepath in dir_path.rglob("*"):
            if not filepath.is_file():
                continue
            
            # Check exclusion
            if is_excluded(filepath, repo_root, exclude_set):
                continue
            
            suffix = filepath.suffix.lower()
            if suffix not in FILE_CONFIGS:
                continue
            
            if any(excl in filepath.parts for excl in EXCLUDE_DIRS):
                continue
            
            style = FILE_CONFIGS[suffix]["style"]
            result[style].append(filepath)
    
    # Add specific target files
    for target_file in TARGET_FILES:
        filepath = repo_root / target_file
        if filepath.exists():
            # Check exclusion
            if is_excluded(filepath, repo_root, exclude_set):
                print(f"⊘ Excluded: {filepath.relative_to(repo_root)}")
                continue
            
            suffix = filepath.suffix.lower()
            if suffix in FILE_CONFIGS:
                style = FILE_CONFIGS[suffix]["style"]
                if filepath not in result[style]:
                    result[style].append(filepath)
            else:
                print(f"⚠️  Unsupported file type: {filepath}")
        else:
            print(f"⚠️  Target file not found: {filepath}")
    
    # Deduplicate and sort
    for style in result:
        seen = set()
        unique = []
        for f in result[style]:
            if f not in seen:
                seen.add(f)
                unique.append(f)
        result[style] = sorted(unique, key=lambda x: str(x))
    
    return result


def update_copyright_header(
    specific_files: Optional[List[str]] = None,
    template_path: str = DEFAULT_TEMPLATE,
    contributors: Optional[List[str]] = None,
    dry_run: bool = False,
    list_only: bool = False
) -> int:
    """
    Main function: Update copyright headers in target files.
    
    Args:
        specific_files: Specific files to process (None = all target files)
        template_path: Path to copyright template file
        contributors: List of contributor entries
        dry_run: If True, only report changes without writing
        list_only: If True, only list files that would be updated
    
    Returns:
        Number of files that would be/were updated
    """
    #if contributors is None:
    #    contributors = []
    
    repo_root = get_repo_root()
    current_year = datetime.now().year
    
    # Load templates
    try:
        template = load_template(template_path, repo_root)
    except FileNotFoundError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        return 1
    
    # Find target files
    files_by_style = find_target_files(repo_root, specific_files)
    
    total_files = len(files_by_style["c_block"]) + len(files_by_style["python"])
    
    if not total_files:
        print("⊘ No target files found.")
        print(f"\nConfigured target directories: {TARGET_DIRS}")
        print(f"Configured target files: {TARGET_FILES}")
        print(f"Configured excluded files: {EXCLUDE_FILES}")
        return 0
    
    if list_only:
        print(f"📁 Target files that would be updated ({total_files} total):")
        print(f"\n   C/C++ Files ({len(files_by_style['c_block'])}):")
        for f in files_by_style["c_block"]:
            print(f"      • {f.relative_to(repo_root)}")
        print(f"\n   Python Files ({len(files_by_style['python'])}):")
        for f in files_by_style["python"]:
            print(f"      • {f.relative_to(repo_root)}")
        print(f"\n⊘ Excluded files ({len(EXCLUDE_FILES)}):")
        for excl in EXCLUDE_FILES:
            print(f"      • {excl}")
        return total_files
    
    print(f"🔄 Processing {total_files} target file(s)...")
    print(f"📁 Target directories: {', '.join(TARGET_DIRS)}")
    print(f"📁 Target files: {', '.join(TARGET_FILES)}")
    print(f"⊘ Excluded files: {len(EXCLUDE_FILES)}")
    print(f"📄 Template: {template_path}")
    print(f"📅 Year: {current_year}")
    print(f"🐍 Python files: @parblock markers will be removed")
    if dry_run:
        print("⚠️  DRY RUN - No files will be modified")
    print()
    
    updated_count = 0
    
    # Update C/C++ files
    if files_by_style["c_block"]:
        print(f"📝 Updating {len(files_by_style['c_block'])} C/C++ file(s)...")
        for filepath in files_by_style["c_block"]:
            if update_file_header(filepath, template, "c_block", False,
                                  current_year, contributors, dry_run):
                updated_count += 1
    
    # Update Python files
    if files_by_style["python"]:
        print(f"🐍 Updating {len(files_by_style['python'])} Python file(s)...")
        for filepath in files_by_style["python"]:
            if update_file_header(filepath, template, "python", True,
                                  current_year, contributors, dry_run):
                updated_count += 1
    
    action = "Would update" if dry_run else "Updated"
    print(f"\n✅ {action} {updated_count} file(s).")
    
    return updated_count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update copyright headers in target C/C++ and Python files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Target directories: {', '.join(TARGET_DIRS)}
Target files: {', '.join(TARGET_FILES)}
Excluded files: {', '.join(EXCLUDE_FILES)}

Examples:
  %(prog)s                          Update all target files
  %(prog)s --dry-run                Preview changes without modifying
  %(prog)s --list                   List target files
  %(prog)s include/options.hpp      Update specific file
  %(prog)s setup.py scripts/run.py  Update specific Python files
  %(prog)s -c "John <john@x.com>"   Add contributor info
  %(prog)s -t custom.tmpl           Use custom template
        """
    )
    
    parser.add_argument("files", nargs="*", help="Specific files to process (optional)")
    parser.add_argument("-t", "--template", default=DEFAULT_TEMPLATE, 
                        help=f"Template file path (default: {DEFAULT_TEMPLATE})")
    parser.add_argument("-c", "--contributor", action="append", dest="contributors",
                        help="Contributor entry (e.g., 'John Doe <john@example.com>')")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Show changes without writing")
    parser.add_argument("--list", action="store_true", dest="list_only",
                        help="List files that would be updated")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    
    args = parser.parse_args()
    
    try:
        count = update_copyright_header(
            specific_files=args.files if args.files else None,
            template_path=args.template,
            contributors=args.contributors,
            dry_run=args.dry_run,
            list_only=args.list_only
        )
        return 0 if count >= 0 else 1
    except FileNotFoundError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"✗ Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())